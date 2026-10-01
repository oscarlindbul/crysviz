// Unit tests for the NCI WASM backend (docs/compiled/nci_backend.c) and its JS
// adapter (docs/math/nci-backend-wasm.js).
//
// Run with `node --test tools/unittest/`. The module is built with the node
// environment, so it loads straight from docs/compiled (rebuild it with
// `cd docs/compiled && make wasm-nci` after changing the C).
//
// Everything below is set up in bohr and converted to Angstrom only where the
// request API wants it (the voxel vectors), so the analytic references can be
// written in atomic units throughout.

import { test } from 'node:test';
import assert from 'node:assert/strict';

import {
  ANGSTROM_PER_BOHR,
  NCI_DEFAULTS,
  computeNci,
  getNciBackend,
  promolecularDensityWithModule,
  sym3EigvalsWithModule,
} from '../../docs/math/nci-backend-wasm.js';
import { parseCube } from '../../docs/io/cubeParse.js';
import { boxedCubeLayout } from '../../docs/io/cubeLayout.js';
import { Field } from '../../docs/model/Field.js';
import { createNciFields } from '../../docs/model/NciField.js';
import { cartToFractional, invert3x3, transpose3x3 } from '../../docs/math/index.js';

const C_S = 1 / (2 * Math.cbrt(3 * Math.PI * Math.PI));
const S_CAP = 2.0;

/* ------------------------------------------------------------------ *
 * Reference helpers
 * ------------------------------------------------------------------ */

/** Symmetric 3x3 eigenvalues (ascending) by cyclic Jacobi rotations. */
function jacobiEigvals(m) {
  const a = m.map((row) => row.slice());
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) off += a[p][q] * a[p][q];
    if (off < 1e-300) break;
    for (let p = 0; p < 3; p++) {
      for (let q = p + 1; q < 3; q++) {
        if (a[p][q] === 0) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;
        for (let k = 0; k < 3; k++) {
          const akp = a[k][p];
          const akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < 3; k++) {
          const apk = a[p][k];
          const aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
    }
  }
  return [a[0][0], a[1][1], a[2][2]].sort((x, y) => x - y);
}

/** Grid point (i,j,k) in Cartesian bohr for voxel rows v (bohr). */
function gridPoint(v, i, j, k) {
  return [0, 1, 2].map((c) => i * v[0][c] + j * v[1][c] + k * v[2][c]);
}

const toAngstrom = (v) => v.map((row) => row.map((x) => x * ANGSTROM_PER_BOHR));

/** Sum of Gaussians A exp(-a |r - c|^2): value, gradient, Hessian. */
function gaussians(list, r) {
  let rho = 0;
  const g = [0, 0, 0];
  const H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const { A, a, c } of list) {
    const d = [r[0] - c[0], r[1] - c[1], r[2] - c[2]];
    const f = A * Math.exp(-a * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
    rho += f;
    for (let p = 0; p < 3; p++) {
      g[p] += -2 * a * d[p] * f;
      for (let q = 0; q < 3; q++) H[p][q] += f * (4 * a * a * d[p] * d[q] - (p === q ? 2 * a : 0));
    }
  }
  return { rho, g, H };
}

/** Sample an analytic density on a grid (x fastest), as float32. */
function sample(nx, ny, nz, v, fn) {
  const values = new Float32Array(nx * ny * nz);
  for (let k = 0; k < nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) values[i + nx * (j + ny * k)] = fn(gridPoint(v, i, j, k)).rho;
    }
  }
  return values;
}

const sOf = (rho, g) => C_S * Math.hypot(g[0], g[1], g[2]) * Math.pow(rho, -4 / 3);

/**
 * Compare an SCF result against an analytic density at every included point.
 * Returns the number of points checked and the sign checks performed.
 */
function checkAgainstAnalytic(result, {
  nx, ny, nz, v, fn, periodic, sTol, sAbs = 0.05, rhoMin = 1e-5, rhoPlot = 0.05,
}) {
  let checked = 0;
  let signChecked = 0;
  let maxErr = 0;
  for (let k = 0; k < nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const idx = i + nx * (j + ny * k);
        const boundary = !periodic && (i === 0 || j === 0 || k === 0
          || i === nx - 1 || j === ny - 1 || k === nz - 1);
        const { rho, g, H } = fn(gridPoint(v, i, j, k));
        const rho32 = Math.fround(rho);
        const inWindow = rho32 >= rhoMin && rho32 <= rhoPlot;
        if (boundary || !inWindow) {
          assert.equal(result.s[idx], S_CAP, `excluded point ${i},${j},${k} should carry sCap`);
          assert.equal(result.sl2rho[idx], 0, `excluded point ${i},${j},${k} should carry 0`);
          continue;
        }
        const sExact = Math.min(sOf(rho, g), S_CAP);
        const err = Math.abs(result.s[idx] - sExact);
        maxErr = Math.max(maxErr, err / (sExact + sAbs));
        assert.ok(err <= sTol * (sExact + sAbs),
          `s at ${i},${j},${k}: got ${result.s[idx]}, expected ${sExact}`);
        assert.ok(Math.abs(Math.abs(result.sl2rho[idx]) - rho) <= 1e-6 * rho + 1e-12,
          `|sign(l2) rho| at ${i},${j},${k}: got ${result.sl2rho[idx]}, expected ${rho}`);
        const ev = jacobiEigvals(H);
        const scale = Math.max(Math.abs(ev[0]), Math.abs(ev[2]));
        // Only where lambda2 is clearly away from zero; the FD Hessian is only
        // good to O(h^2) and the sign of a vanishing lambda2 is not defined.
        if (Math.abs(ev[1]) > 0.05 * scale) {
          assert.equal(Math.sign(result.sl2rho[idx]), Math.sign(ev[1]),
            `sign(l2) at ${i},${j},${k}: got ${result.sl2rho[idx]}, lambda = ${ev}`);
          signChecked++;
        }
        checked++;
      }
    }
  }
  return { checked, signChecked, maxErr };
}

/* ------------------------------------------------------------------ *
 * Eigen-solver
 * ------------------------------------------------------------------ */

test('closed-form symmetric eigen-solver matches Jacobi, including degenerate cases', async () => {
  const module = await getNciBackend();
  const cases = [
    [1, 2, 3, 0, 0, 0],
    [2, 2, 2, 0, 0, 0],
    [2, 2, 2, 1e-12, 0, 0],
    [1, 1, 1, 1, 1, 1],                 // eigenvalues 0, 0, 3
    [-0.1, -0.1, 0.5, 0, 0, 0],         // degenerate pair below a positive one
    [3, 3, 3, -1, -1, -1],              // 1, 4, 4 (degenerate pair on top)
    [1e-3, 2e-3, -4e-3, 5e-4, -2e-4, 7e-4],
    [4.2, -1.3, 0.7, 2.2, -0.9, 1.4],
  ];
  // Plus rotated diag(l, l, l + eps) matrices: near-degenerate, off-diagonal.
  const rot = (t) => {
    const [c1, s1, c2, s2] = [Math.cos(t), Math.sin(t), Math.cos(2 * t), Math.sin(2 * t)];
    const Rz = [[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]];
    const Rx = [[1, 0, 0], [0, c2, -s2], [0, s2, c2]];
    return Rz.map((row) => [0, 1, 2].map((c) => row.reduce((acc, x, k) => acc + x * Rx[k][c], 0)));
  };
  for (const [l, eps] of [[-0.02, 1e-9], [0.3, 1e-7], [-1, 0], [5e-4, 1e-10]]) {
    const R = rot(0.7);
    const D = [l, l + eps, -2 * l];
    const M = [0, 1, 2].map((p) => [0, 1, 2].map((q) => R[p].reduce((acc, x, k) => acc + x * D[k] * R[q][k], 0)));
    cases.push([M[0][0], M[1][1], M[2][2], M[0][1], M[0][2], M[1][2]]);
  }
  for (const a6 of cases) {
    const got = sym3EigvalsWithModule(module, a6);
    const ref = jacobiEigvals([[a6[0], a6[3], a6[4]], [a6[3], a6[1], a6[5]], [a6[4], a6[5], a6[2]]]);
    const scale = Math.max(...ref.map(Math.abs), 1e-300);
    for (let n = 0; n < 3; n++) {
      assert.ok(Math.abs(got[n] - ref[n]) <= 1e-7 * scale, `eig ${a6}: got ${got}, ref ${ref}`);
    }
    assert.ok(got[0] <= got[1] && got[1] <= got[2], `eigenvalues not sorted: ${got}`);
  }
});

/* ------------------------------------------------------------------ *
 * 1. Analytic Gaussian on an orthogonal grid
 * ------------------------------------------------------------------ */

test('SCF: single Gaussian on an orthogonal grid matches the closed form', async () => {
  const h = 0.2;
  const n = 41;
  const v = [[h, 0, 0], [0, h, 0], [0, 0, h]];
  const A = 0.04;
  const a = 0.5;
  const centre = [20.3 * h, 19.7 * h, 20.1 * h];
  const fn = (r) => gaussians([{ A, a, c: centre }], r);
  const values = sample(n, n, n, v, fn);

  const result = await computeNci({
    kind: 'scf', nx: n, ny: n, nz: n, voxel: toAngstrom(v), periodic: false, values,
  });

  // Closed form: s = c * 2 a r rho^(-1/3); lambda2 = -2 a rho < 0 everywhere,
  // so sign(lambda2) rho = -rho at every included point.
  for (let idx = 0; idx < values.length; idx++) {
    if (result.sl2rho[idx] !== 0) assert.ok(result.sl2rho[idx] < 0, 'single Gaussian: lambda2 < 0');
  }
  const { checked, maxErr } = checkAgainstAnalytic(result, {
    nx: n, ny: n, nz: n, v, fn, periodic: false, sTol: 0.03,
  });
  assert.equal(result.included, checked);
  assert.ok(checked > 10000, `too few points checked (${checked})`);
  assert.ok(maxErr < 0.03, `max relative s error ${maxErr}`);
});

/* ------------------------------------------------------------------ *
 * 2. Sheared grid (validates grad = V^-1 g, H = V^-1 H_idx V^-T)
 * ------------------------------------------------------------------ */

test('SCF: sheared grid gives the closed-form s and sign(lambda2) at the same Cartesian points', async () => {
  const h = 0.18;
  const v = [[h, 0, 0], [0.45 * h, 0.9 * h, 0], [0.3 * h, -0.25 * h, 0.85 * h]];
  const nx = 56;
  const ny = 44;
  const nz = 44;
  const mid = gridPoint(v, nx / 2, ny / 2, nz / 2);

  // Single Gaussian: s against the closed form.
  const single = (r) => gaussians([{ A: 0.04, a: 0.5, c: mid }], r);
  let result = await computeNci({
    kind: 'scf', nx, ny, nz, voxel: toAngstrom(v), periodic: false,
    values: sample(nx, ny, nz, v, single),
  });
  let stats = checkAgainstAnalytic(result, { nx, ny, nz, v, fn: single, periodic: false, sTol: 0.03 });
  assert.ok(stats.checked > 10000);

  // Three Gaussians in a triangle: regions with lambda2 > 0 (in-plane ring
  // region) and lambda2 < 0, checked against the analytic Hessian.
  const tri = [0, 1, 2].map((m) => {
    const t = (2 * Math.PI * m) / 3;
    return { A: 0.03, a: 0.4, c: [mid[0] + 2.0 * Math.cos(t), mid[1] + 2.0 * Math.sin(t), mid[2]] };
  });
  const three = (r) => gaussians(tri, r);
  result = await computeNci({
    kind: 'scf', nx, ny, nz, voxel: toAngstrom(v), periodic: false,
    values: sample(nx, ny, nz, v, three),
  });
  // Near the ring centre the gradient is a small difference of large terms,
  // so the O(h^2) error there is absolute rather than relative to s (a larger
  // sAbs); the point of this case is the sign of lambda2.
  stats = checkAgainstAnalytic(result, { nx, ny, nz, v, fn: three, periodic: false, sTol: 0.03, sAbs: 0.3 });
  let positive = 0;
  for (const x of result.sl2rho) if (x > 0) positive++;
  assert.ok(positive > 50, `expected lambda2 > 0 points near the ring centre, got ${positive}`);
  assert.ok(stats.signChecked > 10000);

  // The orthogonal-only treatment (per-axis step length, as Jmol does) must
  // NOT reproduce these numbers — guards against the transform being skipped.
  const naive = await computeNci({
    kind: 'scf', nx, ny, nz, periodic: false, values: sample(nx, ny, nz, v, single),
    voxel: toAngstrom([[h, 0, 0], [0, Math.hypot(...v[1]), 0], [0, 0, Math.hypot(...v[2])]]),
  });
  assert.throws(() => checkAgainstAnalytic(naive, { nx, ny, nz, v, fn: single, periodic: false, sTol: 0.03 }));
});

/* ------------------------------------------------------------------ *
 * 3. Periodic wrap
 * ------------------------------------------------------------------ */

test('SCF: periodic grid wraps neighbours, boundary points are finite and correct', async () => {
  const h = 0.16;
  const v = [[h, 0, 0], [0.3 * h, 0.95 * h, 0], [-0.2 * h, 0.1 * h, h]];
  const nx = 40;
  const ny = 36;
  const nz = 44;
  const n = [nx, ny, nz];
  // Reciprocal vectors of the cell (rows n_i * v_i) times 2 pi.
  const L = v.map((row, i) => row.map((x) => x * n[i]));
  const cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
  const vol = L[0].reduce((acc, x, c) => acc + x * cross(L[1], L[2])[c], 0);
  const b = [cross(L[1], L[2]), cross(L[2], L[0]), cross(L[0], L[1])].map((row) => row.map((x) => (2 * Math.PI * x) / vol));
  const terms = [
    { m: [1, 0, 0], c: 0.006, phase: 0.3 },
    { m: [0, 1, 1], c: 0.004, phase: 1.1 },
    { m: [1, -1, 1], c: 0.003, phase: -0.4 },
  ].map((t) => ({ ...t, G: [0, 1, 2].map((c) => t.m[0] * b[0][c] + t.m[1] * b[1][c] + t.m[2] * b[2][c]) }));
  const fn = (r) => {
    let rho = 0.02;
    const g = [0, 0, 0];
    const H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const { G, c, phase } of terms) {
      const x = G[0] * r[0] + G[1] * r[1] + G[2] * r[2] + phase;
      rho += c * Math.cos(x);
      for (let p = 0; p < 3; p++) {
        g[p] -= c * G[p] * Math.sin(x);
        for (let q = 0; q < 3; q++) H[p][q] -= c * G[p] * G[q] * Math.cos(x);
      }
    }
    return { rho, g, H };
  };

  const result = await computeNci({
    kind: 'scf', nx, ny, nz, voxel: toAngstrom(v), periodic: true, values: sample(nx, ny, nz, v, fn),
  });
  assert.equal(result.included, nx * ny * nz, 'every point of a periodic grid in the window is included');
  const stats = checkAgainstAnalytic(result, { nx, ny, nz, v, fn, periodic: true, sTol: 0.03 });
  assert.equal(stats.checked, nx * ny * nz);
  assert.ok(stats.signChecked > 0.5 * nx * ny * nz);
  for (const x of result.s) assert.ok(Number.isFinite(x));
});

test('SCF: non-periodic boundary layer is excluded, valueScale and abs are applied', async () => {
  const n = 8;
  const values = new Float32Array(n * n * n).fill(-0.2);  // |-0.2 * 0.1| = 0.02
  const result = await computeNci({
    kind: 'scf', nx: n, ny: n, nz: n, voxel: [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]],
    periodic: false, values, valueScale: 0.1,
  });
  assert.equal(result.included, (n - 2) ** 3);
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const idx = i + n * (j + n * k);
        const edge = [i, j, k].some((x) => x === 0 || x === n - 1);
        assert.equal(result.s[idx], edge ? S_CAP : 0);
        // Constant field: Hessian 0, lambda2 = 0 is not < 0, so +rho.
        assert.equal(result.sl2rho[idx], edge ? 0 : Math.fround(0.02));
      }
    }
  }
  // The caller's buffer is untouched.
  assert.equal(values[0], Math.fround(-0.2));
});

/* ------------------------------------------------------------------ *
 * 4. Promolecular
 * ------------------------------------------------------------------ */

// NCIPLOT/Jmol tables, restated independently of the C for the reference.
const PRO = {
  1: { c: [0.2815], z: [0.5288], dMax: 2.982502423 },
  6: { c: [120.2, 1.172], z: [0.0884, 0.5480], dMax: 3.872424373 },
  8: { c: [289.5, 2.879], z: [0.0669, 0.3974], dMax: 3.165369971 },
};

/** Promolecular rho, gradient and Hessian with Jmol's dMax box test. */
function promolecular(atoms, r) {
  let rho = 0;
  const g = [0, 0, 0];
  const H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const { Z, pos } of atoms) {
    const p = PRO[Z];
    const d = [r[0] - pos[0], r[1] - pos[1], r[2] - pos[2]];
    if (d.some((x) => Math.abs(x) > p.dMax)) continue;
    const dist = Math.hypot(...d);
    let f0 = 0;
    let f1 = 0;
    let f2 = 0;
    p.c.forEach((c, n) => {
      const e = c * Math.exp(-dist / p.z[n]);
      f0 += e;
      f1 -= e / p.z[n];            // d rho / dr
      f2 += e / (p.z[n] * p.z[n]); // d2 rho / dr2
    });
    rho += f0;
    if (dist === 0) continue;
    const u = d.map((x) => x / dist);
    for (let a = 0; a < 3; a++) {
      g[a] += f1 * u[a];
      for (let b = 0; b < 3; b++) {
        H[a][b] += f2 * u[a] * u[b] + (f1 / dist) * ((a === b ? 1 : 0) - u[a] * u[b]);
      }
    }
  }
  return { rho, g, H };
}

/** Request for atoms given in Cartesian bohr on an orthogonal grid of step h. */
function promolecularRequest(atoms, n, h, extra = {}) {
  return {
    kind: 'promolecular', nx: n[0], ny: n[1], nz: n[2], periodic: false,
    voxel: toAngstrom([[h, 0, 0], [0, h, 0], [0, 0, h]]),
    atoms: { Z: atoms.map((a) => a.Z), frac: atoms.map((a) => a.pos.map((x, c) => x / (n[c] * h))) },
    ...extra,
  };
}

test('promolecular: single H atom matches the analytic density and s', async () => {
  const module = await getNciBackend();
  const h = 0.1;
  const n = [61, 61, 61];
  const atoms = [{ Z: 1, pos: [30.25 * h, 29.6 * h, 30.1 * h] }];
  const req = promolecularRequest(atoms, n, h);
  const rho = promolecularDensityWithModule(module, req);
  const result = await computeNci(req);
  const v = [[h, 0, 0], [0, h, 0], [0, 0, h]];

  let included = 0;
  for (let k = 0; k < n[2]; k++) {
    for (let j = 0; j < n[1]; j++) {
      for (let i = 0; i < n[0]; i++) {
        const idx = i + n[0] * (j + n[1] * k);
        const ref = promolecular(atoms, gridPoint(v, i, j, k));
        assert.ok(Math.abs(rho[idx] - ref.rho) <= 1e-6 * ref.rho + 1e-12, `rho at ${i},${j},${k}`);
        if (ref.rho < 1e-5 || ref.rho > NCI_DEFAULTS.rhoPlotPromolecular) {
          assert.equal(result.s[idx], S_CAP);
          assert.equal(result.sl2rho[idx], 0);
          continue;
        }
        included++;
        // Single exponential: |rho'| = rho / zeta, so s = c / (zeta rho^(1/3)).
        const sExact = Math.min(C_S / (0.5288 * Math.cbrt(ref.rho)), S_CAP);
        assert.ok(Math.abs(result.s[idx] - sExact) <= 1e-5 * sExact, `s at ${i},${j},${k}`);
        // Radial curvature rho/zeta^2 > 0, two tangential -rho/(zeta r) < 0:
        // lambda2 is tangential and negative.
        assert.ok(Math.abs(result.sl2rho[idx] + ref.rho) <= 1e-6 * ref.rho, `sl2rho at ${i},${j},${k}`);
      }
    }
  }
  assert.equal(result.included, included);
  assert.ok(included > 10000);
  assert.deepEqual(result.clampedElements, []);
});

test('promolecular: H2 bond midpoint has lambda2 < 0, H3 ring centre has lambda2 > 0', async () => {
  const h = 0.1;
  const n = [81, 81, 61];
  const centre = [40 * h, 40 * h, 30 * h];  // a grid point
  const v = [[h, 0, 0], [0, h, 0], [0, 0, h]];
  const at = (idx) => idx[0] + n[0] * (idx[1] + n[1] * idx[2]);
  const rhoH = (r) => 0.2815 * Math.exp(-r / 0.5288);

  // H2, d = 3 bohr, midpoint at the centre. Each atom adds
  // (rho'' - rho'/r) u u^T + (rho'/r) I, so on the axis the curvature is
  // 2 rho_a/zeta^2 > 0 and perpendicular 2 rho_a'/r = -2 rho_a/(zeta r) < 0
  // (twice): lambda2 is the degenerate perpendicular pair, negative.
  let atoms = [
    { Z: 1, pos: [centre[0] - 1.5, centre[1], centre[2]] },
    { Z: 1, pos: [centre[0] + 1.5, centre[1], centre[2]] },
  ];
  let result = await computeNci(promolecularRequest(atoms, n, h));
  let expected = 2 * rhoH(1.5);
  assert.ok(Math.abs(result.sl2rho[at([40, 40, 30])] + expected) <= 1e-6 * expected);
  assert.ok(result.s[at([40, 40, 30])] < 1e-6, 'gradient vanishes at the midpoint');
  const ref = promolecular(atoms, gridPoint(v, 40, 40, 30));
  assert.ok(jacobiEigvals(ref.H)[1] < 0);

  // Jmol's in-loop "rho < rhoMin" exit would drop this point with rhoMin
  // between one atom's and both atoms' density; the sum must be used.
  result = await computeNci(promolecularRequest(atoms, n, h, { params: { rhoMin: 1.2 * rhoH(1.5) } }));
  assert.ok(Math.abs(result.sl2rho[at([40, 40, 30])] + expected) <= 1e-6 * expected);

  // Equilateral H3, centroid at distance R = 2 bohr from each atom, in the
  // xy plane. In-plane: sum u u^T = 3/2 I, so the in-plane eigenvalue is
  // 3/2 (rho'' + rho'/r) = 3/2 rho_a/zeta^2 (1 - zeta/R) > 0 (R > zeta, twice);
  // out of plane 3 rho'/R < 0. lambda2 is in-plane, positive.
  const R = 2;
  atoms = [0, 1, 2].map((m) => {
    const t = (2 * Math.PI * m) / 3 + 0.2;
    return { Z: 1, pos: [centre[0] + R * Math.cos(t), centre[1] + R * Math.sin(t), centre[2]] };
  });
  result = await computeNci(promolecularRequest(atoms, n, h));
  expected = 3 * rhoH(R);
  assert.ok(Math.abs(result.sl2rho[at([40, 40, 30])] - expected) <= 1e-6 * expected,
    `ring centre: got ${result.sl2rho[at([40, 40, 30])]}, expected +${expected}`);
});

test('promolecular: Z > 18 is clamped to Ar and reported, Z < 1 dropped', async () => {
  const h = 0.2;
  const n = [30, 30, 30];
  const base = promolecularRequest([{ Z: 18, pos: [3, 3, 3] }], n, h);
  const clamped = { ...base, atoms: { Z: [26, 0], frac: [base.atoms.frac[0], [0.5, 0.5, 0.5]] } };
  const a = await computeNci(base);
  const b = await computeNci(clamped);
  assert.deepEqual(b.clampedElements, ['Fe']);
  assert.deepEqual(Array.from(b.s), Array.from(a.s));
  assert.deepEqual(Array.from(b.sl2rho), Array.from(a.sl2rho));
});

test('promolecular: periodic images give a translation-invariant field and the brute-force density', async () => {
  const module = await getNciBackend();
  const h = 0.25;
  const n = [20, 24, 18];  // a 5 x 6 x 4.5 bohr cell: images matter on every side
  const v = [[h, 0, 0], [0.1 * h, h, 0], [0, 0.2 * h, h]];
  const frac = [[0.1, 0.95, 0.5], [0.6, 0.4, 0.05]];
  const Z = [8, 1];
  const req = {
    kind: 'promolecular', nx: n[0], ny: n[1], nz: n[2], periodic: true, voxel: toAngstrom(v),
    atoms: { Z, frac },
  };
  const rho = promolecularDensityWithModule(module, req);

  // Brute force: every image in a generous shell, Jmol's box test.
  const L = v.map((row, i) => row.map((x) => x * n[i]));
  const images = [];
  for (let a = 0; a < Z.length; a++) {
    for (let t0 = -3; t0 <= 3; t0++) {
      for (let t1 = -3; t1 <= 3; t1++) {
        for (let t2 = -3; t2 <= 3; t2++) {
          const f = [frac[a][0] + t0, frac[a][1] + t1, frac[a][2] + t2];
          images.push({ Z: Z[a], pos: [0, 1, 2].map((c) => f[0] * L[0][c] + f[1] * L[1][c] + f[2] * L[2][c]) });
        }
      }
    }
  }
  for (let k = 0; k < n[2]; k++) {
    for (let j = 0; j < n[1]; j++) {
      for (let i = 0; i < n[0]; i++) {
        const ref = promolecular(images, gridPoint(v, i, j, k)).rho;
        const got = rho[i + n[0] * (j + n[1] * k)];
        assert.ok(Math.abs(got - ref) <= 1e-6 * ref + 1e-12, `rho at ${i},${j},${k}: ${got} vs ${ref}`);
      }
    }
  }

  // Shift every atom by one grid step along each axis: the field must roll.
  const base = await computeNci(req);
  const shift = [3, 5, 2];
  const moved = await computeNci({
    ...req, atoms: { Z, frac: frac.map((f) => f.map((x, c) => x + shift[c] / n[c])) },
  });
  let maxDiff = 0;
  for (let k = 0; k < n[2]; k++) {
    for (let j = 0; j < n[1]; j++) {
      for (let i = 0; i < n[0]; i++) {
        const src = i + n[0] * (j + n[1] * k);
        const dst = ((i + shift[0]) % n[0]) + n[0] * (((j + shift[1]) % n[1]) + n[1] * ((k + shift[2]) % n[2]));
        maxDiff = Math.max(maxDiff, Math.abs(base.s[src] - moved.s[dst]),
          Math.abs(base.sl2rho[src] - moved.sl2rho[dst]));
      }
    }
  }
  assert.ok(maxDiff < 1e-4, `translation invariance: max diff ${maxDiff}`);
  assert.ok(base.included > 100);
});

/* ------------------------------------------------------------------ *
 * 5. Promolecular vs SCF on the promolecular density
 * ------------------------------------------------------------------ */

test('SCF on the sampled promolecular density agrees with the analytic promolecular NCI', async () => {
  const module = await getNciBackend();
  const h = 0.08;
  const n = [120, 90, 80];
  const atoms = [
    { Z: 1, pos: [2.6, 3.6, 3.2] },
    { Z: 1, pos: [5.4, 3.6, 3.2] },
    { Z: 6, pos: [6.0, 3.9, 3.4] },
  ];
  const req = promolecularRequest(atoms, n, h, { params: { rhoPlot: 0.05 } });
  const pro = await computeNci(req);
  const rho = promolecularDensityWithModule(module, req);
  const scf = await computeNci({
    kind: 'scf', nx: n[0], ny: n[1], nz: n[2], voxel: req.voxel, periodic: false, values: rho,
  });

  // Compare only where the finite-difference stencil does not straddle one of
  // the dMax box faces (the truncated density jumps there) and stays off the
  // nuclei.
  const v = [[h, 0, 0], [0, h, 0], [0, 0, h]];
  let compared = 0;
  let signCompared = 0;
  let maxRel = 0;
  for (let k = 1; k < n[2] - 1; k++) {
    for (let j = 1; j < n[1] - 1; j++) {
      for (let i = 1; i < n[0] - 1; i++) {
        const idx = i + n[0] * (j + n[1] * k);
        if (pro.sl2rho[idx] === 0 || scf.sl2rho[idx] === 0) continue;
        const r = gridPoint(v, i, j, k);
        const nearFace = atoms.some(({ Z, pos }) => r.some((x, c) => Math.abs(Math.abs(x - pos[c]) - PRO[Z].dMax) < 1.01 * h));
        if (nearFace) continue;
        compared++;
        const rel = Math.abs(scf.s[idx] - pro.s[idx]) / (pro.s[idx] + 0.05);
        maxRel = Math.max(maxRel, rel);
        const ev = jacobiEigvals(promolecular(atoms, r).H);
        if (Math.abs(ev[1]) > 0.05 * Math.max(Math.abs(ev[0]), Math.abs(ev[2]))) {
          assert.equal(Math.sign(scf.sl2rho[idx]), Math.sign(pro.sl2rho[idx]), `sign at ${i},${j},${k}`);
          signCompared++;
        }
      }
    }
  }
  assert.ok(compared > 50000, `compared ${compared}`);
  assert.ok(signCompared > 0.8 * compared);
  assert.ok(maxRel < 0.03, `max relative s difference ${maxRel}`);
});

test('computeNci rejects malformed requests', async () => {
  const voxel = [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]];
  await assert.rejects(computeNci({ kind: 'scf', nx: 2, ny: 2, nz: 2, voxel, periodic: false,
    values: new Float32Array(7) }), /needs 8 values/);
  await assert.rejects(computeNci({ kind: 'bogus', nx: 2, ny: 2, nz: 2, voxel, periodic: false }), /unknown kind/);
  await assert.rejects(computeNci({ kind: 'scf', nx: 2, ny: 2, nz: 2, periodic: false,
    voxel: [[1, 0, 0], [2, 0, 0], [0, 0, 1]], values: new Float32Array(8) }), /singular/);
  await assert.rejects(computeNci({ kind: 'promolecular', nx: 2, ny: 2, nz: 2, voxel, periodic: false }), /atoms/);
});

/* ------------------------------------------------------------------ *
 * 7. Cube files: the Field.periodic flag decides the wrap
 *
 * A cube loads periodic by default and as a finite block only when the user
 * says so (io/ReadCubeModule.js buildCubeStructure). NCI follows the flag:
 * a periodic grid wraps its neighbours, a block drops its boundary layer and
 * gets no periodic atom images. The fields below are built from parseCube's
 * output exactly as buildCubeStructure builds them (that module pulls in the
 * renderer and does not load under node; the browser test cubereader.test.js
 * covers it end to end).
 * ------------------------------------------------------------------ */

const CUBE_N = 12;
const CUBE_STEP = 0.5; // bohr
const CUBE_ORIGIN = [-1, -2, -3]; // bohr

/** Periodic sum of a Gaussian density centred at fractional `centre`, as cube text in bohr. */
function periodicGaussianCube(centre, { atoms = [] } = {}) {
  const n = CUBE_N;
  const L = n * CUBE_STEP;
  const sigma = 1.0;
  const rho = (i, j, k) => {
    let sum = 0;
    const p = [i, j, k].map((x) => x * CUBE_STEP);
    for (let a = -1; a <= 1; a++) {
      for (let b = -1; b <= 1; b++) {
        for (let c = -1; c <= 1; c++) {
          const d = [p[0] - (centre[0] + a) * L, p[1] - (centre[1] + b) * L, p[2] - (centre[2] + c) * L];
          sum += 0.03 * Math.exp(-(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) / (2 * sigma * sigma));
        }
      }
    }
    return sum;
  };
  const lines = ['NCI wrap test', 'Electron density from Total SCF Density',
    `${atoms.length} ${CUBE_ORIGIN.join(' ')}`,
    `${n} ${CUBE_STEP} 0 0`, `${n} 0 ${CUBE_STEP} 0`, `${n} 0 0 ${CUBE_STEP}`,
    ...atoms.map(({ Z, pos }) => `${Z} ${Z}.0 ${pos.join(' ')}`)];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const rec = [];
      for (let k = 0; k < n; k++) rec.push(rho(i, j, k).toExponential(8));
      for (let s = 0; s < rec.length; s += 6) lines.push(rec.slice(s, s + 6).join(' '));
    }
  }
  return lines.join('\n');
}

/** The Field buildCubeStructure makes for dataset 0, with the given periodicity and origin. */
function cubeField(cube, periodic, origin = [0, 0, 0]) {
  return new Field({
    nx: cube.grid[0], ny: cube.grid[1], nz: cube.grid[2], origin: [...origin],
    voxel: cube.voxel.map((row) => [...row]), values: cube.values[0], component: 0,
    label: cube.datasetLabels[0], valueUnit: cube.datasetUnits[0], periodic,
  });
}

const onBoundary = (i, j, k, n) => [i, j, k].some((x) => x === 0 || x === n - 1);

test('NCI on a periodic cube wraps its neighbours; on a block cube it drops the boundary layer', async () => {
  const n = CUBE_N;
  const centre = [0.02, 0.5, 0.5]; // straddles the x = 0 face, so the wrap matters there
  const cube = parseCube(periodicGaussianCube(centre));
  assert.equal(cube.datasetUnits[0], 'e/bohr³');

  const periodic = await createNciFields(cubeField(cube, true), { kind: 'scf' });
  const block = await createNciFields(cubeField(cube, false), { kind: 'scf' });
  assert.equal(periodic.sField.periodic, true);
  assert.equal(periodic.colourField.periodic, true);
  assert.equal(block.sField.periodic, false);
  assert.equal(block.colourField.periodic, false);

  // Wrapped: rolling the density by whole grid steps rolls the NCI result,
  // boundary points included, which only holds when neighbours wrap.
  const shift = [5, 3, 7];
  const rolled = parseCube(periodicGaussianCube(centre.map((c, a) => c + shift[a] / n)));
  const moved = await createNciFields(cubeField(rolled, true), { kind: 'scf' });
  let maxDiff = 0;
  let boundaryIncluded = 0;
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const src = i + n * (j + n * k);
        const dst = ((i + shift[0]) % n) + n * (((j + shift[1]) % n) + n * ((k + shift[2]) % n));
        maxDiff = Math.max(maxDiff, Math.abs(periodic.sField.values[src] - moved.sField.values[dst]),
          Math.abs(periodic.colourField.values[src] - moved.colourField.values[dst]));
        if (onBoundary(i, j, k, n) && periodic.sField.values[src] < S_CAP) boundaryIncluded++;
      }
    }
  }
  assert.ok(maxDiff < 1e-4, `periodic NCI is translation invariant: max diff ${maxDiff}`);
  assert.ok(boundaryIncluded > 0, 'the periodic boundary layer has real values');

  // Block: every boundary point is excluded (s = sCap, sign(λ2)ρ = 0), and
  // deep inside the grid, where no stencil reaches the edge, the result is the
  // periodic one.
  let interiorDiff = 0;
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const idx = i + n * (j + n * k);
        if (onBoundary(i, j, k, n)) {
          assert.equal(block.sField.values[idx], S_CAP, `block boundary s at ${i},${j},${k}`);
          assert.equal(block.colourField.values[idx], 0, `block boundary sign(λ2)ρ at ${i},${j},${k}`);
        } else if ([i, j, k].every((x) => x >= 2 && x <= n - 3)) {
          interiorDiff = Math.max(interiorDiff, Math.abs(block.sField.values[idx] - periodic.sField.values[idx]));
        }
      }
    }
  }
  assert.ok(interiorDiff < 1e-5, `block interior matches the periodic result: ${interiorDiff}`);
  assert.ok(block.info.included < periodic.info.included, 'the block drops points the periodic grid keeps');
});

test('promolecular NCI on a cube: periodic adds atom images, a block (padded box, shifted origin) does not', async () => {
  // One H atom near the x = 0 face of the grid, in bohr.
  const atomBohr = [CUBE_ORIGIN[0] + 0.3, CUBE_ORIGIN[1] + 3, CUBE_ORIGIN[2] + 3];
  const cube = parseCube(periodicGaussianCube([0.5, 0.5, 0.5], { atoms: [{ Z: 1, pos: atomBohr }] }));
  const atomFracInGrid = [0, 1, 2].map((a) => (atomBohr[a] - CUBE_ORIGIN[a]) / (CUBE_N * CUBE_STEP));
  const elements = ['H'];
  const fracIn = (lattice, cart) => cartToFractional(cart, lattice, invert3x3(transpose3x3(lattice)));

  // Periodic load: the cell is the grid box and atoms are shifted by the origin.
  const periodicStructure = {
    lattice: cube.lattice, elements,
    atoms: [{ position: fracIn(cube.lattice, cube.atoms[0].position.map((c, k) => c - cube.origin[k])) }],
  };
  const periodic = await createNciFields(cubeField(cube, true), { kind: 'promolecular', structure: periodicStructure });

  // Block load: a padded box around atoms and grid, the field at layout.fieldOrigin.
  const layout = boxedCubeLayout(cube);
  const blockStructure = { lattice: layout.lattice, elements, atoms: [{ position: fracIn(layout.lattice, layout.positions[0]) }] };
  const block = await createNciFields(cubeField(cube, false, layout.fieldOrigin),
    { kind: 'promolecular', structure: blockStructure });

  // References straight from the backend, atoms fractional in the grid's cell.
  const base = {
    kind: 'promolecular', nx: CUBE_N, ny: CUBE_N, nz: CUBE_N, voxel: cube.voxel,
    atoms: { Z: [1], frac: [atomFracInGrid] },
  };
  const refPeriodic = await computeNci({ ...base, periodic: true });
  const refBlock = await computeNci({ ...base, periodic: false });
  const maxDiff = (a, b) => a.reduce((m, x, i) => Math.max(m, Math.abs(x - b[i])), 0);
  assert.ok(maxDiff(periodic.sField.values, refPeriodic.s) < 1e-5, 'periodic cube = periodic request');
  assert.ok(maxDiff(block.sField.values, refBlock.s) < 1e-5,
    'block cube = non-periodic request with the atom placed relative to the field origin');
  assert.ok(maxDiff(block.colourField.values, refBlock.sl2rho) < 1e-7);

  // The image of the atom across the x = 0 face reaches the far side of a
  // periodic grid only: the density there differs between the two.
  const n = CUBE_N;
  let farSideDiffers = 0;
  for (let k = 2; k < n - 2; k++) {
    for (let j = 2; j < n - 2; j++) {
      const idx = (n - 2) + n * (j + n * k);
      if (Math.abs(periodic.colourField.values[idx] - block.colourField.values[idx]) > 1e-8) farSideDiffers++;
    }
  }
  assert.ok(farSideDiffers > 0, 'periodic images change the far side of the grid');
});
