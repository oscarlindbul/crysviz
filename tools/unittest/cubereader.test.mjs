// Gaussian cube parser (docs/io/cubeParse.js) against the cubegen spec, plus
// the pieces around it that load in node: the '.cub' name rule in
// docs/io/formats.js and the valueUnit / periodic propagation of Field and
// CompositeField. Run with `node --test tools/unittest/`.
// (Ported from the NCI branch's original cube parser tests: the same fixtures and
// expectations, read through parseCube's CubeData — `grid`, `units`,
// `datasetIds`, `datasetLabels`, `datasetUnits`, `values`.)
//
// Fixtures in fixtures/cube/ are hand-sized grids whose values encode their
// own file position (100·i + 10·j + k along axes 1, 2, 3), so a wrong index
// order shows up as a wrong number rather than as a plausible-looking field.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import {
  parseCube, cubeDensityUnit, BOHR_TO_ANGSTROM, CUBE_DENSITY_UNIT, GRADIENT_LABELS,
} from '../../docs/io/cubeParse.js';
import { detectFormat, headOf } from '../../docs/io/formats.js';
import { Field } from '../../docs/model/Field.js';
import { combineFields, magnitudeField } from '../../docs/model/CompositeField.js';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const fixture = (name) => fs.readFileSync(path.join(HERE, 'fixtures', 'cube', name), 'utf8');
const B = BOHR_TO_ANGSTROM;

/** Stored index of grid point (i, j, k): axis 1 fastest. */
const at = (cube, i, j, k) => i + cube.grid[0] * (j + cube.grid[1] * k);

function assertClose(actual, expected, message, tol = 1e-9) {
  assert.ok(Math.abs(actual - expected) <= tol, `${message}: ${actual} vs ${expected}`);
}
function assertVecClose(actual, expected, message, tol = 1e-9) {
  assert.equal(actual.length, expected.length, `${message}: length`);
  actual.forEach((x, i) => assertClose(x, expected[i], `${message}[${i}]`, tol));
}

test('(a) bohr cube, NVal = 1: units, atoms, labels and x-fastest ordering', () => {
  const cube = parseCube(fixture('bohr_density.cube'));
  // `label` joins both comment lines; the dataset label is line 2 alone.
  assert.equal(cube.label, 'Water fragment test cube Electron density from Total SCF Density');
  assert.equal(cube.units, 'bohr');
  assert.deepEqual(cube.grid, [2, 3, 4]);
  assert.equal(cube.values.length, 1);
  assert.equal(cube.datasetIds, null);

  assertVecClose(cube.origin, [-B, -B, -B], 'origin in Å');
  assertVecClose(cube.voxel[0], [0.5 * B, 0, 0], 'voxel a');
  assertVecClose(cube.lattice[0], [1.0 * B, 0, 0], 'lattice a = 2 × 0.5 bohr');
  assertVecClose(cube.lattice[1], [0, 1.5 * B, 0], 'lattice b = 3 × 0.5 bohr');
  assertVecClose(cube.lattice[2], [0, 0, 2.0 * B], 'lattice c = 4 × 0.5 bohr');

  assert.equal(cube.atoms.length, 2);
  assert.equal(cube.atoms[0].atomicNumber, 8);
  assert.equal(cube.atoms[1].atomicNumber, 1);
  assertVecClose(cube.atoms[0].position, [0, 0, 0], 'O position');
  // One unit for all three coordinates (the old reader used a per-axis flag).
  assertVecClose(cube.atoms[1].position, [-0.5 * B, -0.5 * B, 0.5 * B], 'H position in Å');

  assert.deepEqual(cube.datasetLabels, ['Electron density from Total SCF Density']);
  assert.deepEqual(cube.datasetUnits, [CUBE_DENSITY_UNIT]);
  const values = cube.values[0];
  assert.equal(values.length, 24);
  for (let i = 0; i < 2; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 4; k++) {
        assert.equal(values[at(cube, i, j, k)], 100 * i + 10 * j + k, `value at ${i},${j},${k}`);
      }
    }
  }
});

test('(b) Å cube (negative N): axes are not flipped, atoms stay in Å', () => {
  const cube = parseCube(fixture('angstrom_potential.cube'));
  assert.equal(cube.units, 'angstrom');
  assert.deepEqual(cube.grid, [2, 2, 3]);
  assertVecClose(cube.voxel[0], [0.5, 0, 0], 'voxel a');
  assertVecClose(cube.lattice[0], [1, 0, 0], 'lattice a');
  assertVecClose(cube.lattice[1], [0, 1, 0], 'lattice b');
  assertVecClose(cube.lattice[2], [0, 0, 1.5], 'lattice c');
  assertVecClose(cube.atoms[0].position, [0.25, 0.5, 0.75], 'atom position');
  // A potential is not a density, whatever its description ends with.
  assert.deepEqual(cube.datasetLabels, ['Electrostatic potential from Total SCF Density']);
  assert.deepEqual(cube.datasetUnits, [null]);
  assert.equal(cube.values[0][at(cube, 1, 1, 2)], 112);
  assert.equal(cube.values[0][at(cube, 0, 1, 2)], 12);
});

test('(c) NVal = 4 gradient cube: four de-interleaved fields', () => {
  const cube = parseCube(fixture('gradient_nval4.cube'));
  assert.equal(cube.values.length, 4);
  assert.deepEqual(cube.datasetLabels, ['Density', '∂ρ/∂x', '∂ρ/∂y', '∂ρ/∂z']);
  assert.deepEqual(cube.datasetLabels, GRADIENT_LABELS);
  assert.deepEqual(cube.datasetUnits, [CUBE_DENSITY_UNIT, null, null, null]);
  for (let v = 0; v < 4; v++) {
    for (let i = 0; i < 2; i++) {
      for (let j = 0; j < 2; j++) {
        for (let k = 0; k < 2; k++) {
          assert.equal(cube.values[v][at(cube, i, j, k)], 100 * i + 10 * j + k + 0.25 * v,
            `value ${v} at ${i},${j},${k}`);
        }
      }
    }
  }
});

test('(d) orbital cube (NAtoms < 0): MO record read, one field per orbital', () => {
  const cube = parseCube(fixture('mo_two_orbitals.cube'));
  assert.equal(cube.atoms.length, 1);
  assert.deepEqual(cube.datasetIds, [5, 6]);
  assert.equal(cube.values.length, 2);
  assert.deepEqual(cube.datasetLabels, ['MO 5', 'MO 6']);
  assert.deepEqual(cube.datasetUnits, [null, null]);
  for (let i = 0; i < 2; i++) {
    for (let k = 0; k < 3; k++) {
      assert.equal(cube.values[0][at(cube, i, 0, k)], 10 * i + k + 1);
      assert.equal(cube.values[1][at(cube, i, 0, k)], -(10 * i + k + 1));
    }
  }
});

test('(e) truncated data throws instead of leaving zeros', () => {
  assert.throws(() => parseCube(fixture('truncated.cube')), /truncated.*expected 24 data values.*found 20/);
  assert.throws(() => parseCube(' t\n d\n 1 0 0 0\n 2 1 0 0\n'), /axis 2/);
  assert.throws(() => parseCube(' t\n d\n 2 0 0 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n 1 1 0 0 0\n'),
    /atom line 2 of 2/);
  assert.throws(() => parseCube(' t\n d\n 0 0 0 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n abc\n'),
    /"abc" is not a number/);
});

test('(f) non-orthogonal axes, CRLF line ends and Fortran D exponents', () => {
  const cube = parseCube(fixture('sheared_spin.cub'));
  assert.ok(cube.label.startsWith('Sheared grid'), cube.label);
  assert.deepEqual(cube.grid, [3, 3, 2]);
  assertVecClose(cube.lattice[0], [1.5 * B, 0, 0], 'lattice a', 1e-6);
  assertVecClose(cube.lattice[1], [0.75 * B, 3 * 0.4330127 * B, 0], 'lattice b', 1e-6);
  assertVecClose(cube.lattice[2], [0.2 * B, 0.2 * B, 1.0 * B], 'lattice c', 1e-6);
  // The atom sits at origin + 1·v1 + 2·v2 (+0·v3) in bohr.
  const shifted = cube.atoms[0].position.map((c, i) => c - cube.origin[i]);
  const expected = [0, 1, 2].map((c) => cube.voxel[0][c] + 2 * cube.voxel[1][c]);
  assertVecClose(shifted, expected, 'atom offset from the origin', 1e-6);
  // Spin densities are densities.
  assert.deepEqual(cube.datasetUnits, [CUBE_DENSITY_UNIT]);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 2; k++) {
        assert.equal(cube.values[0][at(cube, i, j, k)], 100 * i + 10 * j + k);
      }
    }
  }
});

test('tokenizer: exponent forms, long mantissas and a blank title', () => {
  const text = [
    '', ' ', '    0    0.0    0.0    0.0', '    1  1 0 0', '    1  0 1 0', '    8  0 0 1',
    '  1.5E-3 -2.25e+2 3.0D-1 0.12345-100 1.234567890123456789 -7',
    '  .5 +4.',
  ].join('\n');
  const cube = parseCube(text);
  assert.deepEqual(cube.datasetLabels, ['Cube data']);
  assert.deepEqual(cube.datasetUnits, [null]);
  const expected = [1.5e-3, -225, 0.3, 0.12345e-100, 1.234567890123456789, -7, 0.5, 4];
  assert.deepEqual([...cube.values[0]], expected.map(Math.fround));
});

test('cubeDensityUnit: densities vs other cube quantities', () => {
  const cases = [
    ['Electron density from Total SCF Density', CUBE_DENSITY_UNIT],
    ['Spin density from Total SCF Density', CUBE_DENSITY_UNIT],
    ['Alpha density from Total SCF Density', CUBE_DENSITY_UNIT],
    ['Electrostatic potential from Total SCF Density', null],
    ['Laplacian of electron density from Total SCF Density', null],
    ['Norm of density gradient from Total SCF Density', null],
    ['MO coefficients', null],
    ['Contains the selected quantity on a FFT grid', null],
    ['', null],
  ];
  for (const [description, unit] of cases) assert.equal(cubeDensityUnit(description), unit, description);
});

test('formats.js: .cub and .cube are cube files by name, .cubic is not', () => {
  const byName = (name) => detectFormat({ fileName: name, head: headOf('') }).id;
  assert.equal(byName('NCI.cub'), 'cube');
  assert.equal(byName('density.cube'), 'cube');
  assert.equal(byName('x.CUB'), 'cube');
  assert.notEqual(byName('x.cubic'), 'cube');
  // An orbital cube (negative NAtoms) is still sniffed as a cube under any name.
  assert.equal(detectFormat({ fileName: 'orbitals.dat', head: headOf(fixture('mo_two_orbitals.cube')) }).id, 'cube');
});

test('Field: valueUnit / periodic / colorBy defaults and combined-field propagation', () => {
  const grid = (valueUnit, periodic) => new Field({
    nx: 1, ny: 1, nz: 2, values: new Float32Array([1, 2]), valueUnit, periodic, label: 'f',
  });
  const plain = new Field({});
  assert.equal(plain.valueUnit, null);
  assert.equal(plain.periodic, true);
  assert.equal(plain.colorBy, null);

  const a = grid('e/Å³', true);
  const b = grid('e/Å³', true);
  const sum = combineFields([{ field: a, weight: 0.5 }, { field: b, weight: -0.5 }]);
  assert.equal(sum.valueUnit, 'e/Å³');
  assert.equal(sum.periodic, true);
  assert.equal(magnitudeField([a, b]).valueUnit, 'e/Å³');

  const mixed = combineFields([{ field: a, weight: 1 }, { field: grid('e/bohr³', true), weight: 1 }]);
  assert.equal(mixed.valueUnit, null);
  assert.equal(combineFields([{ field: a, weight: 1 }, { field: grid(null, true), weight: 1 }]).valueUnit, null);

  // Two blocks keep the unit and stay blocks.
  const blockSum = combineFields([{ field: grid('e/bohr³', false), weight: 1 }, { field: grid('e/bohr³', false), weight: 1 }]);
  assert.equal(blockSum.valueUnit, 'e/bohr³');
  assert.equal(blockSum.periodic, false);
  assert.equal(magnitudeField([grid('e/bohr³', false), grid('e/bohr³', false)]).periodic, false);
  // A periodic field and a finite block do not share a frame: refused, not
  // silently drawn at the first one's place (feature branch rule).
  assert.throws(() => combineFields([{ field: a, weight: 1 }, { field: grid('e/Å³', false), weight: 1 }]),
    /periodic mismatch/);
  assert.throws(() => magnitudeField([a, grid('e/Å³', false)]), /periodic mismatch/);
});

test('Field: only a literal false makes a block (untrusted share payloads)', () => {
  for (const periodic of ['false', 0, null, undefined, 'no', {}]) {
    assert.equal(new Field({ periodic: /** @type {any} */ (periodic) }).periodic, true, String(periodic));
  }
  assert.equal(new Field({ periodic: false }).periodic, false);
});

test('NVal = 2 (not an orbital cube): "Value i" labels, no unit', () => {
  const text = [' Two values', ' Electron density from Total SCF Density',
    '    0 0 0 0    2', '    1 1 0 0', '    1 0 1 0', '    2 0 0 1', ' 1 2 3 4'].join('\n');
  const cube = parseCube(text);
  assert.deepEqual(cube.datasetLabels, ['Value 1', 'Value 2']);
  assert.deepEqual(cube.datasetUnits, [null, null]);
  assert.deepEqual([...cube.values[0]], [1, 3]);
  assert.deepEqual([...cube.values[1]], [2, 4]);
});

test('NVal = 4 cube of a potential: gradient labels but no density unit', () => {
  const text = [' t', ' Electrostatic potential and gradient', '    0 0 0 0    4', '    1 1 0 0', '    1 0 1 0',
    '    1 0 0 1', ' 1 2 3 4'].join('\n');
  const cube = parseCube(text);
  assert.deepEqual(cube.datasetLabels, GRADIENT_LABELS);
  assert.deepEqual(cube.datasetUnits, [null, null, null, null]);
});

test('labels are parallel to values for every kind of cube', () => {
  for (const name of ['bohr_density.cube', 'angstrom_potential.cube', 'gradient_nval4.cube',
    'mo_two_orbitals.cube', 'sheared_spin.cub']) {
    const cube = parseCube(fixture(name));
    assert.equal(cube.datasetLabels.length, cube.values.length, name);
    assert.equal(cube.datasetUnits.length, cube.values.length, name);
  }
});

test('hostile or malformed headers throw clear errors', () => {
  // Non-string input.
  assert.throws(() => parseCube(/** @type {any} */ (null)), /expected text content/);
  // A grid size far beyond what the text could hold is refused before any
  // array is allocated (1e9 points would otherwise be 4 GB).
  assert.throws(() => parseCube(' t\n d\n 0 0 0 0\n 1000 1 0 0\n 1000 0 1 0\n 1000 0 0 1\n 1 2 3\n'),
    /truncated.*expected 1000000000 data values.*at most/);
  // NVAL must be a positive integer.
  assert.throws(() => parseCube(' t\n d\n 0 0 0 0 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n 1\n'), /NVAL/);
  assert.throws(() => parseCube(' t\n d\n 0 0 0 0 1.5\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n 1\n'), /NVAL/);
  // A zero or fractional grid size names its axis.
  assert.throws(() => parseCube(' t\n d\n 0 0 0 0\n 1 1 0 0\n 0 0 1 0\n 1 0 0 1\n 1\n'), /grid size "0" on axis 2/);
  // A non-numeric origin.
  assert.throws(() => parseCube(' t\n d\n 0 0 x 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n 1\n'), /line 3/);
  // An orbital cube whose id record is cut short.
  assert.throws(() => parseCube(' t\n d\n -1 0 0 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n 1 1 0 0 0\n 3 5 6'),
    /2 of 3 ids/);
  // Empty text.
  assert.throws(() => parseCube(''), /line 3/);
});

test('blank lines between header lines are skipped', () => {
  const text = [' t', ' Electron density', '', '    0 0 0 0', '', '    1 1 0 0', '    1 0 1 0', '    2 0 0 1', '', ' 1', ' 2']
    .join('\n');
  assert.deepEqual([...parseCube(text).values[0]], [1, 2]);
});

test('a boilerplate "OUTER LOOP" second line falls back to the title for the label', () => {
  const text = [' H2 molecule, electron density', ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z',
    '    0 0 0 0', '    1 1 0 0', '    1 0 1 0', '    2 0 0 1', ' 0.1 0.2'].join('\n');
  const cube = parseCube(text);
  assert.deepEqual(cube.datasetLabels, ['H2 molecule, electron density']);
  assert.deepEqual(cube.datasetUnits, [CUBE_DENSITY_UNIT]);
});
