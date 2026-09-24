// Gaussian cube parser (docs/io/CubeParser.js) against the cubegen spec, plus
// the pieces around it that load in node: the '.cub' name rule in
// docs/io/formats.js and the valueUnit / periodic propagation of Field and
// CompositeField. Run with `node --test tools/unittest/`.
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
  parseCubeText, cubeDensityUnit, BOHR_TO_ANGSTROM, CUBE_DENSITY_UNIT,
} from '../../docs/io/CubeParser.js';
import { detectFormat, headOf } from '../../docs/io/formats.js';
import { Field } from '../../docs/model/Field.js';
import { combineFields, magnitudeField } from '../../docs/model/CompositeField.js';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const fixture = (name) => fs.readFileSync(path.join(HERE, 'fixtures', 'cube', name), 'utf8');
const B = BOHR_TO_ANGSTROM;

/** Stored index of grid point (i, j, k): axis 1 fastest. */
const at = (cube, i, j, k) => i + cube.counts[0] * (j + cube.counts[1] * k);

function assertClose(actual, expected, message, tol = 1e-9) {
  assert.ok(Math.abs(actual - expected) <= tol, `${message}: ${actual} vs ${expected}`);
}
function assertVecClose(actual, expected, message, tol = 1e-9) {
  assert.equal(actual.length, expected.length, `${message}: length`);
  actual.forEach((x, i) => assertClose(x, expected[i], `${message}[${i}]`, tol));
}

test('(a) bohr cube, NVal = 1: units, atoms, labels and x-fastest ordering', () => {
  const cube = parseCubeText(fixture('bohr_density.cube'));
  assert.equal(cube.title, 'Water fragment test cube');
  assert.equal(cube.description, 'Electron density from Total SCF Density');
  assert.equal(cube.lengthUnit, 'bohr');
  assert.deepEqual(cube.counts, [2, 3, 4]);
  assert.equal(cube.nval, 1);
  assert.equal(cube.moIndices, null);

  assertVecClose(cube.origin, [-B, -B, -B], 'origin in Å');
  assertVecClose(cube.voxel[0], [0.5 * B, 0, 0], 'voxel a');
  assertVecClose(cube.lattice[0], [1.0 * B, 0, 0], 'lattice a = 2 × 0.5 bohr');
  assertVecClose(cube.lattice[1], [0, 1.5 * B, 0], 'lattice b = 3 × 0.5 bohr');
  assertVecClose(cube.lattice[2], [0, 0, 2.0 * B], 'lattice c = 4 × 0.5 bohr');

  assert.equal(cube.atoms.length, 2);
  assert.equal(cube.atoms[0].Z, 8);
  assert.equal(cube.atoms[1].Z, 1);
  assertVecClose(cube.atoms[0].position, [0, 0, 0], 'O position');
  // One unit for all three coordinates (the old reader used a per-axis flag).
  assertVecClose(cube.atoms[1].position, [-0.5 * B, -0.5 * B, 0.5 * B], 'H position in Å');

  assert.equal(cube.fields.length, 1);
  const [field] = cube.fields;
  assert.equal(field.label, 'Electron density from Total SCF Density');
  assert.equal(field.valueUnit, CUBE_DENSITY_UNIT);
  assert.equal(field.values.length, 24);
  for (let i = 0; i < 2; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 4; k++) {
        assert.equal(field.values[at(cube, i, j, k)], 100 * i + 10 * j + k, `value at ${i},${j},${k}`);
      }
    }
  }
});

test('(b) Å cube (negative N): axes are not flipped, atoms stay in Å', () => {
  const cube = parseCubeText(fixture('angstrom_potential.cube'));
  assert.equal(cube.lengthUnit, 'angstrom');
  assert.deepEqual(cube.counts, [2, 2, 3]);
  assertVecClose(cube.voxel[0], [0.5, 0, 0], 'voxel a');
  assertVecClose(cube.lattice[0], [1, 0, 0], 'lattice a');
  assertVecClose(cube.lattice[1], [0, 1, 0], 'lattice b');
  assertVecClose(cube.lattice[2], [0, 0, 1.5], 'lattice c');
  assertVecClose(cube.atoms[0].position, [0.25, 0.5, 0.75], 'atom position');
  // A potential is not a density, whatever its description ends with.
  assert.equal(cube.fields[0].label, 'Electrostatic potential from Total SCF Density');
  assert.equal(cube.fields[0].valueUnit, null);
  assert.equal(cube.fields[0].values[at(cube, 1, 1, 2)], 112);
  assert.equal(cube.fields[0].values[at(cube, 0, 1, 2)], 12);
});

test('(c) NVal = 4 gradient cube: four de-interleaved fields', () => {
  const cube = parseCubeText(fixture('gradient_nval4.cube'));
  assert.equal(cube.nval, 4);
  assert.deepEqual(cube.fields.map((f) => f.label), ['Density', '∂ρ/∂x', '∂ρ/∂y', '∂ρ/∂z']);
  assert.deepEqual(cube.fields.map((f) => f.valueUnit), [CUBE_DENSITY_UNIT, null, null, null]);
  for (let v = 0; v < 4; v++) {
    for (let i = 0; i < 2; i++) {
      for (let j = 0; j < 2; j++) {
        for (let k = 0; k < 2; k++) {
          assert.equal(cube.fields[v].values[at(cube, i, j, k)], 100 * i + 10 * j + k + 0.25 * v,
            `value ${v} at ${i},${j},${k}`);
        }
      }
    }
  }
});

test('(d) orbital cube (NAtoms < 0): MO record read, one field per orbital', () => {
  const cube = parseCubeText(fixture('mo_two_orbitals.cube'));
  assert.equal(cube.atoms.length, 1);
  assert.deepEqual(cube.moIndices, [5, 6]);
  assert.equal(cube.nval, 2);
  assert.deepEqual(cube.fields.map((f) => f.label), ['MO 5', 'MO 6']);
  assert.deepEqual(cube.fields.map((f) => f.valueUnit), [null, null]);
  for (let i = 0; i < 2; i++) {
    for (let k = 0; k < 3; k++) {
      assert.equal(cube.fields[0].values[at(cube, i, 0, k)], 10 * i + k + 1);
      assert.equal(cube.fields[1].values[at(cube, i, 0, k)], -(10 * i + k + 1));
    }
  }
});

test('(e) truncated data throws instead of leaving zeros', () => {
  assert.throws(() => parseCubeText(fixture('truncated.cube')), /truncated.*expected 24 values.*found 20/);
  assert.throws(() => parseCubeText(' t\n d\n 1 0 0 0\n 2 1 0 0\n'), /axis 2/);
  assert.throws(() => parseCubeText(' t\n d\n 2 0 0 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n 1 1 0 0 0\n'),
    /atom line 2 of 2/);
  assert.throws(() => parseCubeText(' t\n d\n 0 0 0 0\n 1 1 0 0\n 1 0 1 0\n 1 0 0 1\n abc\n'),
    /"abc" is not a number/);
});

test('(f) non-orthogonal axes, CRLF line ends and Fortran D exponents', () => {
  const cube = parseCubeText(fixture('sheared_spin.cub'));
  assert.equal(cube.title, 'Sheared grid');
  assert.deepEqual(cube.counts, [3, 3, 2]);
  assertVecClose(cube.lattice[0], [1.5 * B, 0, 0], 'lattice a', 1e-6);
  assertVecClose(cube.lattice[1], [0.75 * B, 3 * 0.4330127 * B, 0], 'lattice b', 1e-6);
  assertVecClose(cube.lattice[2], [0.2 * B, 0.2 * B, 1.0 * B], 'lattice c', 1e-6);
  // The atom sits at origin + 1·v1 + 2·v2 (+0·v3) in bohr.
  const shifted = cube.atoms[0].position.map((c, i) => c - cube.origin[i]);
  const expected = [0, 1, 2].map((c) => cube.voxel[0][c] + 2 * cube.voxel[1][c]);
  assertVecClose(shifted, expected, 'atom offset from the origin', 1e-6);
  // Spin densities are densities.
  assert.equal(cube.fields[0].valueUnit, CUBE_DENSITY_UNIT);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 2; k++) {
        assert.equal(cube.fields[0].values[at(cube, i, j, k)], 100 * i + 10 * j + k);
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
  const cube = parseCubeText(text);
  assert.equal(cube.fields[0].label, 'Cube data');
  const expected = [1.5e-3, -225, 0.3, 0.12345e-100, 1.234567890123456789, -7, 0.5, 4];
  assert.deepEqual([...cube.fields[0].values], expected.map(Math.fround));
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

  const mixed = combineFields([{ field: a, weight: 1 }, { field: grid('e/bohr³', false), weight: 1 }]);
  assert.equal(mixed.valueUnit, null);
  assert.equal(mixed.periodic, false);
  assert.equal(combineFields([{ field: a, weight: 1 }, { field: grid(null, true), weight: 1 }]).valueUnit, null);
});

test('a boilerplate "OUTER LOOP" second line falls back to the title for the label', () => {
  const text = [' H2 molecule, electron density', ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z',
    '    0 0 0 0', '    1 1 0 0', '    1 0 1 0', '    2 0 0 1', ' 0.1 0.2'].join('\n');
  const [field] = parseCubeText(text).fields;
  assert.equal(field.label, 'H2 molecule, electron density');
  assert.equal(field.valueUnit, CUBE_DENSITY_UNIT);
});
