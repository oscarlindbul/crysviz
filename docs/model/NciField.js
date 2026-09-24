import { Field } from './Field.js';
import { computeFieldStats } from './CompositeField.js';
import { runNci } from '../workers/nciTasks.js';
import { invert3x3, transpose3x3, fracToCartPoint, cartToFractional } from '../math/index.js';

/**
 * NCI (non-covalent interaction) fields from a density grid.
 *
 * The NCI plot (Johnson et al., JACS 2010) is two grids on the density's own
 * grid: the reduced density gradient s, drawn as an isosurface at s = 0.5, and
 * sign(λ₂)ρ, which colours that surface. The number crunching is in
 * workers/nciTasks.js (a WASM port of Jmol's NciCalculation); this module only
 * turns a selected field and a structure into the request and the result into
 * two ordinary `Field`s, built the way `combineFields` builds a derived field —
 * same grid, `derivedFrom` pointing at the source — so isosurfaces, cut planes
 * and the catalog treat them like anything read from a file.
 *
 * The s field carries `colorBy` pointing at the sign(λ₂)ρ field, so selecting
 * it draws the NCI plot straight away; nothing else in the app has to know
 * that the two belong together.
 */

/** Bohr in Å (CODATA 2018). e/Å³ × BOHR³ = e/bohr³. */
const BOHR_IN_ANGSTROM = 0.529177210903;

/** The calculation works in atomic units; this is what turns e/Å³ into e/bohr³. */
const VALUE_SCALE_FROM_UNIT = {
  'e/bohr³': 1,
  'e/Å³': BOHR_IN_ANGSTROM ** 3,
};

/** Colour range of the NCI plot, in e/bohr³: the one Jmol and NCIPLOT use. */
export const NCI_COLOR_RANGE = Object.freeze({ min: -0.04, max: 0.04 });

/** The isosurface level of the NCI plot. */
export const NCI_S_ISOVALUE = 0.5;

/**
 * Element symbols by atomic number (index 0 unused). Kept here rather than
 * imported from io/ReadCubeModule.js's table, which would drag the renderer
 * into the model layer.
 */
const SYMBOLS = ('_ H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu '
  + 'Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr '
  + 'Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac '
  + 'Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split(' ');

/** Element symbol → atomic number. */
const Z_OF = Object.fromEntries(SYMBOLS.map((symbol, z) => [symbol, z]).slice(1));

/** Promolecular densities are tabulated for H–Ar; heavier atoms are treated as Ar. */
const PROMOLECULAR_MAX_Z = 18;

/**
 * An element symbol as the file wrote it ("FE", "fe", "Fe1", "Fe_pv") → its
 * atomic number, or 0 when it is not an element (a dummy "X", a vacancy).
 * @param {string} symbol
 */
export function atomicNumberOf(symbol) {
  const letters = String(symbol ?? '').trim().match(/^[A-Za-z]{1,2}/)?.[0];
  if (!letters) return 0;
  const normalized = letters[0].toUpperCase() + letters.slice(1).toLowerCase();
  // "Fe1" and "Fe_pv" match as two letters; a one-letter symbol followed by
  // another letter ("Cx") falls back to the single letter when that is an
  // element and the pair is not.
  return Z_OF[normalized] ?? Z_OF[normalized[0]] ?? 0;
}

/** Lattice rows the grid spans: n_i × voxel[i]. */
function gridLattice(field) {
  const counts = [field.nx, field.ny, field.nz];
  return field.voxel.map((row, i) => row.map((c) => c * counts[i]));
}

/** True when two lattices agree to within a relative 1e-6. */
function latticesAgree(a, b) {
  if (!Array.isArray(a) || !Array.isArray(b) || a.length !== 3 || b.length !== 3) return false;
  let scale = 0;
  for (const row of a) for (const v of row) scale = Math.max(scale, Math.abs(v));
  const tol = 1e-6 * Math.max(scale, 1);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      if (Math.abs((a[i]?.[j] ?? NaN) - (b[i]?.[j] ?? NaN)) > tol) return false;
    }
  }
  return true;
}

/**
 * The structure's atoms as the promolecular request wants them: atomic numbers
 * and positions fractional with respect to the FIELD's grid cell.
 *
 * Atom positions are stored fractional with respect to `structure.lattice`.
 * For a cube or CHGCAR that is the grid's own cell and they pass straight
 * through; when the field came from another cell (a WAVECAR attached to a
 * structure of a different shape), they go through Cartesian into the grid's
 * frame, which is also the frame the isosurface is drawn in (model/Isosurface.js
 * places vertices with the field's voxel vectors).
 *
 * @returns {{Z: number[], frac: number[][], ignored: string[], clamped: string[]}}
 */
function promolecularAtoms(structure, field) {
  const atoms = structure?.atoms ?? [];
  const fieldLattice = gridLattice(field);
  const sameCell = latticesAgree(structure?.lattice, fieldLattice);
  const toFieldInverse = sameCell ? null : invert3x3(transpose3x3(fieldLattice));

  const Z = [];
  const frac = [];
  const ignored = new Set();
  const clamped = new Set();

  atoms.forEach((atom, index) => {
    const symbol = structure.elements?.[index]
      ?? atom.getRepresentativeElement?.()
      ?? atom.species?.[0]?.element
      ?? '';
    const z = atomicNumberOf(symbol);
    const position = atom.position;
    if (!z || !Array.isArray(position) || position.length < 3 || !position.every(Number.isFinite)) {
      ignored.add(symbol || '?');
      return;
    }
    if (z > PROMOLECULAR_MAX_Z) clamped.add(SYMBOLS[z] ?? symbol);
    Z.push(z);
    frac.push(sameCell
      ? [position[0], position[1], position[2]]
      : cartToFractional(fracToCartPoint(position, structure.lattice), fieldLattice, toFieldInverse));
  });

  return { Z, frac, ignored: [...ignored], clamped: [...clamped] };
}

/** Every label the catalog already uses, loaded or not. */
function labelsInCatalog(catalog) {
  const labels = new Set();
  if (!catalog) return labels;
  for (const node of catalog.walk()) {
    if (node.label) labels.add(node.label);
    const field = node.peek?.();
    if (field?.label) labels.add(field.label);
  }
  return labels;
}

/**
 * The two labels, made unique as a pair. Planes and storage find fields by
 * label, so running NCI twice on the same source must not produce two fields
 * both called the same thing; the suffix is shared so the pair stays
 * recognisable as a pair.
 */
function uniqueLabels(bases, catalog) {
  const taken = labelsInCatalog(catalog);
  for (let n = 1; ; n++) {
    const suffix = n === 1 ? '' : ` [${n}]`;
    const labels = bases.map((base) => `${base}${suffix}`);
    if (labels.every((label) => !taken.has(label))) return labels;
  }
}

/**
 * Build the SCF- or promolecular-NCI fields on `sourceField`'s grid.
 *
 * @param {import('./Field.js').Field} sourceField the density (SCF) or just the
 *   grid to evaluate on (promolecular)
 * @param {object} options
 * @param {'scf'|'promolecular'} options.kind
 * @param {any} [options.structure] needed for 'promolecular'
 * @param {import('./FieldCatalog.js').FieldCatalog | null} [options.catalog]
 *   used only to keep the new labels unique
 * @param {{rhoMin?: number, rhoPlot?: number, sCap?: number}} [options.params]
 * @returns {Promise<{sField: Field, colourField: Field, info: {
 *   kind: string, included: number, gridPoints: number, valueScale: number,
 *   unitWarning: string | null, clampedElements: string[], ignoredElements: string[],
 *   atomCount: number}}>}
 */
export async function createNciFields(sourceField, { kind, structure = null, catalog = null, params } = /** @type {any} */ ({})) {
  if (kind !== 'scf' && kind !== 'promolecular') {
    throw new Error(`createNciFields: unknown kind "${kind}"`);
  }
  if (!sourceField) throw new Error('Select a field first: NCI is computed on its grid.');

  const { nx, ny, nz, voxel } = sourceField;
  const gridPoints = nx * ny * nz;
  if (!(gridPoints > 0) || !Array.isArray(voxel) || voxel.length !== 3) {
    throw new Error(`"${sourceField.label}" has no usable grid.`);
  }
  // Central differences need a neighbour on each side, and a non-periodic grid
  // loses its boundary layer on top of that.
  if (Math.min(nx, ny, nz) < 3) {
    throw new Error(`"${sourceField.label}" is ${nx}×${ny}×${nz}; NCI needs at least 3 points along each axis.`);
  }

  const periodic = sourceField.periodic ?? true;
  const info = {
    kind,
    included: 0,
    gridPoints,
    valueScale: 1,
    unitWarning: /** @type {string | null} */ (null),
    clampedElements: /** @type {string[]} */ ([]),
    ignoredElements: /** @type {string[]} */ ([]),
    atomCount: 0,
  };

  /** @type {any} */
  const request = {
    kind,
    nx,
    ny,
    nz,
    voxel: voxel.map((row) => [...row]),
    periodic,
  };
  if (params) request.params = { ...params };

  if (kind === 'scf') {
    const values = sourceField.values;
    if (!values || values.length < gridPoints) {
      throw new Error(`"${sourceField.label}" has no values loaded.`);
    }
    const unit = sourceField.valueUnit ?? null;
    if (unit && VALUE_SCALE_FROM_UNIT[unit] !== undefined) {
      info.valueScale = VALUE_SCALE_FROM_UNIT[unit];
    } else {
      info.valueScale = 1;
      info.unitWarning = unit
        ? `"${sourceField.label}" is in ${unit}, not a density unit; its values were used as e/bohr³.`
        : `The unit of "${sourceField.label}" is unknown; its values were assumed to be e/bohr³.`;
    }
    request.values = values;
    request.valueScale = info.valueScale;
  } else {
    if (!structure) throw new Error('Promolecular NCI needs a structure with atoms.');
    const atoms = promolecularAtoms(structure, sourceField);
    info.ignoredElements = atoms.ignored;
    info.atomCount = atoms.Z.length;
    if (atoms.Z.length === 0) throw new Error('The structure has no atoms to build a promolecular density from.');
    request.atoms = { Z: atoms.Z, frac: atoms.frac };
    info.clampedElements = atoms.clamped;
  }

  const result = await runNci(request);
  info.included = result.included ?? 0;
  // The backend reports what it clamped too; merge rather than trust one list.
  if (Array.isArray(result.clampedElements)) {
    info.clampedElements = [...new Set([...info.clampedElements, ...result.clampedElements])];
  }

  const sourceLabel = sourceField.label || 'field';
  const bases = kind === 'scf'
    ? [`Reduced density gradient s (SCF-NCI of ${sourceLabel})`, `sign(λ₂)ρ (SCF-NCI of ${sourceLabel})`]
    : ['Reduced density gradient s (promolecular NCI)', 'sign(λ₂)ρ (promolecular NCI)'];
  const [sLabel, colourLabel] = uniqueLabels(bases, catalog);

  const colourField = new Field({
    nx,
    ny,
    nz,
    origin: sourceField.origin,
    voxel: sourceField.voxel,
    values: result.sl2rho,
    component: sourceField.component ?? 0,
    label: colourLabel,
    // Signed by construction: attractive regions are negative, steric ones
    // positive, so the pair of ± surfaces is the natural first view.
    useAbsoluteIsoValue: false,
    isoValue: 0.02,
    valueUnit: 'e/bohr³',
    periodic,
    ...computeFieldStats(result.sl2rho),
  });
  colourField.derivedFrom = [{ field: sourceField, weight: 1 }];
  colourField.derivedOp = 'nci-sl2rho';

  const sField = new Field({
    nx,
    ny,
    nz,
    origin: sourceField.origin,
    voxel: sourceField.voxel,
    values: result.s,
    component: sourceField.component ?? 0,
    label: sLabel,
    useAbsoluteIsoValue: false,
    isoValue: NCI_S_ISOVALUE,
    valueUnit: null, // s is dimensionless
    periodic,
    colorBy: {
      field: colourField,
      colormap: 'bgyor',
      min: NCI_COLOR_RANGE.min,
      max: NCI_COLOR_RANGE.max,
    },
    ...computeFieldStats(result.s),
  });
  sField.derivedFrom = [{ field: sourceField, weight: 1 }];
  sField.derivedOp = 'nci-s';

  return { sField, colourField, info };
}
