// ReadCubeModule.js
// Gaussian .cube / .cub volumetric files → Structure + FieldContainer.
// The text parsing itself lives in CubeParser.js (pure, node-testable).
// Exports: readCubeFile(), PT, Bohr2Angstrom
//
import { Structure } from '../model/index.js';
import { invert3x3, transpose3x3, cartToFractional, normalizeFractional } from '../math/index.js';
import { runPeriodicWrapped } from '../render/index.js';
import { Field } from '../model/index.js';
import { FieldContainer } from '../model/index.js';
import { computeFieldStats } from '../model/index.js';
import { Atom } from '../model/index.js';
import { generateID } from '../utils/index.js';
import { parseCubeText, BOHR_TO_ANGSTROM } from './CubeParser.js';


//------------------------------------------------------------
//  Periodic table (lookup table for cube files as it contains 
//                  only the element number) ! Check if everythig is correct!
//------------------------------------------------------------
export const PT = {
  1: "H",   2: "He",
  3: "Li",  4: "Be",  5: "B",   6: "C",   7: "N",   8: "O",   9: "F",   10: "Ne",
  11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",  17: "Cl", 18: "Ar",
  19: "K",  20: "Ca", 21: "Sc", 22: "Ti", 23: "V",  24: "Cr", 25: "Mn", 26: "Fe",
  27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
  35: "Br", 36: "Kr",
  37: "Rb", 38: "Sr", 39: "Y",  40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru",
  45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te",
  53: "I",  54: "Xe",
  55: "Cs", 56: "Ba",
  // Lanthanides
  57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
  65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
  // Transition continues
  72: "Hf", 73: "Ta", 74: "W",  75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au",
  80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn",
  87: "Fr", 88: "Ra",
  // Actinides
  89: "Ac", 90: "Th", 91: "Pa", 92: "U",  93: "Np", 94: "Pu", 95: "Am", 96: "Cm",
  97: "Bk", 98: "Cf", 99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr",
  // Final row
  104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds",
  111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts",
  118: "Og"
};

export const Bohr2Angstrom = BOHR_TO_ANGSTROM; // conversion factor from Bohr to Angstroms

//------------------------------------------------------------
//  readCubeFile(content, fileName) → { fileName, structure_with_field }
//------------------------------------------------------------
/**
 * Read a Gaussian cube file into a Structure carrying its fields.
 *
 * The text is parsed by `parseCubeText` (io/CubeParser.js, which follows the
 * cubegen specification: Å or bohr axes, NVal values per point, orbital cubes);
 * this turns that plain data into the model objects. One Field per value per
 * point: a plain cube gives one, a Gradient cube four (density and its three
 * derivatives), an orbital cube one per orbital.
 *
 * The cell of the structure is the grid box (N_i · step_i), with the atoms
 * shifted by the grid origin, and every field's origin is [0,0,0] in that
 * frame. Cube grids are molecular boxes rather than periodic cells, so the
 * fields are marked `periodic: false`. Errors (a malformed header, truncated
 * data) are thrown.
 *
 * @param {string} content
 * @param {string} fileName
 */
export function readCubeFile(content, fileName) {
  const cube = parseCubeText(content);
  const structure = buildCubeStructure(cube);

  const [nx, ny, nz] = cube.counts;
  const fields = cube.fields.map(({ label, valueUnit, values }, component) => new Field({
    nx,
    ny,
    nz,
    origin: [0, 0, 0],
    voxel: cube.voxel.map((row) => [...row]),
    values,
    component,
    label,
    valueUnit,
    periodic: false,
    // One pass instead of the four separate `reduce` walks this used to do
    // over an array that runs to millions of entries.
    ...computeFieldStats(values),
  }));

  const container = new FieldContainer({
    fileName: fileName,
    // An orbital cube says so where the field browser shows the source, so a
    // list of "MO n" entries is not mistaken for anything else.
    source: cube.moIndices ? 'Cube (orbitals)' : 'Cube',
    fields: fields,
    fieldCount: fields.length
  });

  structure.volumetricFields = container; // Attach field container to structure for easy access in rendering
  return {
    fileName,
    structure_with_field: structure
  };
}

/**
 * The Structure of a parsed cube: the grid box as the cell, atoms shifted by
 * the grid origin.
 * @param {import('./CubeParser.js').CubeData} cube
 * @returns {Structure}
 */
function buildCubeStructure(cube) {
  const lattice = cube.lattice.map((row) => [...row]);
  const elements = cube.atoms.map((atom) => PT[atom.Z] || "X");
  const positions_cart = cube.atoms.map((atom) => atom.position.map((c, i) => c - cube.origin[i]));

  // --- convert cart → frac
  const latticeInverse = invert3x3(transpose3x3(lattice));
  const positions = (
    positions_cart.map(vec => cartToFractional(vec, lattice, latticeInverse))
  ).map(pos => pos.map(normalizeFractional));

  const atoms = positions.map((pos, i) => new Atom({
    position: pos,
    element: elements[i],
    uuid: generateID([elements[i]])
  }));

  const periodic = runPeriodicWrapped(
    { hash: "None", wrapped: {} },
    positions,
    elements,
    lattice
  );

  return new Structure({
    elements: elements,
    uniqueElements: [...new Set(elements)],
    lattice: lattice,
    atoms: atoms,
    periodic: periodic
  });
}
