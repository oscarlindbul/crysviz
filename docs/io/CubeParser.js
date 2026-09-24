// CubeParser.js
// Pure parser for Gaussian cube files: text in, plain data out. No imports, so
// it loads in node (tools/unittest/cubereader.test.mjs) as well as the browser;
// io/ReadCubeModule.js turns its result into a Structure and Fields.
//
// Format (gaussian.com/cubegen, "Output File Formats"):
//   line 1   title (free text)
//   line 2   description, e.g. "Electron density from Total SCF Density"
//   line 3   NAtoms  X0 Y0 Z0  [NVal]      origin; NVal values per point, default 1
//   line 4-6 N_i  v_ix v_iy v_iz            N_i > 0: bohr, N_i < 0: Å (|N_i| points)
//   |NAtoms| lines  Z  charge  x y z        same length unit as the axes
//   if NAtoms < 0 (orbital cube): NMO MO(1) .. MO(NMO), possibly over several
//            lines; NVal = NMO
//   data     N1·N2 records of N3·NVal values (6 per line, each record starting
//            on a new line). Slowest index is along axis 1, fastest along
//            axis 3, and the NVal values of one point are adjacent.

export const BOHR_TO_ANGSTROM = 0.529177249;

/** Unit of electron densities as a cube stores them (atomic units). */
export const CUBE_DENSITY_UNIT = 'e/bohr³';

/** Field labels for the four values of a Gaussian `Gradient` cube. */
export const GRADIENT_LABELS = ['Density', '∂ρ/∂x', '∂ρ/∂y', '∂ρ/∂z'];

const POW10 = [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
  1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22];

/**
 * Forward-only scanner over the text of a cube file.
 *
 * The data block of a real cube runs to millions of numbers; splitting it into
 * lines and each line into tokens allocated two strings per number and took
 * most of the load time. This walks the string once, parsing plain decimal
 * numbers straight from the character codes (Fortran `D` exponents included)
 * and handing anything unusual (more than 15 significant digits, a Fortran
 * exponent without its `E`, `NaN`) to `Number` on just that token.
 */
class CubeScanner {
  /** @param {string} text */
  constructor(text) {
    this.text = text;
    this.pos = 0;
    /** Set once `number()` / `fill()` has run off the end of the text. */
    this.exhausted = false;
  }

  /**
   * The next line (without its terminator), or null at the end of the text.
   * @returns {string | null}
   */
  line() {
    const text = this.text;
    if (this.pos >= text.length) return null;
    let end = text.indexOf('\n', this.pos);
    if (end === -1) end = text.length;
    let stop = end;
    if (stop > this.pos && text.charCodeAt(stop - 1) === 13) stop--; // CRLF
    const out = text.substring(this.pos, stop);
    this.pos = end + 1;
    return out;
  }

  /**
   * The next whitespace-separated number, or NaN when the text is exhausted
   * (`exhausted` is then set). Throws on a token that is not a number.
   * @returns {number}
   */
  number() {
    const one = this._one || (this._one = new Float64Array(1));
    return this.fill(one, 0, 1) === 1 ? one[0] : NaN;
  }

  /**
   * Read up to `count` numbers into `out` from `offset` on, in file order.
   * One tight loop with everything in locals, since this is where a large
   * file spends its load time.
   * @param {Float32Array | Float64Array} out
   * @param {number} offset
   * @param {number} count
   * @returns {number} how many were read; fewer than `count` only at the end
   *   of the text (and `exhausted` is then set)
   */
  fill(out, offset, count) {
    const text = this.text;
    const len = text.length;
    let p = this.pos;
    let n = 0;
    while (n < count) {
      let c = 0;
      while (p < len && (c = text.charCodeAt(p)) <= 32) p++;
      if (p >= len) {
        this.exhausted = true;
        break;
      }
      const start = p;

      let negative = false;
      if (c === 45 /* - */ || c === 43 /* + */) {
        negative = c === 45;
        c = ++p < len ? text.charCodeAt(p) : 0;
      }
      let mantissa = 0;
      let digits = 0; // significant digits in the mantissa
      let scale = 0;
      let sawDigit = false;
      while (c >= 48 && c <= 57) {
        mantissa = mantissa * 10 + (c - 48);
        if (mantissa !== 0) digits++;
        sawDigit = true;
        c = ++p < len ? text.charCodeAt(p) : 0;
      }
      if (c === 46 /* . */) {
        c = ++p < len ? text.charCodeAt(p) : 0;
        while (c >= 48 && c <= 57) {
          mantissa = mantissa * 10 + (c - 48);
          if (mantissa !== 0) digits++;
          scale--;
          sawDigit = true;
          c = ++p < len ? text.charCodeAt(p) : 0;
        }
      }
      if (sawDigit && (c === 69 || c === 101 || c === 68 || c === 100) /* E e D d */) {
        c = ++p < len ? text.charCodeAt(p) : 0;
        let expNegative = false;
        if (c === 45 || c === 43) {
          expNegative = c === 45;
          c = ++p < len ? text.charCodeAt(p) : 0;
        }
        let exponent = 0;
        let expDigits = 0;
        while (c >= 48 && c <= 57) {
          exponent = exponent * 10 + (c - 48);
          expDigits++;
          c = ++p < len ? text.charCodeAt(p) : 0;
        }
        if (expDigits === 0) sawDigit = false;
        scale += expNegative ? -exponent : exponent;
      }

      let value;
      if (sawDigit && digits <= 15 && c <= 32) {
        // Fast path: a well-formed decimal that a double holds exactly before
        // the final power of ten, ending at whitespace or the end of the text.
        if (scale === 0) value = mantissa;
        else if (scale < 0) value = scale >= -22 ? mantissa / POW10[-scale] : mantissa / Math.pow(10, -scale);
        else value = scale <= 22 ? mantissa * POW10[scale] : mantissa * Math.pow(10, scale);
        if (negative) value = -value;
      } else {
        p = start;
        while (p < len && text.charCodeAt(p) > 32) p++;
        value = parseToken(text.substring(start, p), start);
      }
      out[offset + n++] = value;
    }
    this.pos = p;
    return n;
  }
}

/**
 * The slow path of `CubeScanner.fill`: a token the fast path did not accept.
 * @param {string} token
 * @param {number} at character offset, for the error message
 * @returns {number}
 */
function parseToken(token, at) {
  let value = Number(token.replace(/[dD]/, 'e'));
  // Fortran drops the E of a three-digit exponent: 0.12345-100.
  const bare = Number.isNaN(value) && /^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$/.exec(token);
  if (bare) value = Number(`${bare[1]}e${bare[2]}`);
  if (Number.isNaN(value) && !/^[+-]?nan$/i.test(token)) {
    throw new Error(`Cube file: "${token}" is not a number (at character ${at})`);
  }
  return value;
}

/**
 * Split a header line into numbers; throws naming the line when it is short.
 * @param {string | null} line
 * @param {number} minCount
 * @param {string} what description for the error message
 * @returns {number[]}
 */
function headerNumbers(line, minCount, what) {
  const tokens = line == null ? [] : line.trim().split(/\s+/).filter(Boolean);
  const numbers = tokens.map((t) => Number(t.replace(/[dD]/, 'e')));
  let count = 0;
  while (count < numbers.length && Number.isFinite(numbers[count])) count++;
  if (count < minCount) {
    throw new Error(`Cube file: expected ${what}, got ${line == null ? 'end of file' : `"${line.trim()}"`}`);
  }
  return numbers.slice(0, count);
}

/**
 * The unit of a cube's values, judged from its description line.
 *
 * Gaussian names a cube "<quantity> from <density matrix>" ("Electron density
 * from Total SCF Density", "Electrostatic potential from Total SCF Density"),
 * so only the part before " from " says what the numbers are. Electron, spin,
 * alpha and beta densities are all e/bohr³; potentials, Laplacians, gradients,
 * orbitals and the like are not densities and get no unit.
 *
 * @param {string} description line 2 of the file
 * @returns {string | null}
 */
export function cubeDensityUnit(description) {
  const quantity = String(description || '').toLowerCase().split(/\s+from\s+/)[0];
  if (!/densit/.test(quantity)) return null;
  if (/potential|laplacian|gradient|reduced|orbital|localiz|\belf\b|\bmo\b|norm/.test(quantity)) return null;
  return CUBE_DENSITY_UNIT;
}

/**
 * @typedef {object} CubeFieldData
 * @property {string} label
 * @property {string | null} valueUnit
 * @property {Float32Array} values  x-fastest: index = i + n1*(j + n2*k)
 */

/**
 * @typedef {object} CubeData
 * @property {string} title        line 1, trimmed
 * @property {string} description  line 2, trimmed
 * @property {'bohr' | 'angstrom'} lengthUnit  unit the file was written in
 * @property {number[]} counts     [n1, n2, n3] grid points per axis
 * @property {number[]} origin     grid origin, Å
 * @property {number[][]} voxel    step vector of each axis, Å
 * @property {number[][]} lattice  counts[i] · voxel[i], Å
 * @property {{Z: number, charge: number, position: number[]}[]} atoms
 *   Cartesian positions in Å, as written (not shifted by the origin)
 * @property {number[] | null} moIndices  orbital numbers of an MO cube, else null
 * @property {number} nval         values per grid point
 * @property {CubeFieldData[]} fields  one per value, de-interleaved
 */

/**
 * Parse the text of a Gaussian cube file.
 *
 * @param {string} text
 * @returns {CubeData}
 */
export function parseCubeText(text) {
  if (typeof text !== 'string') throw new Error('Cube file: expected text content');
  const scan = new CubeScanner(text);

  const title = (scan.line() ?? '').trim();
  const description = (scan.line() ?? '').trim();

  const originLine = headerNumbers(scan.line(), 4, 'NAtoms and the origin on line 3');
  const nAtomsSigned = originLine[0];
  if (!Number.isInteger(nAtomsSigned)) throw new Error(`Cube file: NAtoms "${nAtomsSigned}" is not an integer`);
  const isMO = nAtomsSigned < 0;
  const nAtoms = Math.abs(nAtomsSigned);
  let nval = originLine.length > 4 ? originLine[4] : 1;
  if (!Number.isInteger(nval) || nval < 1) throw new Error(`Cube file: NVal "${nval}" must be a positive integer`);

  const signedCounts = [];
  const steps = [];
  for (let axis = 0; axis < 3; axis++) {
    const row = headerNumbers(scan.line(), 4, `the axis ${axis + 1} line (N and a step vector)`);
    if (!Number.isInteger(row[0]) || row[0] === 0) {
      throw new Error(`Cube file: grid count "${row[0]}" on axis ${axis + 1} must be a nonzero integer`);
    }
    signedCounts.push(row[0]);
    steps.push(row.slice(1, 4));
  }
  // Negative counts mean Å. The files in the wild use one sign for all three
  // axes; like Jmol, any negative count switches the whole file to Å, since
  // the atom block has no sign of its own and must share one unit.
  const angstrom = signedCounts.some((n) => n < 0);
  const toAngstrom = angstrom ? 1 : BOHR_TO_ANGSTROM;
  const counts = signedCounts.map(Math.abs);
  const voxel = steps.map((v) => v.map((c) => c * toAngstrom));
  const lattice = voxel.map((v, i) => v.map((c) => c * counts[i]));
  const origin = originLine.slice(1, 4).map((c) => c * toAngstrom);

  const atoms = [];
  for (let a = 0; a < nAtoms; a++) {
    const row = headerNumbers(scan.line(), 5, `atom line ${a + 1} of ${nAtoms} (Z, charge, x, y, z)`);
    atoms.push({
      Z: Math.round(row[0]),
      charge: row[1],
      position: [row[2] * toAngstrom, row[3] * toAngstrom, row[4] * toAngstrom],
    });
  }

  let moIndices = null;
  if (isMO) {
    // NMO, MO(1..NMO), in 10I5 format so it may wrap over several lines.
    const nmo = scan.number();
    if (!Number.isInteger(nmo) || nmo < 1) {
      throw new Error(`Cube file: orbital cube (NAtoms < 0) needs an orbital count after the atoms, got "${nmo}"`);
    }
    moIndices = [];
    for (let m = 0; m < nmo; m++) {
      const mo = scan.number();
      if (!Number.isInteger(mo)) throw new Error(`Cube file: orbital list ends after ${m} of ${nmo} entries`);
      moIndices.push(mo);
    }
    nval = nmo;
  }

  const [n1, n2, n3] = counts;
  const npoints = n1 * n2 * n3;
  const expected = npoints * nval;
  /** @type {Float32Array[]} */
  const values = [];
  for (let v = 0; v < nval; v++) values.push(new Float32Array(npoints));

  // File order: axis 1 slowest, axis 3 fastest, values of a point adjacent.
  // Stored order: axis 1 fastest (index = i + n1*(j + n2*k)), which is what
  // Field.getValueAt and the marching-cubes code expect. Read one record (a
  // fixed i, j: all k and all values) at a time, then scatter it.
  const n12 = n1 * n2;
  const recordLength = n3 * nval;
  const record = new Float32Array(recordLength);
  let read = 0;
  for (let i = 0; i < n1; i++) {
    for (let j = 0; j < n2; j++) {
      const got = scan.fill(record, 0, recordLength);
      read += got;
      if (got < recordLength) {
        throw new Error(`Cube file: data is truncated — expected ${expected} values `
          + `(${n1}×${n2}×${n3} points × ${nval}), found ${read}`);
      }
      const base = i + n1 * j;
      if (nval === 1) {
        const out = values[0];
        for (let k = 0, index = base; k < n3; k++, index += n12) out[index] = record[k];
      } else {
        for (let v = 0; v < nval; v++) {
          const out = values[v];
          for (let k = 0, index = base; k < n3; k++, index += n12) out[index] = record[k * nval + v];
        }
      }
    }
  }

  return {
    title,
    description,
    lengthUnit: angstrom ? 'angstrom' : 'bohr',
    counts,
    origin,
    voxel,
    lattice,
    atoms,
    moIndices,
    nval,
    fields: describeFields(values, { title, description, moIndices }),
  };
}

/**
 * Label each de-interleaved value and give it a unit.
 * @param {Float32Array[]} values
 * @param {{title: string, description: string, moIndices: number[] | null}} meta
 * @returns {CubeFieldData[]}
 */
function describeFields(values, { title, description, moIndices }) {
  if (moIndices) {
    return values.map((v, m) => ({ label: `MO ${moIndices[m]}`, valueUnit: null, values: v }));
  }
  if (values.length === 1) {
    // Many non-Gaussian writers (VMD, CP2K, Quantum ESPRESSO's pp.x, ASE) put
    // the fixed "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" on line 2 and
    // what the file holds on line 1; that line says nothing about the field.
    const informative = description && !/outer\s+loop/i.test(description) ? description : '';
    const label = informative || title || 'Cube data';
    return [{ label, valueUnit: cubeDensityUnit(label), values: values[0] }];
  }
  if (values.length === 4) {
    // A Gaussian `Gradient` cube: the density followed by its gradient. The
    // description names the gradient, so the density unit is judged from
    // whether a density is mentioned at all (and not, say, a potential).
    const text = description.toLowerCase();
    const densityUnit = /densit/.test(text) && !/potential|laplacian/.test(text) ? CUBE_DENSITY_UNIT : null;
    return values.map((v, i) => ({ label: GRADIENT_LABELS[i], valueUnit: i === 0 ? densityUnit : null, values: v }));
  }
  return values.map((v, i) => ({ label: `Value ${i + 1}`, valueUnit: null, values: v }));
}
