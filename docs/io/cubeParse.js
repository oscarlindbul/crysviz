/**
 * Gaussian .cube parser: text in, plain geometry + grid data out.
 *
 * The format-level half of cube loading, split out of `io/ReadCubeModule.js`
 * so it imports nothing and can be unit-tested under node. It builds no model
 * classes — ReadCubeModule turns the result into a Structure and Fields.
 *
 * Layout (gaussian.com/cubegen, "Output File Formats"):
 *   line 1-2   free-text comments (either may be blank); line 2 usually says
 *              what the file holds ("Electron density from Total SCF Density")
 *   line 3     NATOMS  ox oy oz  [NVAL]     NVAL values per point, default 1
 *   line 4-6   Ni  step vector i            (one line per grid axis)
 *   |NATOMS|   Z  charge  x y z              (one line per atom)
 *   if NATOMS < 0:  M  id1 .. idM            (orbital ids, may wrap lines;
 *                                             NVAL = M)
 *   data       x slowest, z fastest; NVAL (or M) values per grid point,
 *              each (x, y) record starting on a new line
 *
 * Units: a positive N1 means the step vectors, origin and atom coordinates are
 * in Bohr; a negative N1 means Angstrom, and the grid size is |N1|. The units
 * are one choice for the whole file, so N1 decides it and the signs of N2/N3
 * are only stripped. Everything returned is in Angstrom. The data values are
 * passed through untouched (a Gaussian density stays in e/Bohr^3); the unit
 * they are in is judged from the description line (`datasetUnits`).
 *
 * The text is untrusted: every header number is checked, a token that is not
 * a number throws naming it, and a data section shorter than the header
 * promises throws instead of leaving a partly zero grid.
 */

/** CODATA 2022 Bohr radius in Angstrom. */
export const BOHR_TO_ANGSTROM = 0.529177210544;

/** Unit of electron densities as a cube stores them (atomic units). */
export const CUBE_DENSITY_UNIT = 'e/bohr³';

/** Dataset labels for the four values of a Gaussian `Gradient` cube. */
export const GRADIENT_LABELS = ['Density', '∂ρ/∂x', '∂ρ/∂y', '∂ρ/∂z'];

/**
 * @typedef {Object} CubeData
 * @property {string} label           the two comment lines, joined
 * @property {'bohr'|'angstrom'} units the units the file was written in
 * @property {number[]} origin        grid origin, Angstrom
 * @property {number[]} grid          [n1, n2, n3] grid points per axis (positive)
 * @property {number[][]} voxel       step vector per axis (rows), Angstrom
 * @property {number[][]} lattice     grid cell rows n_i * step_i, Angstrom
 * @property {{atomicNumber: number, charge: number, position: number[]}[]} atoms
 *                                    positions are absolute Cartesian, Angstrom
 * @property {number[]|null} datasetIds orbital/dataset ids (NATOMS < 0), else null
 * @property {Float32Array[]} values  one array per dataset, index x + n1*(y + n2*z)
 * @property {string[]} datasetLabels one per `values` entry: `MO n` for an
 *                                    orbital cube, GRADIENT_LABELS for NVAL = 4,
 *                                    `Value i` for other NVAL > 1, and for NVAL = 1
 *                                    the description line (the title when that is
 *                                    "OUTER LOOP" boilerplate, 'Cube data' when
 *                                    both are blank)
 * @property {(string|null)[]} datasetUnits one per `values` entry: the unit of
 *                                    its numbers (CUBE_DENSITY_UNIT for a density),
 *                                    null when unknown or not a density
 */

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
    /** @type {Float64Array | null} scratch for `number()` */
    this._one = null;
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
   * The next line that is not blank, or null at the end of the text. Past the
   * two comment lines a blank line carries nothing, so the header skips them.
   * @returns {string | null}
   */
  contentLine() {
    let line = this.line();
    while (line !== null && !line.trim()) line = this.line();
    return line;
  }

  /** Characters not yet consumed. */
  remaining() {
    return this.text.length - this.pos;
  }

  /**
   * The next whitespace-separated number, or NaN when the text is exhausted.
   * Throws on a token that is not a number.
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
   *   of the text
   */
  fill(out, offset, count) {
    const text = this.text;
    const len = text.length;
    let p = this.pos;
    let n = 0;
    while (n < count) {
      let c = 0;
      while (p < len && (c = text.charCodeAt(p)) <= 32) p++;
      if (p >= len) break;
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
 * Split a header line into its leading finite numbers; throws naming the line
 * when there are fewer than `minCount`.
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
    throw new Error(`Cube file: header is incomplete or malformed — expected ${what}, got ${line == null ? 'end of file' : `"${line.trim()}"`}`);
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
 * @param {string} description
 * @returns {string | null}
 */
export function cubeDensityUnit(description) {
  const quantity = String(description || '').toLowerCase().split(/\s+from\s+/)[0];
  if (!/densit/.test(quantity)) return null;
  if (/potential|laplacian|gradient|reduced|orbital|localiz|\belf\b|\bmo\b|norm/.test(quantity)) return null;
  return CUBE_DENSITY_UNIT;
}

/**
 * Parse the text of a Gaussian cube file.
 *
 * @param {string} content
 * @returns {CubeData}
 */
export function parseCube(content) {
  if (typeof content !== 'string') throw new Error('Cube file: expected text content');
  const scan = new CubeScanner(content);

  const title = (scan.line() ?? '').trim();
  const description = (scan.line() ?? '').trim();
  const label = `${title} ${description}`.trim();

  const countLine = headerNumbers(scan.contentLine(), 4, 'NATOMS and the origin on line 3');
  const natomsSigned = countLine[0];
  if (!Number.isInteger(natomsSigned)) throw new Error(`Cube file: atom count "${natomsSigned}" is not an integer`);
  const natoms = Math.abs(natomsSigned);
  let nval = countLine.length > 4 ? countLine[4] : 1;
  if (!Number.isInteger(nval) || nval < 1) throw new Error(`Cube file: NVAL "${nval}" must be a positive integer`);

  const grid = [];
  const rawSteps = [];
  for (let axis = 0; axis < 3; axis++) {
    const row = headerNumbers(scan.contentLine(), 4, `the axis ${axis + 1} line (N and a step vector)`);
    if (!Number.isInteger(row[0]) || row[0] === 0) {
      throw new Error(`Cube file: bad grid size "${row[0]}" on axis ${axis + 1}`);
    }
    grid.push(row[0]);
    rawSteps.push(row.slice(1, 4));
  }
  const units = grid[0] > 0 ? 'bohr' : 'angstrom';
  const unitScale = units === 'bohr' ? BOHR_TO_ANGSTROM : 1;
  for (let axis = 0; axis < 3; axis++) grid[axis] = Math.abs(grid[axis]);
  const origin = countLine.slice(1, 4).map((c) => c * unitScale);
  const voxel = rawSteps.map((step) => step.map((c) => c * unitScale));
  const lattice = voxel.map((step, axis) => step.map((c) => c * grid[axis]));

  const atoms = [];
  for (let j = 0; j < natoms; j++) {
    const row = headerNumbers(scan.contentLine(), 5, `atom line ${j + 1} of ${natoms} (Z, charge, x, y, z)`);
    atoms.push({
      atomicNumber: Math.round(row[0]),
      charge: row[1],
      position: [row[2] * unitScale, row[3] * unitScale, row[4] * unitScale],
    });
  }

  // A negative atom count announces a dataset-id record: M, then M ids. It is
  // written ten ids to a line, so it can wrap.
  let datasetIds = null;
  if (natomsSigned < 0) {
    const m = scan.number();
    if (!Number.isInteger(m) || m < 1) {
      throw new Error(`Cube file: orbital cube (NATOMS < 0) needs a dataset count after the atoms, got "${m}"`);
    }
    datasetIds = [];
    for (let i = 0; i < m; i++) {
      const id = scan.number();
      if (!Number.isInteger(id)) throw new Error(`Cube file: dataset-id record ends after ${i} of ${m} ids`);
      datasetIds.push(id);
    }
    nval = m;
  }

  const [n1, n2, n3] = grid;
  const npoints = n1 * n2 * n3;
  const total = npoints * nval;
  const truncated = (found) => new Error(`Cube file: data is truncated — expected ${total} data values `
    + `(${n1}×${n2}×${n3} points × ${nval}), found ${found}`);
  // Every value takes at least one character and a separator. A header that
  // promises more than the rest of the text can hold is refused before the
  // arrays are allocated, so a hostile grid size cannot exhaust memory.
  if (total > (scan.remaining() + 1) / 2) throw truncated(`at most ${Math.floor((scan.remaining() + 1) / 2)}`);

  const values = Array.from({ length: nval }, () => new Float32Array(npoints));

  // File order: x slowest, z fastest, the values of a point adjacent. Stored
  // order: x fastest (index = x + n1*(y + n2*z)), which is what Field.getValueAt
  // and the marching-cubes code expect. Read one record (a fixed x, y: all z and
  // all values) at a time, then scatter it.
  const n12 = n1 * n2;
  const recordLength = n3 * nval;
  const record = new Float32Array(recordLength);
  let read = 0;
  for (let x = 0; x < n1; x++) {
    for (let y = 0; y < n2; y++) {
      const got = scan.fill(record, 0, recordLength);
      read += got;
      if (got < recordLength) throw truncated(read);
      const base = x + n1 * y;
      if (nval === 1) {
        const out = values[0];
        for (let z = 0, index = base; z < n3; z++, index += n12) out[index] = record[z];
      } else {
        for (let d = 0; d < nval; d++) {
          const out = values[d];
          for (let z = 0, index = base; z < n3; z++, index += n12) out[index] = record[z * nval + d];
        }
      }
    }
  }

  return {
    label, units, origin, grid, voxel, lattice, atoms, datasetIds, values,
    ...describeDatasets(nval, { title, description, datasetIds }),
  };
}

/**
 * Label each dataset and give it a unit.
 * @param {number} count number of datasets
 * @param {{title: string, description: string, datasetIds: number[] | null}} meta
 * @returns {{datasetLabels: string[], datasetUnits: (string | null)[]}}
 */
function describeDatasets(count, { title, description, datasetIds }) {
  if (datasetIds) {
    return { datasetLabels: datasetIds.map((id) => `MO ${id}`), datasetUnits: datasetIds.map(() => null) };
  }
  if (count === 1) {
    // Many non-Gaussian writers (VMD, CP2K, Quantum ESPRESSO's pp.x, ASE) put
    // the fixed "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" on line 2 and
    // what the file holds on line 1; that line says nothing about the field.
    const informative = description && !/outer\s+loop/i.test(description) ? description : '';
    const label = informative || title || 'Cube data';
    return { datasetLabels: [label], datasetUnits: [cubeDensityUnit(label)] };
  }
  if (count === 4) {
    // A Gaussian `Gradient` cube: the density followed by its gradient. The
    // description names the gradient, so the density unit is judged from
    // whether a density is mentioned at all (and not, say, a potential).
    const text = description.toLowerCase();
    const densityUnit = /densit/.test(text) && !/potential|laplacian/.test(text) ? CUBE_DENSITY_UNIT : null;
    return { datasetLabels: [...GRADIENT_LABELS], datasetUnits: [densityUnit, null, null, null] };
  }
  const indices = Array.from({ length: count }, (_, i) => i + 1);
  return { datasetLabels: indices.map((i) => `Value ${i}`), datasetUnits: indices.map(() => null) };
}
