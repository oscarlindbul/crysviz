/* tslint:disable */
/* eslint-disable */

/**
 * Type stub for the Emscripten output of `nci_backend.c`.
 *
 * Same arrangement as `wave_backend.d.ts`: the module is built with
 * `ENVIRONMENT=web,worker,node` (the node branch is what the unit tests in
 * tools/unittest/ use), and the generated glue references `node:module`,
 * `node:fs` and `process`, which this project has no types for. Declaring the
 * shape here makes `tsc` prefer this file over the generated .js.
 *
 * Regenerating the module does not touch this file, so keep the signatures in
 * sync with the EMSCRIPTEN_KEEPALIVE exports in nci_backend.c. All pointer
 * arguments are byte offsets into the module heap; matrices are flat row-major
 * 3x3 with the voxel vectors (bohr) as rows.
 */

/** The instantiated module: Emscripten runtime helpers plus the C exports. */
export interface NciBackendModule {
  _malloc(size: number): number;
  _free(pointer: number): void;

  /**
   * Finite-difference NCI on a density grid (float32, x fastest). Writes the
   * reduced gradient s and sign(lambda2)*rho per point; excluded points get
   * s = sCap and 0. Returns the number of included points, or a negative
   * status (-1 bad input, -2 singular voxel matrix, -3 allocation failure).
   */
  _nci_scf(
    values: number, nx: number, ny: number, nz: number,
    voxelBohr: number, periodic: number, valueScale: number,
    rhoMin: number, rhoPlot: number, sCap: number,
    outS: number, outSl2rho: number,
  ): number;

  /**
   * Promolecular NCI: `Z` is int32[nAtoms], `xyzBohr` float64[3*nAtoms] in the
   * grid frame (relative to the grid origin). Same outputs and return value as
   * `_nci_scf`.
   */
  _nci_promolecular(
    Z: number, xyzBohr: number, nAtoms: number,
    nx: number, ny: number, nz: number, voxelBohr: number,
    rhoMin: number, rhoPlot: number, sCap: number,
    outS: number, outSl2rho: number,
  ): number;

  /** The promolecular density (e/bohr^3) on the same grid, float32 out. */
  _nci_promolecular_density(
    Z: number, xyzBohr: number, nAtoms: number,
    nx: number, ny: number, nz: number, voxelBohr: number,
    outRho: number,
  ): number;

  /** Ascending eigenvalues of [a00, a11, a22, a01, a02, a12] (float64 in/out). */
  _nci_sym3_eigvals(a6: number, out3: number): void;

  // Heap views. Re-read after every allocation on the JS side, because
  // ALLOW_MEMORY_GROWTH can detach the underlying ArrayBuffer.
  HEAPF64: Float64Array;
  HEAPF32: Float32Array;
  HEAP32: Int32Array;
}

/** MODULARIZE=1 EXPORT_ES6=1 factory. */
declare function createNciBackend(moduleOverrides?: Record<string, any>): Promise<NciBackendModule>;

export default createNciBackend;
