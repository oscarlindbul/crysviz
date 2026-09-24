import * as workerPool from './workerPool.js';
import { computeNci } from '../math/nci-backend-wasm.js';

/**
 * Dispatch for the NCI calculation: worker when possible, main thread when
 * not (same arrangement as waveTasks.js).
 *
 * A density grid runs to several million points and the promolecular variant
 * sums dozens of atoms per point, so on the main thread this is a visible
 * freeze. The density goes in and the two Float32 fields come back as
 * transferred ArrayBuffers; the WASM heap copies are allocated and freed inside
 * the worker.
 */

/**
 * Run `computeNci` (see math/nci-backend-wasm.js for the request fields).
 *
 * `req.values` is copied before it is transferred, so the caller's field
 * buffer is never detached.
 *
 * @param {import('../math/nci-backend-wasm.js').NciRequest} req
 * @returns {Promise<import('../math/nci-backend-wasm.js').NciResult>}
 */
export async function runNci(req) {
  if (workerPool.available()) {
    try {
      return await runInWorker(req);
    } catch (error) {
      // A worker that cannot start (module worker unsupported, wasm blocked
      // by a CSP) must not take the feature down with it.
      console.warn('NCI calculation failed in a worker; falling back to the '
        + 'main thread. The UI will block while it runs.', error);
    }
  }
  return computeNci(req);
}

/** @param {import('../math/nci-backend-wasm.js').NciRequest} req */
async function runInWorker(req) {
  const values = req.values ? new Float32Array(req.values) : undefined;
  const payload = {
    kind: req.kind,
    nx: req.nx,
    ny: req.ny,
    nz: req.nz,
    voxel: req.voxel.map((row) => Array.from(row)),
    periodic: !!req.periodic,
    values,
    valueScale: req.valueScale,
    atoms: req.atoms
      ? { Z: Array.from(req.atoms.Z), frac: req.atoms.frac.map((f) => Array.from(f)) }
      : undefined,
    params: req.params ? { ...req.params } : undefined,
  };
  return workerPool.run('nci', payload, values ? [values.buffer] : []);
}
