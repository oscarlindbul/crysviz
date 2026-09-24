#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

/*
 * Credits
 * -------
 * The NCI method is due to its original developers:
 *
 *   E. R. Johnson, S. Keinan, P. Mori-Sanchez, J. Contreras-Garcia,
 *   A. J. Cohen and W. Yang, "Revealing Noncovalent Interactions",
 *   J. Am. Chem. Soc. 132, 6498-6506 (2010), doi:10.1021/ja100936w.
 *
 *   J. Contreras-Garcia, E. R. Johnson, S. Keinan, R. Chaudret,
 *   J.-P. Piquemal, D. N. Beratan and W. Yang, "NCIPLOT: A Program for
 *   Plotting Noncovalent Interaction Regions", J. Chem. Theory Comput. 7,
 *   625-632 (2011), doi:10.1021/ct100641a.
 *
 * The promolecular atomic-density parameters (coef/zeta tables below) are
 * NCIPLOT's. This file is a port of Jmol's NCI implementation,
 * org.jmol.quantum.NciCalculation by Bob Hanson (Jmol, https://jmol.sourceforge.net,
 * LGPL-2.1-or-later), including its discrete "SCF" finite-difference mode.
 * Henry Rzepa's cub2nci page (https://www.ch.ic.ac.uk/rzepa/cub2nci/) was the
 * reference for how Jmol is used on Gaussian cubes.
 */

/*
 * NCI (non-covalent interaction) analysis for WASM.
 *
 * A port of Jmol's `org.jmol.quantum.NciCalculation` (reference copy in
 * NCI_analysis/cub2nci/jmol_src/), which in turn follows NCIPLOT 1.0
 * (Johnson et al., JACS 132, 6498 (2010); Contreras-Garcia et al., JCTC 7,
 * 625 (2011)). Two scalar fields are produced on a grid:
 *
 *   s          = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3))   reduced gradient
 *   sl2rho     = sign(lambda2) * rho                            colouring
 *
 * with lambda2 the middle eigenvalue of the density Hessian. rho is in atomic
 * units (e/bohr^3) and every length is in bohr.
 *
 * Two sources for rho and its derivatives:
 *
 *   nci_scf           - finite differences on a density grid (Jmol's
 *                       "discrete SCF" mode, calcPlane/getPlaneValue).
 *   nci_promolecular  - analytic sum of spherical three-exponential atomic
 *                       densities (NCIPLOT's H-Ar tables, Jmol processAtoms).
 *
 * Differences from Jmol, all deliberate:
 *
 * - Non-orthogonal grids. Jmol differentiates along each grid axis with that
 *   axis' step length, which is only right for an orthogonal grid. Here the
 *   index-space gradient g and Hessian H_idx are transformed to Cartesian:
 *   with V the voxel matrix (rows = voxel vectors, bohr), r = u V for index
 *   coordinates u, so g = V grad and H_idx = V H V^T, i.e.
 *       grad = V^-1 g,   H = V^-1 H_idx V^-T.
 *   (Equivalently, with M = V^T whose columns are the voxel vectors,
 *   grad = M^-T g and H = M^-T H_idx M^-1.)
 * - Periodic grids wrap their neighbours instead of dropping the boundary.
 * - Excluded points get s = sCap and sign(lambda2)rho = 0 instead of NaN /
 *   Jmol's NO_VALUE, because marching cubes and the field statistics do not
 *   tolerate NaN. s is also clipped at sCap.
 * - lambda2 comes from a closed-form symmetric 3x3 eigen-solver.
 * - Promolecular: Jmol leaves the atom loop as soon as the partial sum is
 *   below rhoMin, which is a bug (later atoms still add density). Here a
 *   point is only excluded for rho < rhoMin after all atoms are summed; the
 *   early exit on rho > rhoPlot is kept (the sum only grows).
 *
 * Grid layout matches model/Field.js: x fastest, idx = i + nx*(j + ny*k), and
 * grid point (i,j,k) sits at r = i*v0 + j*v1 + k*v2 relative to the field
 * origin. Matrices are flat row-major length 9, rows are the voxel vectors.
 */

/* c = 1 / (2 (3 pi^2)^(1/3)) */
#define NCI_C 0.16162045967399548
#define NCI_PI 3.14159265358979323846

enum {
  NCI_ERR_BAD_INPUT = -1,
  NCI_ERR_SINGULAR = -2,
  NCI_ERR_ALLOC = -3
};

/* ------------------------------------------------------------------ *
 * Promolecular tables (NCIPLOT 1.0 via Jmol). Index = Z, 1..18.
 * rho_atom(r) = sum_n coef_n * exp(-r / zeta_n), r in bohr.
 * ------------------------------------------------------------------ */

static const double coef1[19] = {
  0, 0.2815, 2.437,
  11.84, 31.34, 67.82, 120.2, 190.9, 289.5, 406.3, 561.3,
  760.8, 1016., 1319., 1658., 2042., 2501., 3024., 3625.
};
static const double coef2[19] = {
  0, 0., 0.,
  0.06332, 0.3694, 0.8527, 1.172, 2.247, 2.879, 3.049, 6.984,
  22.42, 37.17, 57.95, 87.16, 115.7, 158.0, 205.5, 260.0
};
static const double coef3[19] = {
  0, 0., 0.,
  0., 0., 0., 0., 0., 0., 0., 0.,
  0.06358, 0.3331, 0.8878, 0.7888, 1.465, 2.170, 3.369, 5.211
};
static const double zeta1[19] = {
  0, 0.5288, 0.3379,
  0.1912, 0.1390, 0.1059, 0.0884, 0.0767, 0.0669, 0.0608, 0.0549,
  0.0496, 0.0449, 0.0411, 0.0382, 0.0358, 0.0335, 0.0315, 0.0296
};
static const double zeta2[19] = {
  0, 1., 1.,
  0.9992, 0.6945, 0.5300, 0.5480, 0.4532, 0.3974, 0.3994, 0.3447,
  0.2511, 0.2150, 0.1874, 0.1654, 0.1509, 0.1369, 0.1259, 0.1168
};
static const double zeta3[19] = {
  0, 1., 1.,
  1., 1., 1., 1., 1., 1., 1., 1.,
  1.0236, 0.7753, 0.5962, 0.6995, 0.5851, 0.5149, 0.4974, 0.4412
};
/* Jmol's per-element cutoff (bohr) where the atomic density drops to ~0.001.
 * Applied as a box test on each Cartesian component, exactly as Jmol does. */
static const double dMax[19] = {
  0, 2.982502423, 2.635120936,
  4.144887422, 4.105800759, 3.576656363, 3.872424373, 3.497503547,
  3.165369971, 3.204214082, 3.051069564,
  4.251312809, 4.503309314, 4.047465141, 4.666024968, 4.265151411,
  3.955710076, 4.040067606, 3.776022242
};

/* ------------------------------------------------------------------ *
 * Small linear algebra
 * ------------------------------------------------------------------ */

static int wrap_index(int g, int n) {
  int v = g % n;
  if (v < 0) v += n;
  return v;
}

/* Inverse of a row-major 3x3. Returns 0 when singular. */
static int invert3(const double *m, double *out) {
  double c00 = m[4] * m[8] - m[5] * m[7];
  double c01 = m[5] * m[6] - m[3] * m[8];
  double c02 = m[3] * m[7] - m[4] * m[6];
  double det = m[0] * c00 + m[1] * c01 + m[2] * c02;
  double scale = 0.0;
  int i;
  for (i = 0; i < 9; i++) scale = fmax(scale, fabs(m[i]));
  if (!(scale > 0.0) || !isfinite(det) || fabs(det) <= 1e-14 * scale * scale * scale) return 0;
  det = 1.0 / det;
  out[0] = c00 * det;
  out[1] = (m[2] * m[7] - m[1] * m[8]) * det;
  out[2] = (m[1] * m[5] - m[2] * m[4]) * det;
  out[3] = c01 * det;
  out[4] = (m[0] * m[8] - m[2] * m[6]) * det;
  out[5] = (m[2] * m[3] - m[0] * m[5]) * det;
  out[6] = c02 * det;
  out[7] = (m[1] * m[6] - m[0] * m[7]) * det;
  out[8] = (m[0] * m[4] - m[1] * m[3]) * det;
  return 1;
}

/*
 * Eigenvalues of the symmetric 3x3 matrix
 *     | a00 a01 a02 |
 *     | a01 a11 a12 |
 *     | a02 a12 a22 |
 * in ascending order, by the closed-form trigonometric method (Smith 1961,
 * "Eigenvalues of a symmetric 3x3 matrix", CACM 4:168).
 *
 * The matrix is shifted by its mean eigenvalue q and scaled by p, so
 * B = (A - qI)/p has eigenvalues 2cos(phi + 2 pi k/3) with det(B) = 2cos(3phi).
 * Robustness details:
 * - a diagonal matrix is sorted directly (exact);
 * - a (numerically) multiple of the identity returns q three times;
 * - det(B)/2 is clamped to [-1, 1] before acos, which is where rounding
 *   pushes it for (near-)degenerate pairs;
 * - each eigenvalue is evaluated from its own cosine rather than as
 *   3q - l1 - l3, so the middle one does not suffer the cancellation of
 *   subtracting two large eigenvalues.
 */
static void sym3_eigvals(double a00, double a11, double a22,
                         double a01, double a02, double a12, double *out) {
  double p1 = a01 * a01 + a02 * a02 + a12 * a12;
  double q, b00, b11, b22, p2, p, det, r, phi, t;

  if (p1 == 0.0) {
    double e0 = a00, e1 = a11, e2 = a22;
    if (e0 > e1) { t = e0; e0 = e1; e1 = t; }
    if (e1 > e2) { t = e1; e1 = e2; e2 = t; }
    if (e0 > e1) { t = e0; e0 = e1; e1 = t; }
    out[0] = e0; out[1] = e1; out[2] = e2;
    return;
  }

  q = (a00 + a11 + a22) / 3.0;
  b00 = a00 - q;
  b11 = a11 - q;
  b22 = a22 - q;
  p2 = b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * p1;
  p = sqrt(p2 / 6.0);
  if (!(p > 0.0)) {
    out[0] = out[1] = out[2] = q;
    return;
  }
  b00 /= p; b11 /= p; b22 /= p;
  {
    double c01 = a01 / p, c02 = a02 / p, c12 = a12 / p;
    det = b00 * (b11 * b22 - c12 * c12)
        - c01 * (c01 * b22 - c12 * c02)
        + c02 * (c01 * c12 - b11 * c02);
  }
  r = 0.5 * det;
  if (r <= -1.0) phi = NCI_PI / 3.0;
  else if (r >= 1.0) phi = 0.0;
  else phi = acos(r) / 3.0;

  /* phi in [0, pi/3]: cos(phi) is the largest, cos(phi + 2pi/3) the
   * smallest and cos(phi + 4pi/3) the middle of the three. */
  out[2] = q + 2.0 * p * cos(phi);
  out[0] = q + 2.0 * p * cos(phi + 2.0 * NCI_PI / 3.0);
  out[1] = q + 2.0 * p * cos(phi + 4.0 * NCI_PI / 3.0);
  /* Guard the ordering against rounding at the ends of the phi range. */
  if (out[1] < out[0]) { t = out[0]; out[0] = out[1]; out[1] = t; }
  if (out[1] > out[2]) { t = out[2]; out[2] = out[1]; out[1] = t; }
}

/* Exported for the unit tests: a6 = [a00, a11, a22, a01, a02, a12]. */
EMSCRIPTEN_KEEPALIVE
void nci_sym3_eigvals(const double *a6, double *out3) {
  sym3_eigvals(a6[0], a6[1], a6[2], a6[3], a6[4], a6[5], out3);
}

/*
 * Write one point's outputs from rho, the Cartesian gradient and the Cartesian
 * Hessian (h = [xx, yy, zz, xy, xz, yz]). Returns 1 if the point is included.
 */
static int finish_point(double rho, const double *g, const double *h,
                        double rhoMin, double rhoPlot, double sCap,
                        float *outS, float *outSl2rho, size_t idx) {
  double grad, s, ev[3];
  if (!(rho > 0.0) || rho < rhoMin || rho > rhoPlot) {
    outS[idx] = (float)sCap;
    outSl2rho[idx] = 0.0f;
    return 0;
  }
  grad = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
  s = NCI_C * grad * pow(rho, -4.0 / 3.0);
  if (!(s < sCap)) s = sCap; /* also catches NaN */
  sym3_eigvals(h[0], h[1], h[2], h[3], h[4], h[5], ev);
  outS[idx] = (float)s;
  outSl2rho[idx] = (float)(ev[1] < 0.0 ? -rho : rho);
  return 1;
}

/* ------------------------------------------------------------------ *
 * SCF (finite differences on a density grid)
 * ------------------------------------------------------------------ */

/*
 * rho      - density grid, x fastest; rho = |value * valueScale| (Jmol takes
 *            the absolute value too)
 * voxel    - 9 doubles, rows = voxel vectors in bohr
 * periodic - nonzero: neighbours wrap (index n is index 0). Zero: the
 *            boundary layer (any index 0 or n-1) is excluded, like Jmol.
 *
 * Index-space central differences, as Jmol's calcPlane/getPlaneValue:
 *   g_a     = (f(+a) - f(-a)) / 2
 *   H_aa    = f(+a) - 2 f + f(-a)
 *   H_ab    = (f(+a+b) - f(+a-b) - f(-a+b) + f(-a-b)) / 4
 * then transformed to Cartesian with V^-1 (see the header comment).
 *
 * Returns the number of included points, or a negative NCI_ERR_* code.
 */
EMSCRIPTEN_KEEPALIVE
int nci_scf(const float *values, int nx, int ny, int nz,
            const double *voxel, int periodic, double valueScale,
            double rhoMin, double rhoPlot, double sCap,
            float *outS, float *outSl2rho) {
  double B[9];
  int *xp, *xm, *yp, *ym, *zp, *zm;
  int i, j, k, count = 0;
  const size_t sx = 1, sy = (size_t)nx, sz = (size_t)nx * (size_t)ny;

  if (!values || !voxel || !outS || !outSl2rho || nx < 1 || ny < 1 || nz < 1)
    return NCI_ERR_BAD_INPUT;
  if (!invert3(voxel, B)) return NCI_ERR_SINGULAR;

  xp = (int *)malloc(sizeof(int) * (size_t)(2 * (nx + ny + nz)));
  if (!xp) return NCI_ERR_ALLOC;
  xm = xp + nx;
  yp = xm + nx;
  ym = yp + ny;
  zp = ym + ny;
  zm = zp + nz;
  for (i = 0; i < nx; i++) { xp[i] = wrap_index(i + 1, nx); xm[i] = wrap_index(i - 1, nx); }
  for (j = 0; j < ny; j++) { yp[j] = wrap_index(j + 1, ny); ym[j] = wrap_index(j - 1, ny); }
  for (k = 0; k < nz; k++) { zp[k] = wrap_index(k + 1, nz); zm[k] = wrap_index(k - 1, nz); }

#define RHO(ii, jj, kk) fabs((double)values[(size_t)(ii) * sx + (size_t)(jj) * sy + (size_t)(kk) * sz] * valueScale)

  for (k = 0; k < nz; k++) {
    int kb = !periodic && (k == 0 || k == nz - 1);
    int k1 = zp[k], k0 = zm[k];
    for (j = 0; j < ny; j++) {
      int jb = kb || (!periodic && (j == 0 || j == ny - 1));
      int j1 = yp[j], j0 = ym[j];
      size_t row = (size_t)j * sy + (size_t)k * sz;
      for (i = 0; i < nx; i++) {
        size_t idx = row + (size_t)i;
        double rho = fabs((double)values[idx] * valueScale);
        int i1, i0;
        double gi[3], hi[6], g[3], h[6], T[9];
        int a, b;

        if (jb || (!periodic && (i == 0 || i == nx - 1))
            || !(rho > 0.0) || rho < rhoMin || rho > rhoPlot) {
          outS[idx] = (float)sCap;
          outSl2rho[idx] = 0.0f;
          continue;
        }
        i1 = xp[i];
        i0 = xm[i];
        {
          double fxp = RHO(i1, j, k), fxm = RHO(i0, j, k);
          double fyp = RHO(i, j1, k), fym = RHO(i, j0, k);
          double fzp = RHO(i, j, k1), fzm = RHO(i, j, k0);
          gi[0] = 0.5 * (fxp - fxm);
          gi[1] = 0.5 * (fyp - fym);
          gi[2] = 0.5 * (fzp - fzm);
          hi[0] = fxp - 2.0 * rho + fxm;
          hi[1] = fyp - 2.0 * rho + fym;
          hi[2] = fzp - 2.0 * rho + fzm;
          hi[3] = 0.25 * (RHO(i1, j1, k) - RHO(i1, j0, k) - RHO(i0, j1, k) + RHO(i0, j0, k));
          hi[4] = 0.25 * (RHO(i1, j, k1) - RHO(i1, j, k0) - RHO(i0, j, k1) + RHO(i0, j, k0));
          hi[5] = 0.25 * (RHO(i, j1, k1) - RHO(i, j1, k0) - RHO(i, j0, k1) + RHO(i, j0, k0));
        }
        /* grad = B g */
        for (a = 0; a < 3; a++)
          g[a] = B[3 * a] * gi[0] + B[3 * a + 1] * gi[1] + B[3 * a + 2] * gi[2];
        /* H = B H_idx B^T. T = B H_idx first (H_idx symmetric, packed). */
        {
          double Hf[9];
          Hf[0] = hi[0]; Hf[4] = hi[1]; Hf[8] = hi[2];
          Hf[1] = Hf[3] = hi[3];
          Hf[2] = Hf[6] = hi[4];
          Hf[5] = Hf[7] = hi[5];
          for (a = 0; a < 3; a++)
            for (b = 0; b < 3; b++)
              T[3 * a + b] = B[3 * a] * Hf[b] + B[3 * a + 1] * Hf[3 + b] + B[3 * a + 2] * Hf[6 + b];
        }
#define HB(a_, b_) (T[3 * (a_)] * B[3 * (b_)] + T[3 * (a_) + 1] * B[3 * (b_) + 1] + T[3 * (a_) + 2] * B[3 * (b_) + 2])
        h[0] = HB(0, 0);
        h[1] = HB(1, 1);
        h[2] = HB(2, 2);
        h[3] = HB(0, 1);
        h[4] = HB(0, 2);
        h[5] = HB(1, 2);
#undef HB
        count += finish_point(rho, g, h, rhoMin, rhoPlot, sCap, outS, outSl2rho, idx);
      }
    }
  }
#undef RHO

  free(xp);
  return count;
}

/* ------------------------------------------------------------------ *
 * Promolecular (analytic atomic densities)
 * ------------------------------------------------------------------ */

typedef struct {
  double x, y, z;
  double dmax;          /* box half-width (Jmol's dMax[Z]) */
  double reach;         /* sqrt(3)*dmax: farthest point the box test passes */
  int nterm;
  double c[3], iz[3];   /* coefficients and 1/zeta for the nonzero terms */
} nci_atom;

/*
 * Shared driver. For every grid point, sums rho, grad and Hessian over the
 * atoms whose dMax box contains the point (Jmol's box test on |dx|,|dy|,|dz|).
 *
 * Work is organised per x-row: the slab (fixed k) keeps only atoms within
 * `reach` of the slab plane, and each row then solves for the i-range where
 * the atom can pass the box test, so every exp() evaluated is one Jmol would
 * evaluate too. In NCI mode a point whose partial sum already exceeds rhoPlot
 * skips the remaining atoms (the sum only grows, so it stays excluded).
 *
 * outRho != NULL selects density-only mode: the full rho is written and no
 * early exit is taken.
 */
static int promolecular_driver(const int *Z, const double *xyz, int nAtoms,
                               int nx, int ny, int nz, const double *voxel,
                               double rhoMin, double rhoPlot, double sCap,
                               float *outS, float *outSl2rho, float *outRho) {
  const double *v0 = voxel, *v1 = voxel + 3, *v2 = voxel + 6;
  nci_atom *atoms = NULL;
  int *slab = NULL, *rowList = NULL;
  int *rowLo = NULL, *rowHi = NULL;
  double *acc = NULL;
  double nrm[3], nlen, v0sq;
  int nUse = 0, a, i, j, k, count = 0;
  const int densityOnly = outRho != NULL;

  if (!voxel || nx < 1 || ny < 1 || nz < 1 || nAtoms < 0 || (nAtoms > 0 && (!Z || !xyz)))
    return NCI_ERR_BAD_INPUT;
  if (!densityOnly && (!outS || !outSl2rho)) return NCI_ERR_BAD_INPUT;
  {
    double B[9];
    if (!invert3(voxel, B)) return NCI_ERR_SINGULAR;
  }

  /* Unit normal of the (v0, v1) plane: a slab is the plane through k*v2. */
  nrm[0] = v0[1] * v1[2] - v0[2] * v1[1];
  nrm[1] = v0[2] * v1[0] - v0[0] * v1[2];
  nrm[2] = v0[0] * v1[1] - v0[1] * v1[0];
  nlen = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
  nrm[0] /= nlen; nrm[1] /= nlen; nrm[2] /= nlen;
  v0sq = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];

  atoms = (nci_atom *)malloc(sizeof(nci_atom) * (size_t)(nAtoms > 0 ? nAtoms : 1));
  slab = (int *)malloc(sizeof(int) * (size_t)(nAtoms > 0 ? nAtoms : 1));
  rowList = (int *)malloc(sizeof(int) * (size_t)(nAtoms > 0 ? nAtoms : 1));
  rowLo = (int *)malloc(sizeof(int) * (size_t)(nAtoms > 0 ? nAtoms : 1));
  rowHi = (int *)malloc(sizeof(int) * (size_t)(nAtoms > 0 ? nAtoms : 1));
  /* per-row accumulators: rho, g[3], h[6] */
  acc = (double *)malloc(sizeof(double) * 10 * (size_t)nx);
  if (!atoms || !slab || !rowList || !rowLo || !rowHi || !acc) {
    count = NCI_ERR_ALLOC;
    goto done;
  }

  for (a = 0; a < nAtoms; a++) {
    int z = Z[a];
    nci_atom *at = &atoms[nUse];
    if (z < 1) continue;
    if (z > 18) z = 18;
    at->x = xyz[3 * a];
    at->y = xyz[3 * a + 1];
    at->z = xyz[3 * a + 2];
    at->dmax = dMax[z];
    at->reach = sqrt(3.0) * dMax[z] * (1.0 + 1e-12);
    at->nterm = 0;
    if (coef1[z] != 0.0) { at->c[at->nterm] = coef1[z]; at->iz[at->nterm++] = 1.0 / zeta1[z]; }
    if (coef2[z] != 0.0) { at->c[at->nterm] = coef2[z]; at->iz[at->nterm++] = 1.0 / zeta2[z]; }
    if (coef3[z] != 0.0) { at->c[at->nterm] = coef3[z]; at->iz[at->nterm++] = 1.0 / zeta3[z]; }
    nUse++;
  }

  for (k = 0; k < nz; k++) {
    int nSlab = 0;
    double p0[3];
    p0[0] = k * v2[0]; p0[1] = k * v2[1]; p0[2] = k * v2[2];
    for (a = 0; a < nUse; a++) {
      double d = (atoms[a].x - p0[0]) * nrm[0] + (atoms[a].y - p0[1]) * nrm[1]
               + (atoms[a].z - p0[2]) * nrm[2];
      if (fabs(d) <= atoms[a].reach) slab[nSlab++] = a;
    }

    for (j = 0; j < ny; j++) {
      double r0[3];
      int nRow = 0;
      size_t row = (size_t)nx * ((size_t)j + (size_t)ny * (size_t)k);
      r0[0] = p0[0] + j * v1[0];
      r0[1] = p0[1] + j * v1[1];
      r0[2] = p0[2] + j * v1[2];

      /* Atoms whose reach sphere meets the row line r0 + t*v0, and the
       * t-range where it does: |r0 - A + t v0|^2 <= reach^2. */
      for (a = 0; a < nSlab; a++) {
        const nci_atom *at = &atoms[slab[a]];
        double dx = r0[0] - at->x, dy = r0[1] - at->y, dz = r0[2] - at->z;
        double bq = dx * v0[0] + dy * v0[1] + dz * v0[2];
        double cq = dx * dx + dy * dy + dz * dz - at->reach * at->reach;
        double disc = bq * bq - v0sq * cq;
        double sq, t0, t1;
        int lo, hi;
        if (disc < 0.0) continue;
        sq = sqrt(disc);
        t0 = (-bq - sq) / v0sq;
        t1 = (-bq + sq) / v0sq;
        if (t1 < 0.0 || t0 > nx - 1) continue;
        lo = t0 <= 0.0 ? 0 : (int)ceil(t0);
        hi = t1 >= nx - 1 ? nx - 1 : (int)floor(t1);
        if (lo > hi) continue;
        rowList[nRow] = slab[a];
        rowLo[nRow] = lo;
        rowHi[nRow] = hi;
        nRow++;
      }

      for (i = 0; i < 10 * nx; i++) acc[i] = 0.0;

      for (a = 0; a < nRow; a++) {
        const nci_atom *at = &atoms[rowList[a]];
        const double dm = at->dmax;
        int n;
        for (i = rowLo[a]; i <= rowHi[a]; i++) {
          double *ac = acc + 10 * (size_t)i;
          double x, y, z, r, rhoA = 0.0, f1 = 0.0, f2 = 0.0;
          if (!densityOnly && ac[0] > rhoPlot) continue; /* saturated */
          x = r0[0] + i * v0[0] - at->x;
          y = r0[1] + i * v0[1] - at->y;
          z = r0[2] + i * v0[2] - at->z;
          if (fabs(x) > dm || fabs(y) > dm || fabs(z) > dm) continue;
          r = sqrt(x * x + y * y + z * z);
          for (n = 0; n < at->nterm; n++) {
            double ce = at->c[n] * exp(-r * at->iz[n]);
            rhoA += ce;
            f1 += ce * at->iz[n];               /* -d rho / dr    */
            f2 += ce * at->iz[n] * at->iz[n];   /*  d2 rho / dr2  */
          }
          ac[0] += rhoA;
          if (densityOnly) continue;
          /* At the nucleus the gradient vanishes by symmetry and the cusp
           * Hessian is undefined; rho there is far above any rhoPlot. */
          if (r < 1e-12) continue;
          {
            /* grad rho = rho'(r) rhat = -f1 * x / r
             * H = rho'' rhat rhat^T + (rho'/r)(I - rhat rhat^T)
             *   = (f2 + f1/r) rhat rhat^T - (f1/r) I          (Jmol's fr2, fac1r) */
            double fac1r = f1 / r;
            double fr2 = (fac1r + f2) / (r * r);
            ac[1] -= fac1r * x;
            ac[2] -= fac1r * y;
            ac[3] -= fac1r * z;
            ac[4] += fr2 * x * x - fac1r;
            ac[5] += fr2 * y * y - fac1r;
            ac[6] += fr2 * z * z - fac1r;
            ac[7] += fr2 * x * y;
            ac[8] += fr2 * x * z;
            ac[9] += fr2 * y * z;
          }
        }
      }

      if (densityOnly) {
        for (i = 0; i < nx; i++) outRho[row + (size_t)i] = (float)acc[10 * (size_t)i];
        count += nx;
      } else {
        for (i = 0; i < nx; i++) {
          const double *ac = acc + 10 * (size_t)i;
          count += finish_point(ac[0], ac + 1, ac + 4, rhoMin, rhoPlot, sCap,
                                outS, outSl2rho, row + (size_t)i);
        }
      }
    }
  }

done:
  free(atoms);
  free(slab);
  free(rowList);
  free(rowLo);
  free(rowHi);
  free(acc);
  return count;
}

/*
 * Promolecular NCI on the grid r(i,j,k) = i*v0 + j*v1 + k*v2 (bohr).
 *
 * Z        - atomic numbers; Z < 1 is skipped, Z > 18 is treated as 18 (Ar),
 *            as in Jmol (the caller reports the clamping)
 * xyz      - 3*nAtoms Cartesian bohr in the grid frame (relative to the grid
 *            origin). Periodic images are the caller's job.
 *
 * Returns the number of included points, or a negative NCI_ERR_* code.
 */
EMSCRIPTEN_KEEPALIVE
int nci_promolecular(const int *Z, const double *xyz, int nAtoms,
                     int nx, int ny, int nz, const double *voxel,
                     double rhoMin, double rhoPlot, double sCap,
                     float *outS, float *outSl2rho) {
  return promolecular_driver(Z, xyz, nAtoms, nx, ny, nz, voxel,
                             rhoMin, rhoPlot, sCap, outS, outSl2rho, NULL);
}

/*
 * The promolecular density itself (e/bohr^3) on the same grid, with the same
 * dMax box truncation. Used by the tests (and available for a "promolecular
 * density" field). Returns the number of points written or a negative code.
 */
EMSCRIPTEN_KEEPALIVE
int nci_promolecular_density(const int *Z, const double *xyz, int nAtoms,
                             int nx, int ny, int nz, const double *voxel,
                             float *outRho) {
  if (!outRho) return NCI_ERR_BAD_INPUT;
  return promolecular_driver(Z, xyz, nAtoms, nx, ny, nz, voxel,
                             0.0, 0.0, 0.0, NULL, NULL, outRho);
}
