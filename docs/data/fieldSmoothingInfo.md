# Surface Smoothing (experimental)

Marching cubes cuts the isosurface out of the grid one cube at a time, so on a coarse grid the
surface shows facets and a faint staircase. The methods here smooth it in one of two ways:

- **Mesh methods** (Laplacian, Taubin, HC, Loop, Catmull–Clark) move or add vertices of the
  finished triangle surface.
- **Field methods** (tricubic refinement, Gaussian pre-filter) change the grid *before* it is
  meshed.

The ray/path tracers render the raw field directly, so mesh methods only change the rasterised
surface. The open rims of the surface — where it meets the cell faces or the NCI density window —
never move, so the periodic copies of a cell keep meeting seamlessly.

**Project onto isovalue** (mesh methods): after smoothing, each vertex is pulled back along the
field gradient onto the exact isovalue. The surface keeps its smooth triangulation but stays
quantitatively where the isovalue says it is. Turn it off for purely visual smoothing.

## Laplacian

Moves every vertex a fraction **λ** of the way toward the average of its neighbours, **Iterations**
times. Very effective against noise, but it shrinks the surface as it goes — use projection, or
Taubin/HC instead.

## Taubin λ|μ

Each iteration is a Laplacian step with **λ** followed by a negative (inflating) step μ, chosen
from the **pass-band** k_PB as 1/λ + 1/μ = k_PB. Small bumps are removed while the overall shape
and volume are kept (G. Taubin, 1995). Lower pass-band = stronger smoothing.

## HC-Laplacian

A Laplacian step followed by a correction that pushes each vertex back toward its original position
(**α**) and toward agreement with its neighbours' corrections (**β**). Keeps volume and sharper
features (Vollmer, Mencl & Müller, 1999).

## Loop subdivision

Splits every triangle into four per **level** and places the vertices by Loop's weights
(C. Loop, 1987), giving a finer, rounder surface. The triangle count grows 4× per level.

## Catmull–Clark

The classic quad subdivision scheme (Catmull & Clark, 1978) applied to the triangle mesh: every
triangle becomes three quads, each drawn as two triangles, so the count grows 6× per level. Loop is
its triangle-native counterpart and usually gives the more even result; Catmull–Clark is offered
for comparison.

## Tricubic refinement

Upsamples the grid by **factor** along each axis with a cubic kernel before meshing, so marching
cubes works on a finer grid:

- **Catmull-Rom** interpolates — the original grid values are reproduced exactly.
- **B-spline** approximates — it smooths slightly as it refines.

Memory and triangles grow with factor³ (factor² for triangles); very large results are capped by
lowering the factor automatically. This is also the most accurate method.

## Gaussian pre-filter

Blurs the grid with a Gaussian of width **σ** (in voxels) before meshing. Cheap and effective on
noisy data, but peaks flatten, so the same isovalue encloses a slightly different volume than on
the raw field.
