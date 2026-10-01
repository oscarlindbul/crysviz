# Non-covalent interactions (NCI)

The NCI method (Johnson et al., *J. Am. Chem. Soc.* **2010**, 132, 6498) finds weak interactions —
hydrogen bonds, van der Waals contacts, steric clashes — directly from the electron density ρ. Where
two fragments interact without a covalent bond, the density is low *and* flat: its gradient
nearly vanishes. NCI looks for exactly those regions.

## The two fields it creates

Both buttons add two fields under **Derived** in the list above, on the same grid as the selected
field, and select the s field:

- **Reduced density gradient s** — `s = |∇ρ| / (2 (3π²)^(1/3) ρ^(4/3))`, a dimensionless measure of
  how far the density is from locally uniform. Low s at low ρ marks an interaction region.
- **sign(λ₂)ρ** — the density multiplied by the sign of λ₂, the middle eigenvalue of the density's
  Hessian. λ₂ < 0 where the density is concentrated between two atoms (a bond, however weak),
  λ₂ > 0 where it is depleted (a ring or a steric clash). In e/bohr³.

## Reading the plot

The s field is shown as its **s = 0.5 isosurface, coloured by sign(λ₂)ρ** over −0.04 … 0.04 a.u.:

- **Blue** (sign(λ₂)ρ strongly negative) — strong, attractive: hydrogen bonds, halogen bonds.
- **Green to yellow** (close to 0) — weak, delocalized van der Waals contacts.
- **Red** (positive) — repulsive: steric clashes, ring and cage centres.

The colour saturates at the ends of the range. The isovalue, the colouring field and its range can
be changed like any other field's.

## The two buttons

- **Create SCF-NCI field** differentiates the **selected field** by finite differences on its grid.
  It must be a **total electron density** (e.g. a Gaussian `density` cube). On anything else — a
  magnetization, an orbital, ELF — the result has no meaning.
- **Create Promolecular NCI field** ignores the field's values and builds a density from the atom
  positions alone: a sum of spherical free-atom densities (the three-exponential fits of NCIPLOT),
  evaluated on the selected field's grid so the two variants line up. It is exact for H–Ar;
  heavier elements are approximated by the Ar density, and the status line says which ones were.
  Useful for large systems, or when no density was computed.

## Density window

Only points with ρ_min < ρ < ρ_plot are kept, with ρ_min = 1e-5 a.u. and ρ_plot = 0.05 a.u.
for both buttons (Jmol uses 0.07 for promolecular densities, which leaves wide sheets around the
bonds). s is also small near nuclei and inside covalent bonds, where the
density is high; the upper cut removes those, and the lower cut removes the vacuum, where s is
dominated by numerical noise. Points outside the window get s = 2 (no surface) and
sign(λ₂)ρ = 0. The status line reports how many grid points fell inside it.

## Units

The calculation works in atomic units (e/bohr³). Gaussian cube densities already are; a CHGCAR's
density (e/Å³ after CrysViz's volume scaling) is converted automatically. For a field of unknown
unit the values are used as they are, and the status line warns about it.

## Caveats

- A **CHGCAR** holds the PAW **valence pseudo-density**, which is smooth near the cores but not the
  real density. The sum of **AECCAR0 + AECCAR2** (use *Combine fields*) is much closer to the
  all-electron density and is the better input for SCF-NCI.
- λ₂ comes from **second derivatives by finite differences**, so the grid needs to be fine:
  a spacing of about **0.1 Å or less** is recommended. Coarse grids give noisy, broken surfaces.
- A periodic grid (a CHGCAR, or a cube loaded as **Periodic**, the default) wraps around: the
  outermost layer of grid points takes its neighbours from the opposite face, and the promolecular
  density includes the atoms' periodic images. A cube loaded as **Not periodic** (a finite block,
  offered when its atoms lie outside the grid) has no neighbour beyond its outermost layer, so that
  layer is left out, and only the atoms themselves contribute.
- Intra/intermolecular filtering is not available; all interactions are shown.

## Credits

The algorithm follows Jmol's NCI implementation (Bob Hanson), including the promolecular atomic
densities from NCIPLOT (Contreras-García et al., *J. Chem. Theory Comput.* **2011**, 7, 625), and
Henry Rzepa's cub2nci page, which showed how to get NCI plots from Gaussian cubes in the browser.
