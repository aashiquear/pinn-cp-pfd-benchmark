# Workflow: perforated polycrystal plate for the PFD and CP-PFD stages

How a plate with a stress concentrator is meshed and handed to PINNACLE in
Part II of the benchmark. Companion documents: `09_pfd_benchmark_plan.md`
(the plan and the reasoning), `workflow_polycrystal_plate.md` (Part I:
grain growth, labelling, orientations, which are reused unchanged).

## Stage 1 -- mesh (hexaform `plate`)

The mesh is NOT a carved grid: deleting background elements inside the
hole leaves a staircase edge whose steps are spurious stress
concentrations exactly where damage must initiate. hexaform's plate
generator (added 2026-10-13, branch `feat/plate-hole`) builds a uniform
background of cell size h with one HEX8 through the thickness and a
block-structured O-grid ring around the hole:

* ring nodes lie exactly on the circle or rotated ellipse;
* the ring meets the background on a grid-aligned rectangle (half-extents
  a whole number of cells, by default the hole's bounding box plus two
  cells), so no transition element exists and the background outside the
  rectangle is identical to the unperforated plate of the same h;
* radial layers are graded so that the innermost spacing equals the
  hole-edge arc length (near-square elements at the hole) and the
  outermost equals h; `n_r` and `grading` can be fixed instead.

Deck-driven: `scripts/plate_hole_mesh.py decks/pfd_study/mesh/<case>.i`.
The `[Mesh]` block carries width, height, thickness, h, hole (circle |
ellipse | none), diameter (one value, or 'major minor' full axes for an
ellipse), angle, centre, ring_cells, n_r, grading, the Part I orientation
point file (for the grain-id figure and `elem_grain.csv`) and the figure
path. Outputs in `out_dir`: `mesh.e` (Exodus II; node and side sets `x-`,
`x+`, `y-`, `y+`, `z-`, `z+`, `hole`, node set `outer`; blocks `plate` = 1
and `ring` = 2), `mesh.vtk`, `mesh.json` (hexaform's summary and
validation report), `elem_grain.csv` (element id, grain id).

Check before any run: `validation.passed` true in `mesh.json`, no
inverted element, and the figure shows the ring where the hole is.

## Stage 2 -- grains and orientations

Unchanged from Part I: the 30-grain structure is the phase-field growth
run labelled on its own 0.49 um grid, and the orientation point file
(`orientations_points.csv` of the calibrated case) assigns each element
the grain of the nearest point to its centroid. The hole removes the
elements inside it and nothing else; the grain map is independent of the
mesh, so every mesh level and every hole shape sees the same grains.

## Stage 3 -- PINNACLE deck

`[Mesh] type = file path = mesh.e format = exodus` with the sets above;
`[OrientationField] source point_file` as in Part I; the staggered
executioner over the mechanics sub-solve (`fem_quasi_static` for J2,
`cp_quasi_static` for CP) and `pfd_steady` for the damage; damage
Dirichlet d = 0 on `x-` and `x+` only (the loaded faces); the engineering
stress is the reaction on `x+` over the gross section. The deck template
and its pipeline script are written with the first J2+PFD case.

## Surrogate comparison (Part II)

After the surrogate decks have run on the host, pull `h40_curve/`, `h40_field/`
(with `fields.vtk` for the field maps) and the recollected `h40_data/samples.csv`,
then from the study root:

```
python3 scripts/surrogate_compare.py decks/pfd_study/cp/surrogate/compare_h40.i
```

It writes `results/pfd_study/cp/surrogate/cp_surrogate_h40_summary.csv` and the
held-out curve, held-out field-map and cost-per-query figures under `figures/pfd/`.
