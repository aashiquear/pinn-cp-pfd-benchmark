# Workflow: phase-field polycrystal plate to a CPFE run

How a pseudo-3D polycrystal plate is built and run in this study. The
companion documents are `reproduce_polycrystal_plate.md` (the exact commands
and inputs) and `06_polycrystal_plate_results.md` (what came out). This one
explains the *why* of each stage and what each stage owes the next.

## Reference geometry

The plate is the reference MOOSE deck's (`moose/CPFE.i`): 100 x 50 x 1 um,
50 x 25 x 1 HEX8 elements (2 um in-plane), pulled along +y from the top face
with u_x = 0 on the left, u_y = 0 on the bottom and u_z = 0 on the back face.
Loading is 0.15 um/s for 10 s: 3.0% nominal strain at 3e-3 /s in 40 steps of
0.25 s. The MOOSE output (`moose/CPFE_out.e`, single crystal, Euler (0,0,0))
is summarised in `data/reference/moose_cpfe_plain_plate_stress_strain.csv`
by `scripts/tools/exodus_summary.py`.

## Stage 1 -- grain structure by phase-field grain growth

A Voronoi tessellation has straight boundaries and arbitrary junction angles.
The structure used here is the Voronoi seeding **evolved** by the Moelans
multi-order-parameter Allen-Cahn model (`gb_evolution` material, `ac_bulk` +
`ac_interface` kernels, `pf_fem_transient` executioner) on its own fine 2-D
grid of the plate's in-plane extent. One element through the thickness makes
in-plane 2-D growth the exact analogue of the pseudo-3D plate.

Settings that matter, with the measurement behind each:

* **Time scheme: `bdf1`, not `convex_splitting`.** The convex-splitting
  stabiliser for the 8-order-parameter Moelans energy is S = 147.5 mu
  (Gershgorin bound over the band [-1.5, 1.5]), so even at dt = 0.25/(L mu)
  the implicit lag factor dt L S is about 37 and every step advances a fixed
  amount of the explicit force whatever dt or the mobility is. That is the
  "inert mobility" the earlier campaign recorded. Measured on the 16-grain
  deck over one interface time: convex splitting dropped the free energy by
  29.5 (per-step drop 0.92 at dt = 0.113 and 1.04 at dt = 0.227), bdf1 by
  107 and dt-converged to 0.2% (E vs F in `results/graingrowth_verify/`),
  with M_GB = 2 at half the steps reproducing M_GB = 1 to eight digits.
* **Step: dt = 0.5/(L mu)** where 1/(L mu) = w_GB^2/(8 M_GB sigma) is the
  bulk relaxation time. Newton takes 4-5 iterations per step and no
  cutbacks.
* **Grid: h = w_GB/3, w_GB = grain width/8.** The material's own guidance.
* **Horizon** is chosen from the pilot (`decks/graingrowth/pilot/`), which
  ran the same seeding for 40, 80 and 160 steps (2.5, 5 and 10 interface
  times). The results document records the morphology at each and which one
  the mechanics used.
* **Boundary:** natural (no-flux), so boundaries meet the plate edge at 90
  degrees, as in MOOSE. No padding or sectioning.

Output: `fields.vtk` with the final nodal `eta_k`, and `solver_history.csv`
with the free energy per step.

## Stage 2 -- labelling, orientations, deck (`scripts/polycrystal_plate.py`)

Driven by one pipeline deck (`decks/polycrystal_plate/*.i`, PINNACLE `.i`
syntax):

1. `label`: argmax colour per grid node, connected components within each
   colour (graph colouring puts non-adjacent grains on one order parameter),
   specks under `min_cells` dropped, diffuse-boundary nodes filled from the
   nearest labelled node. Morphology metrics from
   `scripts/grain_morphology.py` are written beside the labels.
2. `orient`: one Bunge triple per grain, uniform on SO(3) from the recorded
   seed (Phi drawn as arccos of a uniform variate). Written as a PINNACLE
   **point file** -- one row per growth-grid node with its grain's
   orientation, z at mid-thickness -- and as `grains.csv`.
3. `deck`: `decks/polycrystal_plate/cpfe_plate_template.i` filled from the
   pipeline deck's `[Mechanics]` block; the point file is referenced from the
   deck through `[OrientationField] source = point_file` (PINNACLE ADR 0026),
   so the run is reproducible from `cpfe.i` plus the CSV, with nothing
   injected by a script.

## Stage 3 -- the CPFE run and its curve

Optional `[Mechanics] snapshots_every = N` (default 0) makes the generated
deck record the fields every N accepted load steps (PINNACLE `cp_quasi_static`
`snapshots_every`, commit 6517a31): the run then also writes
`fields_step{NNNN}.vtk` + `fields.pvd`, which `scripts/plastic_strain_animation.py`
turns into the equivalent-plastic-strain animation of the report. Each frame
costs one copy of the per-point fields in memory and one VTK file on disk
(about 10 MB per frame at h160), so the study uses N = 4.

### Further grain instances

`scripts/make_growth_instance.py TEMPLATE.i --seed S --out NEW.i --output-dir D`
writes a growth deck identical to the pilot's except for a new uniform
Voronoi seeding drawn with `default_rng(S)`; the coordinates are written into
the new deck, so the instance is reproducible from the deck alone. A pipeline
deck then points `[Growth] run_dir` at that run and `[Orientations] seed` at a
new draw (`decks/cp_study/grain_s{A,B}.i`).

`pinnacle run cpfe.i`. The `cp_quasi_static` executioner records, per
accepted load step, the ramp time and the quadrature-point mean stress and
strain tensors (PINNACLE commit 14dd45c), so `solver_history.csv` is the
stress-strain curve. `polycrystal_plate.py ... curve` turns it into
`results.csv` with the nominal strain from the prescribed displacement.

## What is deliberately not done here

* No damage: the first full run is CP only, as the plan requires. The
  template has no `d` variable; the `pfd_*` blocks of the earlier plate
  template are not carried.
* No parameter fitting: the mechanics values are the reference deck's. CP
  parameters are model-, framework- and mesh-sensitive, so the published set
  is a starting point for calibration, not a transferable truth.
