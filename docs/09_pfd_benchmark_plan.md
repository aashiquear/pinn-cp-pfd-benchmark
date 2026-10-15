# Phase-field ductile fracture (PFD) and coupled CP-PFD benchmark: plan

Plan of 2026-10-13, written after the author's go-ahead for Part II of the
benchmark (PFD on J2 plasticity, then CP-PFD) and before any simulation.
Status 2026-10-15: executed; the author's revision round of 2026-10-14
changed three things against this plan (h160 dropped and h20 added to the
mesh levels; the J2 configuration carries W_c = 200 MPa from a threshold
scan, docs/10 §4.10; the shape study of §6 ran on CP+PFD, docs/10 §9), the
results are in `docs/10_pfd_study.md` and the report §12-§20.
Nothing in this document is a result; results go to `docs/10_*.md` once
they exist and to the report (`docs/reports/cp_pfd_pinnacle_benchmark_report.md`)
once agreed.

## 1. Author's guidelines (2026-10-13)

1. Benchmark the PFD on J2 plasticity first, then the coupled CP-PFD.
2. Same plate (100 x 50 x 1 um) and mesh resolution as Part I, with a
   central hole as stress concentrator. The hole diameter follows the
   crack-tip plastic-zone concept: the plastic zone must be unaffected by
   the other boundaries; a 5 um hole on the 100 x 50 surface satisfies it.
   Irwin's model gives the plastic zone. Polycrystal again.
3. Meshes from hexaform, hexahedral only; the mesh is shown to the author
   before any run.
4. Calibrate G_c and l; eta is fixed (quasi-static). The brittle rule of
   thumb l > 2h has no J2- or CP-coupled counterpart, so the mesh
   sensitivity reproduces the curve of (chakraborty 2016 multiscale) Fig. 8:
   stress-strain, peak stress and fracture energy against l (the mesh rule
   they used was h = l/3, from their Fig. 5). Calibration target: the last
   point of the (vaughan 2024 mechanistic) stress-strain curve is the
   fracture initiation point. Fracture initiation = damage d = 0.8 at any
   node or element average; the strain at which this first happens is the
   fracture strain, and G_c and l are calibrated so it equals the last
   experimental strain.
5. After calibration and mesh sensitivity: a stress-concentrator shape
   study (hole, ellipse, ...), since the Irwin plastic zone depends on the
   shape. Research and propose (section 6).
6. One report for Part I and Part II (renamed
   `cp_pfd_pinnacle_benchmark_report.md`), final version only; everything
   else in separate documents. Abstract and introduction stay empty.

Hardware (author, 2026-10-13): use the best available hardware for every
run except the hardware comparison. Section 2 records what "best" is.

## 2. Hardware decision

The torch CPU arm is the production arm for the finite-element solves.
The 2026-10-13 thread-scaling addendum to the hardware benchmark (report
§10, Table 12; `results/cp_study/hardware/hardware_h40.csv`):

| FE arm | threads | Newton time per iteration | memory |
|---|---|---|---|
| numpy | 1-8 | 0.45 s | 175 MB |
| torch CPU | 1 | 0.142 s | 1.18 GB |
| torch CPU | 2 | 0.092 s | 1.19 GB |
| torch CPU | 4 | 0.076 s | 1.20 GB |
| torch CPU | 8 | 0.067 s | 1.20 GB |
| torch GPU | -- | 0.079 s | 1.61 GB |

The gain is the batched kernel of ADR 0027 (3.2x on one thread), not the
threads (a further 2x from 1 to 8 threads with diminishing returns). The
direct linear solver (2026-10-09, `linear_method = direct`) cured the CG
wall on every backend and is unrelated to the torch gain. The host has 16
physical cores and 62 GB: campaigns run **seven cases side by side at two
threads each** (0.092 s per iteration, 8.4 GB); single runs use eight
threads. The GPU is reserved for network training (the neural routes and
the surrogate) and for the hardware comparison.

## 3. Plastic zone at the hole (Irwin) and what it implies

Calibrated Part I plate: 0.2%-offset yield 364 MPa, flow stress 593 MPa at
the last experimental point (17.2% nominal strain). Hole radius a = 2.5 um.
The plate treated as a centre-cracked panel of half-length a with K_I =
sigma sqrt(pi a) (the Creager-Paris equivalence between a blunt notch of
root radius rho and a crack, (creager 1967 blunt): the root stress
2 K_I / sqrt(pi rho) equals the Kirsch 3 sigma for the circle), Irwin's
plastic-zone radius r_p = K_I^2 / (2 pi sigma_y^2) (plane stress) and
K_I^2 / (6 pi sigma_y^2) (plane strain):

| remote stress | r_p, plane stress | r_p, plane strain |
|---|---|---|
| sigma_y / 3 = 121 MPa (first yield at the hole edge, Kirsch factor 3) | 0.14 um | 0.05 um |
| sigma_y / 2 = 182 MPa | 0.31 um | 0.10 um |
| sigma_y = 364 MPa | 1.25 um | 0.42 um |
| 593 MPa (flow stress at the fracture strain) | 3.3 um | 1.1 um |

The elastic Kirsch field with a von Mises check gives the same picture:
at sigma = 182 MPa the yielded zone reaches 0.56 um from the hole edge;
at 291 MPa 2.3 um. The ligament from the hole edge to the plate edge is
22.5 um and to the loaded faces 47.5 um, so the plastic zone of the hole
(at most 3.3 um even at the fracture stress, and the Irwin estimate is
only valid below about sigma_y / 2) is 7 to 14 times smaller than the
nearest boundary distance: the author's condition holds with margin.
Two consequences the plan must respect:

- **The plastic zone is smaller than a grain.** The 30-grain structure has
  a median equivalent grain diameter of 13.7 um; the hole (5 um) sits at a
  junction of four grains (Fig. `figures/pfd/mesh_h40_hole5.png`). Damage
  initiation at the hole is therefore governed by the orientation of the
  grains touching it, which is the microstructural effect the CP-PFD stage
  is meant to show, and the grain-instance check of Part I (§7) must be
  repeated for the CP-PFD case.
- **Beyond sigma_y / 3 the whole net section yields** (the plate reaches
  17% nominal strain): the hole is a strain concentrator in a fully plastic
  plate, not a small-scale-yielding crack. The Irwin numbers size the
  process zone for the mesh (section 4), they do not predict the fracture
  strain.
- **The hole-edge element must resolve the process zone and l.** The mesh
  of section 4 has 0.65 um elements at the hole edge at h40 and 0.49 um at
  h80; the Irwin zone at yield (1.25 um) is two to three elements, and l
  (section 5) is resolved by the ring whatever the background h.

## 4. Mesh (hexaform `plate`, feat/plate-hole 12a29d5)

hexaform could only mesh cubes with closed inclusions; the plate generator
was added on 2026-10-13 (`hexaform plate`, `hexaform.plate.plate_mesh`,
tests, user guide, example 12): a uniform in-plane background of cell size
h, one HEX8 through the thickness, and a block-structured O-grid ring whose
nodes lie exactly on the hole curve (circle or rotated ellipse), graded so
the innermost radial spacing equals the hole-edge arc length (near-square
elements where damage starts) and the outermost equals h; the ring meets
the background on a grid-aligned rectangle with no transition elements.
Sets: `x-`, `x+`, `y-`, `y+`, `z-`, `z+`, `hole`, `outer` (the four
in-plane edges without the hole, for the damage Dirichlet set).

Study script `scripts/plate_hole_mesh.py` (deck-driven, `decks/pfd_study/mesh/`)
writes `mesh.e` (Exodus, for PINNACLE `[Mesh] type = file`), `mesh.vtk`,
`mesh.json` (summary and hexaform validation), `elem_grain.csv` (grain id
per element from the Part I point file, nearest point at the element
centroid, exactly PINNACLE's `point_file` rule) and the figure.

| mesh | h (um) | elements | nodes | ring layers x segments | hole-edge element (um) | min scaled Jacobian | max edge ratio |
|---|---|---|---|---|---|---|---|
| h40_hole5 | 2.5 | 860 | 1,864 | 4 x 24 | 0.65 | 0.71 | 3.3 |
| h80_hole5 | 1.25 | 3,232 | 6,736 | 3 x 32 | 0.49 | 0.71 | 2.9 |
| h160_hole5 | 0.625 | 12,800 | 26,128 | 3 x 48 | 0.33 | 0.71 | -- |

Both pass hexaform's validation (conforming, no inverted element). The
0.71 Jacobian is the 45-degree skew of the ring elements at the rectangle
corners. Figures: `figures/pfd/mesh_h40_hole5.png`, `mesh_h80_hole5.png`
(plate coloured by grain id, zoom on the ring). The same 30-grain
structure and orientations as Part I (point file of the calibrated case)
apply unchanged: the grain structure is independent of the mesh.

Open for the author: whether the hole at the grain junction is the wanted
configuration (it is where the Part I structure puts it; moving the hole
or drawing another grain instance are both one deck line).

## 5. Study design

### 5.1 Stage order

1. **J2 fit** (0-D, no FE): fit `j2_isotropic_hardening` (E, nu, sigma_y0,
   K, n; or a piecewise curve) to the digitized Vaughan curve; verify on the
   plain h40 plate against the experiment and the Part I CP curve.
2. **J2+PFD calibration** at the Part I mesh (h40 background, ring at the
   hole): for each l of the mesh-rule set, bisection on G_c so that d = 0.8
   is first reached at 17.2% nominal strain (within one load step, 7.8e-4).
   One scalar target cannot fix two parameters: the deliverable is the
   calibrated G_c(l) curve and the crack path for each pair; l is then
   chosen by the mesh study.
3. **Mesh and length-scale sensitivity** (the Fig. 8 analogue): stress-
   strain curve, peak stress, fracture strain and fracture energy against
   l at fixed l/h (h40, h80, h160 background with the ring scaled) and
   against h at fixed l. Deliverable: the l/h at which the J2-PFD and
   CP-PFD answers stop moving (the rule the author asked for), and the
   production mesh.
4. **Parameter sensitivity** of G_c and l (and the fixed eta, once, to
   show it does not matter at the quasi-static rate) at the chosen mesh,
   one at a time, FE only.
5. **CP-PFD**: steps 2 to 4 with the Part I crystal-plasticity set, plus
   the grain-instance check.
6. **Neural routes** on the calibrated CP-PFD case: PINN+FEM (linear
   extrapolation predictor on the mechanics sub-solve of the staggered
   loop), PINN as solver (grain-aware verified), parametric surrogate over
   (G_c, l, and the CP set) with the fracture strain and the crack path as
   targets; accuracy, time and hardware against the FE.
7. **Hardware and time benchmark** on the calibrated case and mesh.
8. **Shape study** (section 6).

### 5.2 Definitions fixed now

- Fracture initiation: first load step at which max over nodes of d >= 0.8
  (and, reported alongside, the element-average criterion); the crack path
  is the set of elements with average d >= 0.8 at the end of the run.
- Fracture strain: nominal strain at that step; target 17.2% (the last
  point of the digitized curve, `al2219_fig3_calibration_stress_strain.csv`;
  the last marker is 17.0%).
- Stress reported: the reaction force on the loaded face divided by the
  gross section (50 x 1 um^2), the engineering stress of the perforated
  plate; the volume-mean stress of Part I is not the nominal stress of a
  plate with a hole.
- Damage Dirichlet set: `outer` is NOT pinned (a crack from the hole must
  reach the free edges); the loaded faces `x-`, `x+` carry d = 0 (grips).
- eta fixed; AT2 unless the mesh study says otherwise; degradation
  quadratic; the plastic-work driving force with the threshold decided in
  the development (section 7).

## 6. Stress-concentrator shape study (proposal)

Inglis/Kirsch: an elliptical hole with semi-axis a transverse to the load
and b along it has root radius rho = b^2 / a and stress-concentration
factor K_t = 1 + 2 a / b = 1 + 2 sqrt(a / rho); the circle is K_t = 3. In
the Creager-Paris form the Irwin plastic-zone radius r_p = (a / 2)
(sigma / sigma_y)^2 depends on the transverse half-extent a only, while the
elastic-plastic zone shape and the strain at the root depend on rho. The
study therefore separates the two:

| series | shapes (full axes 2a x 2b, um) | what is held | K_t | r_p at sigma_y (um) |
|---|---|---|---|---|
| A, transverse extent | 5 x 5 (circle), 5 x 2.5, 5 x 1.25 | a = 2.5 | 3, 5, 9 | 1.25 for all |
| B, area | 5 x 5, 7.07 x 3.54, 10 x 2.5 | hole area 19.6 um^2 | 3, 5, 9 | 1.25, 1.77, 2.5 |
| C, orientation | 5 x 2.5 at 0, 45, 90 degrees | shape | 5 (0), 2 (90) | 1.25 .. 0.63 |

Nine meshes (three shared), all from the same generator (`hole = ellipse`,
`diameter = 'major minor'`, `angle`), same background h and ring rule.
Predictions to test: with r_p fixed (series A) the fracture strain should
fall with K_t and the crack path stay transverse; with K_t fixed (series B)
the fracture strain should fall with r_p; rotating the ellipse (C) should
move initiation from the root to the grain boundary nearest the root and
tilt the path. Outputs per case: fracture strain, crack path, peak
engineering stress, fracture energy, plastic-zone extent at first yield
(von Mises >= sigma_y) against the Irwin value. Cost: nine CP-PFD runs at
the chosen mesh, about 9 x the calibrated run. Resolution caveat: the
sharpest root (rho = 0.31 um for 5 x 1.25) needs a hole-edge element of
0.3 um or less; at h40 the ring gives 0.6 um, so series A's sharpest case
runs at the h80 ring (0.49 um with the default rule, finer with a wider
transition rectangle) or is dropped to a/b = 3. The author decides after
the mesh study.

## 7. PINNACLE: what exists, what is built before the first run

Survey of 2026-10-13 (feat/pinn_benchmark fdc911d). Exists: 3-D hex8
staggered quasi-static J2+PFD and CP+PFD (`staggered` over
`fem_quasi_static` / `cp_quasi_static` and `pfd_steady`; ADR 0024), AT1/AT2,
quadratic / cubic / Lorentz / parsed g(d), Miehe spectral driving force with
the history field H, the ductile term H = max(H, psi+ + zeta W_p)
(`pfd_driving_force`, `plastic_work_prefactor`), stress degradation
(spectral split for J2, volumetric Je split for CP), Exodus import, point-
file orientations, restart checkpoints, VTK series with damage and H,
PINN+FEM hybrid with damage on the PINN, PINN-as-solver for the CP
mechanics (ADR 0028/0029), data-source surrogate (ADR 0030).

To build (PINNACLE rules: general, deck-driven, tested, documented, ADR):

1. **Sparse damage system**: `phasefield/damage_field.py` assembles a dense
   n x n matrix and `pfd_steady` solves it densely; CSR assembly and the
   existing sparse solvers (`linear_method`), required for any plate.
2. **Damage boundary set**: `pfd_steady` pins every node when the named set
   is empty and its default `all_boundary` includes the hole surface; an
   explicit `none` and a validated named set.
3. **Reaction force per step** on a named node set for `fem_quasi_static`
   and `cp_quasi_static` (the engineering stress of a perforated plate) and
   de-duplicated staggered history rows (last alternation per step: load,
   reaction, mean stress, max d, fracture-initiation flag).
4. **Viscous regularisation eta** in the damage equation (the author fixes
   it, it must exist to be fixed), **plastic-work threshold W_c** (Borden-
   style, so purely elastic steps do not drive damage) and the choice of
   the driving-force split documented against the stress split (the
   survey found the driving force uses an isotropic Voigt-average Miehe
   split while the J2 stress uses a spectral split and the CP stress a
   volumetric split).
5. **Exodus reader**: element sets keyed by block name (several hex blocks
   collide today), element variables read, element map kept, so hexaform's
   `ring` block and a `grain_id` element variable arrive intact;
   `elem_grain_file` from the mesh script already works as the fallback.
6. **PINN+FEM and PINN-as-solver inside the staggered loop**: the mechanics
   sub-solve of `staggered` must reach `load_step_services` (today only the
   stand-alone `cp_quasi_static` does), so the neural routes are per-step
   comparable with the FE; the hybrid coupled-physics executioner has no
   outer load loop (ADR 0018) and is not the route of this benchmark.
7. **Surrogate on PFD**: `PFDKernel` flagged `history_dependent`, staggered
   history rows with curve columns, the damage field mapped as a field
   target.
8. **Examples, goldens, docs**: a perforated-plate J2+PFD and CP+PFD example
   on an external hex mesh with the plastic-work driving force (none today),
   ADR for the sparse damage solve and the reaction/threshold additions,
   CHANGELOG, syntax regeneration.

Order: 1-3 and 5 first (needed by the first J2+PFD run), 4 with the
calibration, 6-7 before the neural stages.

## 8. Files of this stage

- Meshes and decks: `decks/pfd_study/mesh/{h40,h80}_hole5.i`,
  `results/pfd_study/mesh/*/mesh.{e,vtk,json}`, `elem_grain.csv`,
  `figures/pfd/mesh_*_hole5.png`.
- Scripts: `scripts/plate_hole_mesh.py` (this study),
  hexaform `feat/plate-hole` (12a29d5): `src/hexaform/plate.py`,
  `tests/test_plate.py`, `examples/12_plate_with_hole.py`, user guide.
- Workflow: `docs/workflow_pfd_plate.md`.
