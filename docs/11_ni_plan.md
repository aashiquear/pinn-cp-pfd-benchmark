# Ni single-crystal benchmark (Part III): plan

Plan of 2026-10-15, written after the author's instructions for the third
study and before any simulation. Nothing in this document is a result;
results go to `docs/12_ni_study.md` once they exist and to the Part III
report (`docs/reports/ni_single_crystal_pinnacle_benchmark_report.md`)
once agreed.

## 1. Author's guidelines (2026-10-15)

1. Reference: (rezwan 2026 predicting), `data/literature/Rezwan_2026_Ni_single_crystal.pdf`;
   the earlier extraction had mistakes, so the data are re-extracted from
   the paper (§2).
2. A new report in the format of the Al2219 one, without abstract and
   introduction for now; everything documented in `docs/` first, the report
   written after the results are confirmed.
3. Same campaign as for Al2219 but on single crystals (no polycrystal data):
   the ground truth is one single-crystal experimental curve (Fig. 3,
   specimen L) from a specimen with a given orientation and thickness
   profile; each of the five orientations has three thickness variants, so
   15 specimens in all.
4. Specimens are dog-bone coupons. If hexaform can build them, every
   specimen gets a dog-bone mesh whose ends have the same thickness across
   the specimens and whose 8 mm gage carries the measured thickness
   variation.
5. CP only first: calibrate the same parameter set as for Al2219
   (rho_0, k_m, y_c, tau_c0); elastic constants, Burgers vector, shear
   modulus and the other constants stay at the paper's Ni values. Then mesh
   and parameter sensitivity (no grain-instance check), then the neural
   routes and the surrogate.
6. Then CP+PFD (no J2: Part II showed its inaccuracy for this class of
   model): calibrate G_c and l_0, report h against l_0; compare fracture
   initiation (strain) and fracture shape with the paper's experiments
   (Fig. 4(b) experiment axis, Fig. 5), not with the paper's MOOSE results;
   our FEM is the reference for the PINN, PINN+FEM and surrogate routes,
   with hardware and time cost.
7. No stress concentrator (the thickness variation is the only
   inhomogeneity).

## 2. Data extracted from the paper (2026-10-15)

| item | file | status |
|---|---|---|
| Table 1: five orientations (Bunge angles, max Schmid factor) | `data/extracted/ni_table1_orientations.csv` | transcribed from the text layer, exact |
| Table 2: thickness at five stations of the 8 mm gage, 15 specimens | `data/extracted/ni_table2_thickness.csv` | transcribed, exact |
| Table 3: CP parameters and elastic constants (mm, MPa, s); G_c 125 kJ/m^2 from the text; failure criterion c > 0.85 across the section | `data/extracted/ni_table3_cp_parameters.csv` | transcribed, exact |
| Fig. 2: planform (25.4 mm long, 8 mm gage, 4.8 mm grips, R 3.175 fillets; nominal thickness 0.35 mm); gage width 1.85 +- 0.05 mm measured on the plan view with the grip width as scale (not dimensioned in the paper) | `data/extracted/ni_fig2_planform.csv`, stations `data/digitized/ni_fig2_thickness_stations.csv` | digitized |
| Fig. 3: specimen L calibration curve, experiment and the paper's simulation | `data/digitized/ni_fig3_specimenL_calibration.csv` (earlier digitization, checked against the figure: ranges and the elastic branch agree) | digitized |
| Fig. 4(b): experimental and simulated failure strain of the 15 specimens | `data/digitized/ni_fig4b_failure_strain_parity.csv` (earlier digitization; the 15 pairs re-checked point by point against the figure; L.4 and 10.4 occluded, flagged) | digitized |
| Fig. 5: fracture location along the gage and mode (inclined 45-degree crack or necking) of the 15 experiments | `data/digitized/ni_fig5_fracture_descriptors.csv` | visual estimate (+-0.08 of the gage), no scale in the figure |
| Fig. 4(a): the paper's 15 simulated curves | not digitized (MOOSE results, not compared by decision 6) | -- |

What the paper does not give: the gage width (measured, above), the
thickness outside the five stations (spline inside the gage as the paper
did; the author's common 0.35 mm in the grips), the identity of the Fig. 3
specimen among L.1/L.3/L.4 (its experiment reaches 0.513 strain while
Fig. 4(b) gives 0.461-0.481 for the three L specimens; the Fig. 3 curve is
taken as the hardening reference and L.1 as the calibration geometry),
the PFD length scale l_0 and the mesh size (9,880 HEX8 in the gage).

## 3. Geometry and mesh (hexaform `dogbone`, 2026-10-15)

hexaform could not build a dog-bone, so a generator was added
(`hexaform.dogbone`, CLI `hexaform dogbone`): a single mapped hexahedral
block whose width follows the planform (grips 4.8 mm, circular fillets
R 3.175 tangent to the 1.85 x 8 mm gage, total length 25.4 mm) and whose
thickness follows the measured stations by a natural cubic spline in the
gage, blends linearly to the common grip thickness across each fillet and
is constant in the grips. Node sets `x-`, `x+`, `grips`, `gage-`, `gage+`
(the first and last station planes, for the gauge strain), `station_1..5`,
`edge-/edge+`, `z-/z+`; blocks gage / shoulder / grip. The study script
`scripts/ni_specimen_mesh.py decks/ni/mesh/<level>.i` builds all 15 from
Table 2 and draws the planform, the thickness profile and the 15 profiles
(`figures/ni/mesh_h020.png`). Levels: gage cell 0.4, 0.2 (reference) and
0.1 mm, one element through the thickness (pseudo-3-D as in Parts I-II;
two layers are a mesh-study case), grips at twice the gage cell.

| level | gage cell (mm) | elements | nodes |
|---|---|---|---|
| h040 | 0.4 | about 200 | -- |
| h020 | 0.2 | 792 | 1,780 |
| h010 | 0.1 | about 3,200 | -- |

The meshes are shown to the author before any run (rule 3 of Part II).

## 4. Units, loading and measures

- Unit system of this study: mm, MPa, s (the paper's), so Table 3 enters
  the decks unchanged (b = 2.5e-7 mm, rho_0 = 2.35e6 mm^-2, y_c = 0.79e-6
  mm). Parts I-II used um; the CP law is scale-free under a consistent
  system.
- Loading: the grip faces, `x-` fixed in x and `x+` pulled in x by a ramp,
  rigid-body motion removed by one node (y, z); everything else free; the
  nominal grip strain is the ramp over the full length and the
  displacement-controlled rate is set so the gauge strain rate matches the
  experiment's quasi-static order (1e-3 /s on the gage).
- Engineering strain: the gauge strain = (mean u_x on `gage+` minus mean
  u_x on `gage-`) / 8 mm, the FE counterpart of the DIC gage average
  (PINNACLE `gauge_node_sets`, ADR 0032). Engineering stress: the reaction
  on `x+` divided by the mean initial gage section (1.85 mm x the mean of
  the five station thicknesses), as the paper defines it.
- Fracture initiation (CP+PFD): the paper's c > 0.85 across the cross-section
  (every node of one station column above 0.85) is the primary criterion,
  the Part II d >= 0.8 at any node the secondary; the fracture strain is the
  gauge strain of that step; the fracture location is the station column
  that crosses first, the mode the inclination of the d >= 0.8 band.

## 5. Study design (stage order)

1. Meshes at three levels (§3), shown to the author.
2. CP only on specimen L.1, paper's parameter set: a check run, then the
   Gauss-Newton calibration of (rho_0, k_m, y_c, tau_c0) against the Fig. 3
   experiment markers (strains 0.05 to 0.45) with `scripts/cp_gauss_newton.py`
   on the new pipeline (`scripts/ni_specimen.py`). The calibrated set
   is the Ni set of this study.
3. Mesh sensitivity on L.1 (h040, h020, h010, h020 with two layers) and the
   parameter sensitivity (four values per parameter, as the author asked in
   Part II).
4. CP on all 15 specimens with the calibrated set (curves; the
   orientation ordering of the flow stress against Fig. 4(a)'s experiment-
   free ordering is a qualitative check only).
5. Neural routes (linear predictor, verified grain-aware solver on the single
   crystal = one grain) and the surrogate over (rho_0, k_m, y_c, tau_c0, and
   the orientation as three Euler angles) on L.1; hardware benchmark.
6. CP+PFD: G_c(l_0) calibration on L.1 to its Fig. 4(b) experimental
   failure strain (0.463), l_0 sweep, h against l_0 (three mesh levels),
   the 15 specimens; parity of fracture strain against Fig. 4(b)'s
   experiment axis and the location/mode table against Fig. 5.
7. Neural routes and the surrogate on the coupled problem; hardware.
8. Report Part III after the author confirms the results.

Costs (implicit backward-Euler update, grip increment 1e-3): 10 s per load
step on the host's torch arm (2 threads) for the h040 bar (220 elements,
1,760 points; 40 s locally on numpy), so the 360-step path of the
calibration deck (`decks/ni/cp/base_L1_h040.i`, grip strain 36%, gauge
about 54%) is about 1 h per run; a Gauss-Newton iteration (base + 4
perturbations, 5 side by side) is about 1 h; h020 (4x the points) about 4 h
per run; 15 specimens at 7 side by side 2-8 h per stage, and the mesh study
decides whether h040 can be the production level. The first step from the
virgin state takes 24 Newton iterations (the clamped grip corners), the rest
5-7.

## 6. PINNACLE needs (ADR 0032, ADR 0033)

0. A consistent CP state update (ADR 0033, found by the pipeline check of
   docs/12 §1): the explicit scheme of Parts I-II is first order in the load
   increment and would need 30,000 steps on the Ni path; Part III runs
   `state_update = implicit` (backward-Euler slip update) at increments of
   about 1e-3 of gauge strain, about 500 load steps to 50%.
1. `gauge_node_sets` on `cp_quasi_static` and `fem_quasi_static`: the
   gauge strain per step from two node sets (mean displacement difference
   over mean initial distance).
2. The single-crystal orientation from the crystal block's `euler_angles`
   (exists) without an orientation field.
3. The through-thickness damage criterion of the paper (c > 0.85 across a
   section): a study-script postprocessor on the field series first; a
   framework criterion only if the staggered stop needs it.
