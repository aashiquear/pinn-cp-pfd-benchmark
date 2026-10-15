# PINNACLE CP–PFD benchmark study

Comparative study of FEM, PINN+FEM hybrid, and neural-constitutive solvers for
coupled crystal plasticity and phase-field damage, validated against
(Rezwan 2025 intragranular) — Al2219 single crystal with an intracrystalline
void — and (Rezwan 2026 nickel) — single-crystal Ni tensile specimens.

**Study files live here. Framework changes live in
`/mnt/extended/projects/pinnacle` on branch `feat/pinn_benchmark`. Meshing
tooling lives in `/mnt/extended/projects/hexaform`.**

Read `docs/PROJECT_GUIDELINES.md` first — it binds every conversation in this
project.

## Layout

| Path | Contents |
|---|---|
| `bib/references.bib` | The single central bibliography, keys `author-year-significantword` |
| `data/extracted/`, `data/digitized/`, `data/reference/` | Provenance-tagged dataset, per-figure digitized data, the MOOSE reference curve |
| `moose/` | The MOOSE reference decks; their Exodus outputs stay on disk, ignored by git |
| `decks/graingrowth/`, `decks/polycrystal_plate/`, `decks/cp_study/` | Part I: grain growth, the CP plate pipeline, the CP calibration / mesh / sensitivity / neural / surrogate studies |
| `decks/pfd_study/` | Part II: meshes (`mesh/`), the J2 fit, J2+PFD (`j2/`) and CP+PFD (`cp/`) studies, the surrogate (`cp/surrogate/`), the report field decks (`report/`), templates |
| `decks/benchmark/`, `decks/report/` | Hardware benchmarks; report figure, comparison and shipping decks |
| `scripts/` | Deck-driven drivers: `polycrystal_plate.py`, `cp_parameter_study.py` (Part I); `plate_hole_mesh.py`, `pfd_plate.py`, `pfd_study.py`, `pfd_calibrate.py`, `pfd_*_figures.py`, `pfd_fields.py` (Part II); `hardware_benchmark.py`, `surrogate_design.py`, `surrogate_compare.py`, `ship_results.py` |
| `results/` | Every run's deck, history, curve and summary; the field series, meshes and networks the report presents (`results/SHIPPED.md`, 239 MB); other field files are not committed |
| `figures/report/`, `figures/pfd/` | Part I and Part II figures (the animations `figures/pfd/*_evolution.gif`) |
| `archive/` | The first campaign (Sept 2026); see `archive/README.md` |

## Reports and notebooks

| Document | Contents |
|---|---|
| `docs/reports/cp_pfd_pinnacle_benchmark_report.md` | The technical report: Part I (CP, §2-§11) and Part II (PFD on J2 and CP, §12-§20); abstract and introduction pending |
| `docs/PROJECT_GUIDELINES.md` | Project-wide rules |
| `docs/run_ledger.md` | One row per run: deck, scripts, revision, host, outcome |
| `docs/07_cp_parameter_study.md` | Part I lab notebook (CP calibration, mesh, sensitivity, neural routes, hardware, surrogate) |
| `docs/09_pfd_benchmark_plan.md`, `docs/10_pfd_study.md` | Part II plan and lab notebook (PINNACLE readiness, J2 fit, W_c scan, calibrations, sweeps, fields, neural routes, hardware, surrogate, shape study) |
| `docs/08_pinn_route_assessment.md` | What PINNACLE has for CP on a PINN / PINN+FEM and the development plan |
| `docs/00_capability_audit.md`, `docs/02_backend_envelope.md`, `docs/03_material_point_verification.md`, `docs/05_orientation_notation.md`, `docs/06_polycrystal_plate_results.md` | Phase 0 and the plate's first runs |
| `docs/workflow_polycrystal_plate.md`, `docs/workflow_pfd_plate.md`, `docs/reproduce_polycrystal_plate.md` | How the pipelines work and the exact commands on the ROCm host |
| `docs/reconciliation.md`, `docs/digitization_uncertainty.md` | Source extraction and its error |

## Status (2026-10-15)

Part I (crystal plasticity on the 100 x 50 x 1 um, 30-grain Al2219 plate):
calibrated CP set at mesh h40, mesh and parameter sensitivity, grain
instances, the PINN+FEM and verified PINN-as-solver routes, the hardware
benchmark and the parametric surrogate (report §2-§11).

Part II (phase-field ductile fracture on the same plate with a 5 um hole):
J2+PFD carried at W_c = 200 MPa, Gc = 385.5 MPa um, l0 = 5 um (W_c scan);
CP+PFD at Gc = 1106 MPa um, l0 = 5 um, W_c = 0; both with length-scale,
four-value parameter and h20/h40/h80 mesh sensitivities, the J2 against
CP field comparison and animation, the stress-concentrator shape study,
the neural routes inside the staggered loop, the coupled hardware
benchmark and the surrogates (report §12-§20). PINNACLE changes for Part
II are ADR 0031 on `feat/pinn_benchmark`.

Carried limits: 46 of 93 Al2219 orientations recoverable, loading-axis
only; Ni failure stress not recovered; gaps C.15 and C.16 open; the J2
initiation strain with the threshold is not mesh-converged at h40
(report §13.5), kept at h40 for the side-by-side comparison.
