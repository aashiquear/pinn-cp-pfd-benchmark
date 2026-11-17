# PINNACLE CP-PFD benchmark study

Comparative study of the finite-element, PINN+FEM, PINN-as-solver and
surrogate routes of PINNACLE for crystal plasticity and phase-field ductile
fracture, on two problems: the Al2219 polycrystal plate of
(Rezwan 2025 intragranular) with the tensile curve of
(Vaughan 2024 mechanistic), and the single-crystal Ni dog-bone specimens of
(Rezwan 2026 nickel).

**Study files live here. Framework changes live in
`/mnt/extended/projects/pinnacle` on branch `feat/pinn_benchmark`. Meshing
tooling lives in `/mnt/extended/projects/hexaform`.**

Read `PROJECT_GUIDELINES.md` first: it binds every conversation in this
project.

## What is shipped

The repository ships what is needed to reproduce the results, figures and
analysis: input decks, scripts, digitized data, figures and documents. It
does not ship PINNACLE output (`results/`, `archive/results/`) or any
published paper (`data/literature/`, `*.pdf`). The reproduction documents
give the commands that regenerate the output.

## Reports

| Document | Contents |
|---|---|
| `docs/reports/cp_pfd_pinnacle_benchmark_report.md` | Al2219 plate. Part I: crystal plasticity, calibration, mesh and parameter sensitivity, grain instances, neural routes, hardware benchmark, surrogate (§2-§11). Part II: phase-field ductile fracture on J2 and on crystal plasticity, shape study, neural routes, hardware benchmark, surrogate (§12-§20) |
| `docs/reports/ni_single_crystal_pinnacle_benchmark_report.md` | Ni single crystal: one parameter set, 15 specimens on four routes against the experiment |

## Layout

| Path | Contents |
|---|---|
| `bib/references.bib` | The single bibliography, keys `author-year-significantword` |
| `data/extracted/`, `data/digitized/`, `data/reference/` | Provenance-tagged dataset, per-figure digitized data, the MOOSE reference curve |
| `moose/` | The MOOSE reference decks |
| `decks/graingrowth/`, `decks/polycrystal_plate/` | Grain growth and the plate templates |
| `decks/cp_study/bunge/` | Al2219 Part I: calibration, mesh, sensitivity, grain instances, neural routes (`hybrid/`), surrogate (`surrogate/`) |
| `decks/pfd_study/` | Al2219 Part II: meshes (`mesh/`), templates, the J2 fit and J2+PFD studies (`j2/`), the CP+PFD studies and surrogate (`cp/bunge/`), the field decks (`report/bunge/`) |
| `decks/benchmark/bunge/`, `decks/report/bunge/` | Hardware benchmarks; report figure and comparison decks |
| `decks/ni/` | Ni specimens: meshes, templates, crystal plasticity, fracture and the final four-route decks (`*/bunge/`) |
| `scripts/` | Deck-driven drivers and figure scripts; `scripts/tools/` holds the host chains |
| `figures/report/`, `figures/pfd/`, `figures/ni/` | Figures of the two reports |
| `archive/` | Superseded campaigns, decks, notebooks and report versions; see `archive/README.md` |

## Documents

| Document | Contents |
|---|---|
| `PROJECT_GUIDELINES.md` | Project-wide rules |
| `docs/run_ledger.md` | One row per run: deck, scripts, revision, host, cores, outcome |
| `docs/reproduce_al2219_plate.md`, `docs/reproduce_ni_specimen.md` | The commands that rebuild each report |
| `docs/workflow_polycrystal_plate.md`, `docs/workflow_pfd_plate.md` | How the plate pipelines work and why |
| `docs/02_backend_envelope.md`, `docs/03_material_point_verification.md` | Backend timing envelope and the material-point verification of the crystal-plasticity law |
| `docs/05_orientation_notation.md` | Orientation notation of the source paper |
| `docs/13_history_aware_pinn_review.md` | Review of history-aware network routes and what PINNACLE implements of it |
| `docs/reconciliation.md`, `docs/digitization_uncertainty.md` | Source extraction and its error |
| `docs/literature/` | Literature research kept as evidence (no published paper is shipped) |

## Status (2026-11-17)

Both reports are complete on PINNACLE 9548123: Euler angles in Bunge's
convention (ADR 0036), total-Lagrangian balance (ADR 0035), implicit
constitutive update (ADR 0033). Earlier versions read two Euler angles in
exchanged order; they are kept under `archive/docs/` with their decks and
figures.

Al2219 plate (100 x 50 x 1 um, 30 grains, mesh h40): calibrated set
rho0 53.93 um^-2, k1 0.1768, y_c 0.00271 um, tau_c0 5.39 MPa, within 1% of
the experiment from 0.6% to 17% strain. CP+PFD on the plate with a 5 um
hole: G_c 1,123 MPa um at l0 5 um, fracture at 17.1%. J2+PFD as in the first
study (W_c 200 MPa, G_c 385.5 MPa um).

Carried limits: h40 is 5.3% above the finest mesh at 17% strain and is kept
by the author's decision (report §5.3); the perforated CP+PFD plate lies 2 to
7% below the plain-specimen markers (report §14.1); the J2 initiation strain
with the threshold is not mesh-converged at h40 (report §13.5).
