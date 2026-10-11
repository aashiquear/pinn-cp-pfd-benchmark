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
| `data/extracted/` | `benchmark_dataset.csv` (provenance-tagged rows), `gap_register.csv` (C.1–C.16) |
| `data/digitized/` | Per-figure extracted data |
| `data/reference/` | `moose_cpfe_plain_plate_stress_strain.csv`, the MOOSE test deck's curve |
| `moose/` | The MOOSE reference decks (CP-only plate, coupled CP-PFD test); their Exodus outputs stay on disk, ignored by git |
| `decks/graingrowth/` | Phase-field grain-growth decks: scheme verification set and the 32-grain pilots |
| `decks/polycrystal_plate/` | The CP plate template and the pipeline decks that drive `scripts/polycrystal_plate.py` |
| `scripts/` | `polycrystal_plate.py` (label → orient → deck → curve), `grain_morphology.py`, plotting, `tools/exodus_summary.py`, the material-point and backend-envelope harnesses |
| `meshes/plate/poly_d3_n32.npz` | The 32 Voronoi seeds of the pilot |
| `results/`, `figures/` | Outputs of the current pipeline only (field files are not committed) |
| `env/` | `pinnacle_syntax_dump.json` — the registry snapshot audited |
| `archive/` | The first campaign (Sept 2026): its scripts, decks, results, figures and two reports; see `archive/README.md` |

## Reports

| Document | Contents |
|---|---|
| `docs/PROJECT_GUIDELINES.md` | Project-wide rules |
| `docs/00_capability_audit.md` | Solver routes available; gap register G1–G27 |
| `docs/02_backend_envelope.md` | Measured backend/device policy |
| `docs/03_material_point_verification.md` | Constitutive law verified at one quadrature point; gap C.15 |
| `docs/05_orientation_notation.md` | Miller-direction and Euler conventions |
| `docs/06_polycrystal_plate_results.md` | Phase-field grain growth verified; the 30-grain plate and its first CPFE run |
| `docs/07_cp_parameter_study.md` | Lab notebook of the CP parameter study: base case, solver diagnosis, calibration, mesh sensitivity |
| `docs/reports/cp_fe_calibration_report.md` | Technical report (draft): CP FE methodology, numerical method, calibration results; figures from `scripts/report_figures.py` + `decks/report/` |
| `docs/workflow_polycrystal_plate.md` | How the plate pipeline works and why each setting |
| `docs/reproduce_polycrystal_plate.md` | Exact commands and environment on the ROCm host |
| `docs/reconciliation.md`, `docs/digitization_uncertainty.md` | Source extraction and its error |

## Status (2026-10-09)

The phase-field grain-growth route works under the fully implicit scheme, a
30-grain pseudo-3D plate exists, and a full crystal-plasticity run on it
completes from a PINNACLE deck plus two CSV companions. The reference CP set
overshoots the experiment by about 2x, so the next stage is calibration on
the x-loaded 100 x 50 x 1 um plate with 30 grains, the paper as ground truth.

Carried limits from Phase 0 (unchanged): 46 of 93 Al2219 orientations
recoverable, loading-axis only; Ni failure stress not recovered; gaps C.15
(Taylor-term density convention) and C.16 (Gc transferability) open.
