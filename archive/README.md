# Archive

Nothing here is run by the current workflow. Three generations are kept for
provenance.

## Superseded versions of the two reports (2026-11)

| Path | What it holds |
|---|---|
| `docs/al2219_report_exchanged_angles_2026-11-16.md`, `figures/al2219_exchanged_angles/` | The Al2219 report before the recheck of 2026-11-16: Euler angles read with phi1 and phi2 exchanged, the earlier balance of linear momentum, the explicit constitutive update |
| `decks/al2219_exchanged_angles/` | Its decks, moved from `decks/cp_study/`, `decks/pfd_study/cp/`, `decks/pfd_study/report/`, `decks/benchmark/` and `decks/report/`. Relative paths inside them refer to the old locations |
| `docs/ni_report_first_attempt_2026-11-15.md`, `figures/ni_first_attempt/` | The Ni report before the restart of 2026-11-15 |
| `decks/ni_exchanged_angles/`, `figures/ni_exchanged_angles/` | The Ni decks and figures before the Bunge rerun of 2026-11-16 |
| `docs/notebooks/` | Lab notebooks and plans of those campaigns: `00_capability_audit.md`, `06_polycrystal_plate_results.md`, `07_cp_parameter_study.md`, `08_pinn_route_assessment.md`, `09_pfd_benchmark_plan.md`, `10_pfd_study.md` (the J2 notes of report §13 are here), `11_ni_plan.md`, `12_ni_study.md`, and the reproduction document of the first Al2219 study |
| `scripts/ship_results.py`, `results/SHIPPED.md` | The result-shipping tool and its manifest, retired by the shipping rule of 2026-11-16 (no solver output in the repository) |
| `scripts/ni_diag_newton.py`, `scripts/ni_pfd_search_figures.py` | Ni diagnostics of the first attempt |

## First campaign (2026-09)

Material from the study's first campaign (2026-09-08 to 2026-09-14),
superseded by the phase-field polycrystal plate pipeline of 2026-10-09. It is
kept for provenance and for the findings its documents record; nothing here
is run by the current workflow. Binary outputs (Exodus meshes, VTK fields,
pickles) and files over 2 MB were deleted when it was archived.

| Folder | What it holds |
|---|---|
| `docs/` | `01_bundled_deck_assessment.md` (why the bundled Rezwan decks are a demonstration), `04_seventeen_percent_blocker.md` (the stress-concentrator sub-step wall and the negative parameter sweeps) |
| `scripts/` | The earlier harnesses: void/plate/polycrystal mesh generators (hexaform), `run_plate.py` / `run_polycrystal.py` (coupled CP-PFD plate runs with injected orientations), `run_convergence_level.py` and the Modal job scripts, `calibrate_cp_polycrystal.py` (Taylor material-point fit), `grow_grains.py` / `round_voronoi.py` / `map_grains.py` / `plot_morphology.py` (the convex-splitting growth attempt and the threshold-dynamics substitute), diagnostics |
| `decks/al2219/` | Void-cube and plate templates with placeholders, isochoric BC deck |
| `meshes/` | Void-cube and plate mesh metadata (`.npz`, family JSON); the `.e` files were deleted |
| `results/` | Every earlier run family: convergence, calibration, Gc sweeps, BC comparisons, rate and rho0 sweeps, the schedule/prediction study that suspended the FE plasticity numbers, the convex-splitting grain-growth runs and the rounded (threshold-dynamics) structures |
| `figures/` | Figures of those runs |

`docs/03_material_point_verification.md` and the archived notebooks still cite
some of these scripts by their old `scripts/` path; read `archive/scripts/`.

Archived 2026-11-17 with the last cleanup: the three `polycrystal_plate_grains_*.png`
figures and `plot_grain_structures.py` of notebook 06, the Phase 0 syntax dump
(`env/pinnacle_syntax_dump.json`, notebook 00) and the seed file
`meshes/plate/poly_d3_n32.npz`. No current document, deck or script uses them.
