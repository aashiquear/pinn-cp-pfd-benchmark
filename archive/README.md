# Archive

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

Current documents (`docs/00_capability_audit.md`, `docs/03_material_point_verification.md`)
still cite some of these scripts by their old `scripts/` path; read `archive/scripts/`.
