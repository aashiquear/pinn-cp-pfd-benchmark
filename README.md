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
| `bib/references.bib` | The single central bibliography, 184 entries, keys `author-year-significantword` |
| `bib/citation_keys.md` | Key → full citation → prose form, grouped by topic; verification flags |
| `data/extracted/` | `benchmark_dataset.csv` (314 provenance-tagged rows), `gap_register.csv` (C.1–C.16) |
| `data/digitized/` | Per-figure extracted data, one CSV per panel |
| `data/reference/` | Reference simulation and experimental data |
| `decks/al2219/`, `decks/ni/` | Study input decks |
| `meshes/` | Generated meshes and generation scripts |
| `env/` | `pinnacle_syntax_dump.json` — the 249-component registry snapshot audited |
| `docs/` | Reports, listed below |
| `results/`, `figures/`, `runs/` | Outputs |
| `scripts/` | Harnesses and analysis code |

## Reports

| Document | Contents |
|---|---|
| `docs/PROJECT_GUIDELINES.md` | Project-wide rules; includes the project-description block |
| `docs/00_capability_audit.md` | Solver routes available; gap register G1–G11 |
| `docs/01_bundled_deck_assessment.md` | Why the bundled Rezwan decks are a demonstration, not a reproduction |
| `docs/02_backend_envelope.md` | Measured backend/device policy (research question R2) |
| `docs/reconciliation.md` | Independent extraction vs the companion data sheet, discrepancy by discrepancy |
| `docs/digitization_uncertainty.md` | Per-figure digitization method, calibration, and error |

## Status

**Phase 0 complete.** Sources extracted and reconciled, bibliography built,
framework audited, backend policy measured.

Carried limits, to be restated wherever they affect a result:

- **46 of 93** Al2219 orientations are recoverable, and only as loading-axis
  directions (two of three rotational degrees of freedom). Every IPF-based
  comparison reports `n = 46` and declares the third-angle convention it fixed.
- Simulated failure **stress** per Ni specimen is not recovered; failure
  **strain** is.
- Six blocking framework gaps (G1–G6) and G11 must close before a faithful
  reproduction is possible. G1 (no plastic-work term in the fracture driving
  force) and G2 (Taylor prefactor α₀ absent) would otherwise produce wrong
  answers silently rather than errors.
- Gap **C.15** — whether the Taylor term uses the per-system self density
  (τ_c(0) = 43.72 MPa) or the latent-hardened sum over 12 systems
  (176.60 MPa, a factor √(1+11q) = 4.05) — must be settled against the digitized
  calibration curve before any calibration run. PINNACLE implements the latent
  sum; the paper's own reported value corresponds to the self term.
- Gap **C.16** — the reported Gc differs between the two materials by a factor
  of 5,495 (22.75 J/m² for Al2219 against 1.25 × 10⁵ J/m² for Ni). The Al2219
  value is far below the range expected for a metal, so it is treated as an
  effective regularized value pending sensitivity analysis, not as a
  transferable material property.
