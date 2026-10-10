# Run ledger

One row per run (or run family launched together). Inputs are the deck and
the scripts that produced or consumed it; "rev" is the PINNACLE commit on
`feat/pinn_benchmark`. Host `rocm` = rocm-server, numpy backend, float64,
one physical core per run unless stated. Results and their interpretation
are in the numbered documents; this table only says how each run was made
and what it returned in one line.

| id | date | deck(s) | scripts | rev | host / cores | wall | outcome | record |
|---|---|---|---|---|---|---|---|---|
| gg-verify-D | 2026-10-09 | `decks/graingrowth/verify/D_M1_dt8.00_n1.i` | `pinnacle run` | 3e964c6^ | rocm / 1 | 9 s | convex_splitting, dt = 8/(L mu), one step: energy drop 0.04 | 06 §1 |
| gg-verify-ABC | 2026-10-09 | `.../verify/{A,B,C}_*.i` | `pinnacle run` x3 concurrent | 3e964c6^ | rocm / 1 each | ~3 min | convex_splitting: drop 29.5 / 16.7 / 16.7; B = C to 8 digits | 06 §1 |
| gg-verify-EFG | 2026-10-09 | `.../verify/{E,F,G}_bdf1_*.i` | `pinnacle run` x3 concurrent | 3e964c6^ | rocm / 1 each | 584 / 325 / 320 s | bdf1: drop 107.0 / 107.3 / 107.3; F = G to 8 digits | 06 §1 |
| gg-pilot-h40/80/160 | 2026-10-09 | `decks/graingrowth/pilot/n32_bdf1_h{40,80,160}.i` | `pinnacle run` x3 concurrent | 3e964c6^ | rocm / 1 each | 28 / 56 / 112 min | 30 / 30 / 27 grains; junction rms 14.3 / 14.1 / 15.3 deg | 06 §2 |
| pipe-smoke | 2026-10-09 | scratch deck from `plate_n32_h40.i` on the one-step n16 field, 2 steps | `polycrystal_plate.py label orient deck`, `pinnacle run` | 14dd45c then 141ebd8 | rocm / 1 | 10 s | before the translator fix: 47.8 / 95.7 MPa (= single crystal); after: 54.1 / 108.4 | 06 §4 |
| cpfe-n32-h40-y | 2026-10-09 | `decks/polycrystal_plate/plate_n32_h40.i` -> `results/polycrystal_plate/n32_h40/cpfe.i` | `polycrystal_plate.py`, `pinnacle run`, `plot_plate_curve.py` | 141ebd8 | rocm / 1 | 10 min 22 s | y-loaded, 40 steps to 3%, 0 cutbacks; 961.5 MPa at 3% | 06 §5 |
| moose-ref-csv | 2026-10-09 | `moose/CPFE_out.e` | `scripts/tools/exodus_summary.py --csv` | -- | rocm | 5 s | 40 steps to 3.0%, 808.8 MPa | 06 §3 |
| backend-4step | 2026-10-09 | scratch copies of `results/cp_study/h40_base/cpfe.i`, 4 steps | `pinnacle run` | 141ebd8 | rocm / 1 (CPU), GPU | 32 / 145 / 549 s | numpy / torch-CPU / torch-GPU, identical stresses | 06 §6 |
| backend-scaling | 2026-10-09 | same 4-step numpy deck at OMP 1/4/16 | `pinnacle run` | 141ebd8 | rocm | 32.2 / 32.5 / 32.4 s | no thread scaling; ~100% of one core each | 06 §6 |
| gpu-profile | 2026-10-09 | same deck, 1 step, torch-GPU | `python -m cProfile -m pinnacle run` | 141ebd8 | rocm / GPU | 130 s | 108 s in the per-point tangent (224k evaluations) | 06 §6 |
| cp-h40-base | 2026-10-09 | `decks/cp_study/h40_base.i` -> `results/cp_study/h40_base/cpfe.i` | `polycrystal_plate.py`, `pinnacle run` | 141ebd8 | rocm / 1 | 38 min | x-loaded; sub-step budget exhausted at step 100 (7.77%); 125/258 rows converged; curve valid to 7.8% | 07 §1 |
| cp-h20-base | 2026-10-09 | `decks/cp_study/h20_base.i` | `polycrystal_plate.py`, `pinnacle run` | 141ebd8 | rocm / 1 | 9 min | h20 twin: walls at step 126 (9.83%), 1,227 MPa; agrees with h40 to 3% at 5% strain | 07 §2 |
| cp-h20-diag | 2026-10-09 | copy of `results/cp_study/h20_base/cpfe.i` | `pinnacle run` | 6271f71 | rocm / 1 | 9 min | logs `CG did not converge (info=1386)` at step 126 on every cutback: the wall is the linear solver | 07 §2 |
| cp-h20-direct | 2026-10-09 | `results/cp_study/h20_base/cpfe.i` with `linear_method = direct` | `pinnacle run` | bb774d0 | rocm / 1 | 15 min | 256/256 steps to 20%, 0 cutbacks, 1,618 MPa: the wall was CG | 07 §2 |
| cp-h40-base (direct) | 2026-10-09 | `decks/cp_study/h40_base.i` (linear_method = direct) | `polycrystal_plate.py`, `pinnacle run` | 47a41e0 | rocm / 1 | (running) | base case to 20% | 07 |
| cp-h40-rescaled | 2026-10-09 | `decks/cp_study/h40_rescaled.i` | `polycrystal_plate.py`, `pinnacle run` | 6271f71+ | rocm / 1 | killed at 17 min (CG), relaunched with direct | paper's set in latent-sum form (rho0 13.41, k1 0.08475) | 07 §3 |
| taylor-screen | 2026-10-09 | `decks/cp_study/taylor_screen.i` | `taylor_surrogate.py oat`, `fit` | 47a41e0 | local / 12 | 4 s + ~30 s | surrogate optimum rho0 37.2, k1 0.126, y_c 0.00267, tau_c0 5.4; RMSE 13.2 | 07 §4a |
| cal-gn1 | 2026-10-09 | `decks/cp_study/cal_gn1.i` (5 cases) | `cp_parameter_study.py prepare, run` | 47a41e0 | rocm / 1 each, 5 concurrent | (running) | FE Gauss-Newton iteration 1 at h40 | 07 §4b |
| cp-oat-h40 | prepared, not launched | `decks/cp_study/oat_h40.i` | `cp_parameter_study.py prepare` | 141ebd8 | -- | -- | 9 cases + base written under `results/cp_study/` | -- |
