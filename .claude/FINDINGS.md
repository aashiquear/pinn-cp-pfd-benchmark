# Findings and working rules for agents on this study (as of 2026-10-09)

Read `README.md` and `PROJECT_GUIDELINES.md` first. This file holds what is
NOT derivable from the code: measured facts, decisions, and traps.

## Decisions by the author (Rezwan)
- Study plate: 100 x 50 x 1 um, ONE element through thickness, loaded along
  **x** (the 100 um length). Mesh level hN = N elements along x, N/2 along y.
- 30 grains for every mesh, from phase-field grain growth of a Voronoi seeding.
- Ground truth is the paper (Rezwan 2025 intragranular) and its digitized
  experiment markers; the MOOSE decks in `moose/` are initial tests only.
- Runs must reach >= 20% strain.
- Fixed CP parameters: C11, C12, C44, b, mu, q, alpha0, eta, c_thr. gamma0 and
  m (= 1/n) change only for a convergence reason and only with the author.
  Free, fit in this order: rho0, k_m (k1), y_c, tau_c0. Report what each does
  to the curve, then fit the experimental shape, then error analysis.
- Taylor term uses the LATENT-HARDENED SUM (PINNACLE's convention), to be
  discussed as a model limitation later.
- Phase-field (damage) parameters are revisited after CPFE is settled.
- Grain-instance check: a 30-grain pseudo-3D realization can be a severe
  outlier (one was 17.5% stronger than the calibrated structure, mostly from
  grain arrangement). The author drops such outliers and draws another seed;
  an instance within a few percent is acceptable.
- Compute: CPU for now. The author wants GPU eventually; see below.
- Never launch a multi-run campaign without (1) one checked single case and
  (2) the author's agreement on the design.
- Every run is reproducible from a PINNACLE `.i` deck; pipeline scripts read a
  deck. Workflow, reproduction steps and results are separate documents.
  Each run is logged in `docs/run_ledger.md`.

## Measured facts
- Phase-field grain growth: use `scheme = bdf1`. `convex_splitting` is
  lag-dominated for multi-order-parameter energies (stabiliser 147.5 mu for
  8 components) and makes mobility look inert. dt = 0.5/(L mu) is converged
  to 0.2%. Junctions equilibrate by 2.5 interface times; longer only rounds
  grains and consumes small ones. Details: `docs/06_polycrystal_plate_results.md` §1-2.
- CP plate run cost (numpy): ~8 s per load step on 800 elements; one core per
  run regardless of OMP threads (Python-bound). Host has 16 physical cores:
  run cases concurrently, ~14 at a time.
- torch backend on the same deck: CPU 4.5x slower, GPU (Radeon 8060S) 17x
  slower than numpy. Cause: per-point single evaluations in the tangent
  (224k calls per step). The APU's fp64 GPU throughput equals its CPU's.
  Working wheel: `torch==2.13.0+rocm7.2` (the rocm7.1 wheel segfaults).
- With the reference parameter set, PINNACLE's 30-grain plate (y-loaded,
  3%) gives 961 MPa at 3% and MOOSE's single crystal 809; the experiment is
  ~412 MPa at 2%. Calibration is required; parameters are framework-sensitive.

- C.15 is a 4.05x effect on tau_c: the paper's Table 2 is a SELF-term set;
  under PINNACLE's latent-hardened sum the same numbers give a 2,340 MPa
  plateau (no saturation before ~28% shear). The latent-sum equivalent is
  rho0 = 13.41, k1 = 0.08475 (y_c, tau_c0 unchanged). Neither reading
  matches the experiment at 1e-3 /s; calibration from the rescaled set.
  Rate sensitivity (gamma_dot0 = 1e-5, n = 8.33) multiplies the flow stress
  by ~2 at 1e-3 /s, so the assumed test rate matters.
- The explicit CP sub-stepper walls at 7.8% (h40) / 9.8% (h20) on the plate
  with the reference set; cause being diagnosed (PINNACLE now logs the
  Newton or budget reason per cutback). Curves are valid only to the first
  failed step; scripts/polycrystal_plate.py truncates them.

## Traps found
- `[OrientationField]` in a deck was silently dropped by the translator until
  PINNACLE 141ebd8 (fixed). Always check `pinnacle.log` for the line
  "[OrientationField] ... resolved on N elements".
- A FLAT `config["material"]` dict in CrystalPlasticity takes SI defaults in a
  um/MPa deck and produces stresses ~1e10 MPa silently (earlier campaign).
- Mid-run counters (cutbacks, steps) are not results; only a written
  `solver_history.csv` is.
- MOOSE deck `moose/CPFE.i` has y_c = 0.00286; the paper table says 0.00268.

## Technical report (docs/reports/)
- `docs/reports/cp_fe_calibration_report.md` is the paper-style record; the
  numbered `docs/0N_*.md` stay the lab notebook. The report carries ONLY the
  final condition: no tests, no verification, no solver diagnosis, unless the
  author decides otherwise. Abstract and introduction stay blank until all
  studies are done; mesh sensitivity is held out for now.
- Report figure rules (author, 2026-10-10): prepend (0, 0) to every
  stress-strain curve; parameters as math symbols; no grid; closed axes box;
  box from (0, 0) to one tick block above the largest stress (700 MPa) and
  20.5% strain. All of it is in `decks/report/cp_calibration_figures.i`,
  rendered by `scripts/report_figures.py` into `figures/report/`.
- Mesh sensitivity with the calibrated set (2026-10-11): 20% stress 636.4 /
  603.9 / 586.4 / 573.5 MPa at h20/h40/h80/h160; change per halving -5.1 /
  -2.9 / -2.2%, ratio rising to 0.74 at 20% (order ~0.44): slip localization
  without an internal length keeps the plastic branch mesh-dependent. The
  converged-mesh choice is the author's; h160 costs 17.7 h per run.
- h40 has 1,722 nodes / 5,166 displacement DOFs (41 x 21 x 2 nodes); an
  earlier note of 7,749 was wrong.

## Where things are
- Host: `ssh rocm-server`, study mirror `~/projects/pinn-cp-pfd-benchmark`
  (rsync, not git), PINNACLE `~/projects/pinnacle` on `feat/pinn_benchmark`,
  venv `~/venvs/pinnacle`. Use `OMP_NUM_THREADS=1` per run.
- Pipeline: `scripts/polycrystal_plate.py <pipeline deck> label orient deck | curve`;
  study driver `scripts/cp_parameter_study.py <study deck> prepare|run|collect`.
- Earlier campaign (Sept 2026) is in `archive/`; its reports still hold
  valid findings (G20 strain-increment limit 8e-4 per step, G25 cutback
  behaviour, the stress-concentrator sub-step wall).
