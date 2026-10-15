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
- FIELD VARIABLES ARE PART OF EVERY RESULT (author, 2026-10-12): every run
  (FE, PINN, PINN+FEM) writes the field series (`snapshots_every = 4` in the
  pipeline deck, i.e. every 0.3125% strain) so fields can be compared
  between routes, not only the stress-strain curve.
- Newton criterion of every run: `newton_rel_tol = 0`, `newton_tol = 1e-3`
  (docs/07 §10.0): the relative-to-first-residual test penalises a warm
  start, so cold and predicted starts must share one absolute criterion.
- Results are produced with the fixed solver only (PINNACLE a311a8b or
  later); nothing computed before it is kept in results/ or in the report.
- Never launch a multi-run campaign without (1) one checked single case and
  (2) the author's agreement on the design.
- Every run is reproducible from a PINNACLE `.i` deck; pipeline scripts read a
  deck. Workflow, reproduction steps and results are separate documents.
  Each run is logged in `docs/run_ledger.md`.

## Part II (PFD, CP-PFD) decisions by the author (2026-10-13)
- Order: PFD on J2 first, then CP-PFD; same strategy as Part I (calibration,
  mesh sensitivity, parameter sensitivity, PINN+FEM, PINN, surrogate,
  hardware). Same plate and mesh resolution, central hole of diameter 5 um
  (Irwin plastic zone clear of the boundaries). Polycrystal again (Part I
  grains and point file unchanged).
- Meshes from hexaform, hexahedral only, shown to the author before any run.
  hexaform `plate` (branch feat/plate-hole) is the generator;
  `scripts/plate_hole_mesh.py` + `decks/pfd_study/mesh/*.i`.
- Calibrate G_c and l, eta fixed. Fracture initiation = d >= 0.8 at any node
  or element average; fracture strain target = last point of the Vaughan
  curve (17.2% nominal). Mesh sensitivity reproduces (chakraborty 2016
  multiscale) Fig. 8 (curve, peak stress, fracture energy vs l); the brittle
  l > 2h rule is untested for J2/CP coupling and is what the study measures.
- Shape study after calibration and mesh sensitivity (plan in docs/09 §6).
- One report for Part I and II: `docs/reports/cp_pfd_pinnacle_benchmark_report.md`
  (renamed 2026-10-13); final results only; abstract and introduction empty.
- Hardware: best available for every run except the hardware comparison =
  torch CPU arm (see Measured facts); GPU for network training.
- Plan and PINNACLE gap list: `docs/09_pfd_benchmark_plan.md`.
- J2+PFD configuration (author, 2026-10-14): AT2 with the plastic-work
  threshold W_c = 20 MPa (restores the peak stress without moving the
  initiation strain); the W_c = 0 pass stays in docs/10 as the comparison.
  CP+PFD keeps W_c = 0 (no early softening there) unless the author says
  otherwise. The first CP calibration (Gc 519/847/1302) is VOID (coupled
  window dt bug, PINNACLE f3eeae7); the corrected one gives Gc of the J2 order.

## Measured facts
- Phase-field grain growth: use `scheme = bdf1`. `convex_splitting` is
  lag-dominated for multi-order-parameter energies (stabiliser 147.5 mu for
  8 components) and makes mobility look inert. dt = 0.5/(L mu) is converged
  to 0.2%. Junctions equilibrate by 2.5 interface times; longer only rounds
  grains and consumes small ones. Details: `docs/06_polycrystal_plate_results.md` §1-2.
- CP plate run cost (numpy): ~8 s per load step on 800 elements; one core per
  run regardless of OMP threads (Python-bound). Host has 16 physical cores:
  run cases concurrently, ~14 at a time.
- torch backend on the same deck: BEFORE ADR 0027, CPU 4.5x slower, GPU
  (Radeon 8060S) 17x slower than numpy. AFTER (batched, stacked state,
  2026-10-11): 4-step h40 plate numpy 6.0 s, torch-CPU 4.6 s, torch-GPU
  4.55 s, identical results; CP field solves may run on torch at no penalty. Cause: per-point single evaluations in the tangent
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
- Killing by pattern (`pkill -f`, `pgrep -f`) matches the shell that runs the
  kill when the SAME command line also contains the pattern text (an ssh
  remote string, or a launch of the same program later in the command): the
  job kills itself (exit 144/255). Kill by PID from a separate command, and
  never combine a pattern kill with a launch of the same program.
- PINNACLE before a311a8b: on an ORIENTATION-FIELD deck the batched CP
  sub-stepper passed the per-point slip tensors but not the per-point
  rotated stiffness C, so the state evolved under the material-global
  elastic law while the residual used the rotated one (single crystals
  unaffected). Fixed in a311a8b (ADR 0027). MEASURED EFFECT (docs/07 §9): the
  calibrated h40 curve is 5-10% lower with the fix (545.6 vs 603.9 MPa at
  20%); every pre-fix polycrystal result (calibration, mesh, sensitivity,
  instances, field series) is superseded and the campaign must be rerun on
  the fixed framework (decks unchanged; runs now 4-5x faster).
- `[OrientationField]` in a deck was silently dropped by the translator until
  PINNACLE 141ebd8 (fixed). Always check `pinnacle.log` for the line
  "[OrientationField] ... resolved on N elements".
- A FLAT `config["material"]` dict in CrystalPlasticity takes SI defaults in a
  um/MPa deck and produces stresses ~1e10 MPa silently (earlier campaign).
- Mid-run counters (cutbacks, steps) are not results; only a written
  `solver_history.csv` is.
- MOOSE deck `moose/CPFE.i` has y_c = 0.00286; the paper table says 0.00268.

## PINN / hybrid routes for CP (assessed 2026-10-11, docs/08)
- Neither route can run the CP plate today: the PINN residual for
  stress_divergence_cp is displacement-only at the virgin state (elastic) and
  ignores [OrientationField]; pinn_fem_hybrid's corrector never updates the
  CP state, hardcodes finite_strain = False and does not ramp Dirichlet
  values. ADR 0017: a coordinates-in network cannot hold CP history. Plan in
  docs/08: GPU-resident batched CP evaluation first, then a load-path hybrid
  built on cp_quasi_static, then a per-step "PINN as solver".

- Status 2026-10-12: both routes RUN on the plate (PINNACLE ADR 0027-0029,
  `cp_quasi_static` + `[Utilities] load_step_predictor`). Measured at h40:
  linear extrapolation saves 37% of Newton iterations under the absolute
  criterion; strong-form collocation (`mode = physics`) is unreachable on a
  30-grain plate (network within 0.55% L2 of the FE field has FE residual
  0.73 of the reaction); FE-residual training (`mode = fem_residual`, exact
  consistent-tangent gradient, 0.28 s per assembly epoch) improves the
  extrapolation's residual ratio by only 20-30% per plastic step even with
  450 L-BFGS assemblies per step and never reaches the 0.05 tolerance in the
  plastic regime: the global smooth network of (x, lambda) is the limit, not
  the loss or the optimizer. Adam lr 1e-2 diverges in this mode; use L-BFGS.
  Full h40 runs (docs/07 §10.8, report §8): hybrid 945 iterations vs 948
  with the extrapolation alone (no gain, +349 s); verified solver 132/256
  steps without Newton, curve within 0.85%, fields within 0.6%, 6-8x the
  FE cost. Grain-aware inputs (`node_features`, 9ae5ff3): residual 3-4x
  lower per step, solver fields 15x closer (curve within 0.05%), same
  acceptance (131/256), cost still 5-6x FE; hybrid 923 iterations.
- Trap: `predictor_residual_ratio_initial` on later steps is the baseline
  PLUS the carried-over correction, not the pure extrapolation; the pure
  extrapolation's ratio is on the `h40_linear`-style run (guard columns).

- Surrogate route (PINNACLE ADR 0030, 2026-10-13): physics-only training works
  on conduction (0.3%, certified) but NOT on the 2-D elastic plate (70-90%)
  despite traction term, physical scales, 512 collocation points; open item.
  CP surrogates are data-only by construction (history dependent).

- Plate surrogate (2026-10-13): 24 LHS FE samples (623 s each); curve surrogate
  held-out 2-3% (6.6% max, data-limited), field surrogate 1%; training 28-121 s;
  0.2-0.4 ms per query (CPU faster than GPU at this network size); break-even
  14 FE runs. Report §11.

- Hardware benchmark h40 (2026-10-13): the numpy CP assembly ignores threads
  (0.45 s/iteration at 1-8 threads); torch CPU is 6.7x faster per iteration
  (0.067 s) and torch GPU 5.7x (0.079 s) at 7-9x memory. Use the torch CPU arm
  for future absolute timings; route ratios are unaffected. The solver route's
  cost is its training assemblies (42/step), not the network's device.
  Thread scaling of the torch CPU arm (2026-10-13): 0.142 / 0.092 / 0.076 /
  0.067 s per iteration at 1 / 2 / 4 / 8 threads, 1.2 GB each; campaigns run
  seven cases at two threads side by side. The direct linear solver
  (2026-10-09) is unrelated to this gain (it removed the CG wall on all arms).
- PINNACLE PFD readiness (ADR 0031, 2026-10-13, docs/10 §1): sparse damage
  solve (identical to dense to 4e-15; h80 0.05 s / 24 MB vs 4.3 s / 1.8 GB),
  `boundary_node_set = none` and strict named sets, `pfd_fracture viscosity`,
  `pfd_driving_force plastic_work_threshold`, `reaction_node_set` on both
  quasi-static executioners, per-step staggered rows with `max_damage` /
  `nodes_above_threshold` / `initiated`, netCDF Exodus reader (meshio had
  misaligned the hexaform node sets), analytic J2 split tangent (the FD one
  doubled the diagonal columns), einsum contraction path, `line_search =
  backtracking` (a full Newton step diverges on a degraded state) and
  `fem_quasi_static max_cutbacks`.
- Newton tolerances on the plate: use Part I's criterion (newton_tol 1e-3
  absolute, rel 0). An absolute tolerance of 1e-6 or PINNACLE's 1e-10 default
  sits below the force-scale noise floor (forces ~1e4 MPa um^2) and reports a
  false "did not converge"; this cost an hour of diagnosis.
- Variables default to `order = 1`; set `order = 2` on u and d in every plate
  deck (the templates do).
- J2 law for Part II: the piecewise table of `results/pfd_study/j2_fit/j2_fit.json`
  (within 0.3% of the markers to 17%); the power law misses by 6% at 17%.

## Technical report (docs/reports/)
- `docs/reports/cp_pfd_pinnacle_benchmark_report.md` is the paper-style record; the
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

- 2026-10-14: `hardware_benchmark.py` stages each case from the `out_dir` of the
  pipeline deck it names; a `[Plate]` deck that was never run by itself has no
  staged `pfd.i` there (the study decks stage their own copies under
  `results/pfd_study/<study>/<case>/`). Point hardware cases at a staged
  `pipeline.i` (`neural_h40/base/pipeline.i`), not at the calibration deck.
  The coupled routes' Newton counts to compare are `first_mech_iterations`
  (what a predictor replaces) and `mech_iterations_total`; the plain
  `mech_iterations` column is the last alternation only.

- 2026-10-15 (author's revision round, memory `study-revision-2026-10-14`): the J2+PFD
  configuration is W_c = 200 MPa, Gc = 385.5 MPa um, l0 = 5 um (W_c scan, docs/10 §4.10);
  every J2 result at W_c 0 or 20 is superseded (kept as the calibration record). CP+PFD
  unchanged (Gc 1106, l0 5, W_c 0). Sensitivities carry four values per parameter; the
  mesh studies carry h20/h40/h80; the shape study (docs/10 §9) ran all three series on
  CP+PFD. Figure rules: colour per parameter and line style per value in categorical
  studies, colour bars on field maps, crack maps = equivalent plastic strain with the
  d >= 0.8 elements solid, a GIF of initiation/propagation (`pfd_fields.py gif = true`).
- The J2 VTK series carry `equivalent_plastic_strain` only from PINNACLE 863b787 (ADR 0031
  item 16); older J2 series (W_c 0/20 calibrations) do not. Local PINNACLE runs need
  netCDF4 in the python3.11 user site (`python3.11 -m pip install --user netCDF4`,
  done 2026-10-15); the Exodus reader refuses to run without it.
- Study-deck case names are lowercased by the translator (`A_x2.5_y5` -> `a_x2.5_y5`):
  result directories and `[Labels]` keys follow the lowercase form.
- 2026-10-15: the author keeps the J2+PFD benchmark at h40 despite the 1.2% mesh
  shift of the thresholded source (side-by-side comparison with the CP h40 cases
  takes precedence); reported as a limitation in report §13.5. The study is pushed to
  `origin/main` (github.com/aashiquear/pinn-cp-pfd-benchmark, public); local branch
  `master` tracks `origin/main`; git uses `gh` as the credential helper.
- 2026-10-15: shipping policy reverted by the author: the result files the report
  presents ARE in git (field series behind the figures and the GIFs, meshes, networks;
  `scripts/ship_results.py decks/report/ship_results.i`, manifest `results/SHIPPED.md`,
  239 MB). Everything else under results/ stays summaries only (csv/json/logs); the
  superseded J2 W_c 0/20 sweep and sensitivity outputs, the check/probe runs and the
  hybrid check_fem diagnostics were removed from git and disk. When a new result enters
  the report, add its series to the ship deck and rerun the script before pushing.
