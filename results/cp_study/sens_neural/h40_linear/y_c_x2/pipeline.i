# Case 'y_c_x2' of study sens_h40_linear.i: {'y_c': 0.00687516}
# NEURAL-ROUTE CHECK at h40 (author's mesh decision 2026-10-12) with the
# calibrated set of cal_center2.i, the study's solver criterion and field series.
# The cold FE reference is the calibrated run itself (results/cp_study/cal/gn2/base).
# Only the [Predictor] block and `corrector` differ between the variants.
# Variant: PINN+FEM reference predictor: linear extrapolation of the last two accepted steps.

[Growth]
  run_dir = /home/aashiquear/projects/pinn-cp-pfd-benchmark/results/graingrowth_pilot/n32_bdf1_h40
  deck    = /home/aashiquear/projects/pinn-cp-pfd-benchmark/decks/graingrowth/pilot/n32_bdf1_h40.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  seed = 20261009
[]

[Mechanics]
  out_dir      = /home/aashiquear/projects/pinn-cp-pfd-benchmark/results/cp_study/sens_neural/h40_linear/y_c_x2
  load_axis    = x
  strain       = 0.20
  t_end        = 200.0
  n_steps      = 256
  nx           = 40
  ny           = 20
  thickness    = 1.0
  newton_tol   = 1.0e-3
  newton_rel_tol = 0
  max_cutbacks = 10
  linear_method = direct
  corrector    = newton
  snapshots_every = 4
  C11     = 108.0e3
  C12     = 61.3e3
  C44     = 28.5e3
  mu      = 25.7e3
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 5.29147
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 47.9872
  k1      = 0.198833
  y_c     = 0.00687516
[]

[Predictor]
  kind         = load_step_predictor
  mode         = linear
  warmup_steps = 2
[]
