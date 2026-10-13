# Case 'y_c_x0p8' of study sens_h40_grain_solver.i: {'y_c': 0.00275006}
# NEURAL-ROUTE CHECK at h40 (author's mesh decision 2026-10-12) with the
# calibrated set of cal_center2.i, the study's solver criterion and field series.
# The cold FE reference is the calibrated run itself (results/cp_study/cal/gn2/base).
# Only the [Predictor] block and `corrector` differ between the variants.
# Variant (GRAIN-AWARE INPUTS, node_features = orientation grain_boundary): PINN as solver, verified (PINNACLE ADR 0029): mode = fem_residual with
# residual_target = corrector_tol and corrector = auto. Training stops once the
# residual over the free degrees of freedom is below corrector_tol of the
# reaction; the step is then accepted without Newton (corrector_used = none) or
# corrected by Newton when the epoch budget runs out (corrector_used = newton).

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
  out_dir      = /home/aashiquear/projects/pinn-cp-pfd-benchmark/results/cp_study/sens_neural/h40_grain_solver/y_c_x0p8
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
  corrector    = auto
  corrector_tol = 0.05
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
  y_c     = 0.00275006
[]

[Predictor]
  kind            = load_step_predictor
  mode            = fem_residual
  baseline        = linear
  warmup_steps    = 2
  backend         = torch
  device          = auto
  use_gpu         = true
  precision       = float64
  hidden_layers   = 3
  neurons         = 64
  activation      = tanh
  fourier_features = false
  epochs_initial  = 5
  epochs_per_step = 4
  lr              = 1.0
  optimizer       = lbfgs
  residual_target = 0.05
  reset_each_step = false
  seed            = 0
  node_features   = 'orientation grain_boundary'
[]
