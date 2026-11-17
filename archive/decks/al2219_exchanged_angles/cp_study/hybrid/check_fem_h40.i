# NEURAL-ROUTE CHECK at h40 (author's mesh decision 2026-10-12) with the
# calibrated set of cal_center2.i, the study's solver criterion and field series.
# The cold FE reference is the calibrated run itself (results/cp_study/cal/gn2/base).
# Only the [Predictor] block and `corrector` differ between the variants.
# Variant: 8-step TIMING CHECK of mode = fem_residual at h40 (cost of one
# assembly epoch on the host, residual ratio before/after training); not a study result.

[Growth]
  run_dir = ../../../results/graingrowth_pilot/n32_bdf1_h40
  deck    = ../../graingrowth/pilot/n32_bdf1_h40.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  seed = 20261009
[]

[Mechanics]
  out_dir      = ../../../results/cp_study/hybrid/check_fem_h40
  load_axis    = x
  strain       = 0.00625
  t_end        = 6.25
  n_steps      = 8
  nx           = 40
  ny           = 20
  thickness    = 1.0
  newton_tol   = 1.0e-3
  newton_rel_tol = 0
  max_cutbacks = 10
  linear_method = direct
  corrector    = newton
  corrector_tol = 0.05
  snapshots_every = 0
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
  y_c     = 0.00343758
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
  epochs_initial  = 40
  epochs_per_step = 20
  lr              = 2.0e-3
  optimizer       = adam
  residual_target = 0.0
  reset_each_step = false
  seed            = 0
[]
