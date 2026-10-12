# HYBRID FIRST CHECK (2026-10-11), h40 plate, calibrated-set parameters of
# cal_center2.i, PINNACLE a54715a (ADR 0027 batched CP path, ADR 0028
# load-step predictor services). Same plate, structure, orientations, loading
# and solver as the FE run; the [Predictor] block below is the only
# difference. The FE arm runs on numpy (one core); the network trains on its
# own torch backend (GPU when present). Compared against the cold run of the
# same deck (h40_cold.i) step by step: Newton iterations, newton_time_s,
# predictor_time_s, predictor_rel_error, and the training cost.
#
# Variant: PINN AS SOLVER -- the physics-trained network (frozen quadrature-point
# state, first Piola-Kirchhoff residual, soft Dirichlet, no data term) IS the
# step solver: corrector = none; its equilibrium residual norm is recorded.

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
  out_dir      = ../../../results/cp_study/hybrid/h40_pinn_solver
  load_axis    = x
  strain       = 0.20
  t_end        = 200.0
  n_steps      = 256
  nx           = 40
  ny           = 20
  thickness    = 1.0
  newton_tol   = 1.0e-6
  max_cutbacks = 10
  linear_method = direct
  snapshots_every = 4
  corrector    = none
  C11     = 108.0e3
  C12     = 61.3e3
  C44     = 28.5e3
  mu      = 25.7e3
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 5.34925
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 41.1521
  k1      = 0.16771
  y_c     = 0.00319744
[]

[Predictor]
  kind            = load_step_predictor
  mode            = physics
  warmup_steps    = 2
  backend         = torch
  device          = auto
  use_gpu         = true
  precision       = float64
  hidden_layers   = 3
  neurons         = 64
  activation      = tanh
  fourier_features = false
  epochs_initial  = 1500
  epochs_per_step = 300
  lr              = 2.0e-3
  window          = 0
  seed            = 0
  w_physics       = 1.0
  w_dirichlet     = 10.0
  w_data          = 0.0
  physics_scale   = 0.0
[]
