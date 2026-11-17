# NEURAL-ROUTE CHECK at h40 (author's mesh decision 2026-10-12) with the
# recalibrated set of ../cal_center2.i (Bunge angles, total-Lagrangian balance,
# implicit update, torch backend on the CPU; network settings unchanged), the study's solver criterion and field series.
# The cold FE reference is the calibrated run itself (results/cp_study/bunge/cal/gn2/base).
# Only the [Predictor] block and `corrector` differ between the variants.
# Variant: PINN+FEM reference predictor: linear extrapolation of the last two accepted steps.

[Growth]
  run_dir = ../../../../results/graingrowth_pilot/n32_bdf1_h40
  deck    = ../../../graingrowth/pilot/n32_bdf1_h40.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  seed = 20261009
[]

[Mechanics]
  out_dir      = ../../../../results/cp_study/bunge/hybrid/h40_linear
  backend      = torch
  kinematics   = total_lagrangian
  state_update = implicit
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
  tau_c0  = 5.38956
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 53.9306
  k1      = 0.176771
  y_c     = 0.00271423
[]

[Predictor]
  kind         = load_step_predictor
  mode         = linear
  warmup_steps = 2
[]
