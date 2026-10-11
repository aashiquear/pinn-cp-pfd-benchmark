# FRAMEWORK-FIX CHECK (2026-10-11): the calibrated h40 run (cal_center2.i)
# repeated with PINNACLE a311a8b, whose batched sub-stepper now passes the
# per-point rotated stiffness C to the state update (ADR 0027). Before the
# fix every orientation-field run evolved its CP state under the material-
# global elastic stiffness while the residual used the rotated one. This run
# measures the effect on the reported curve; everything else is identical.
[Growth]
  run_dir = ../../results/graingrowth_pilot/n32_bdf1_h40
  deck    = ../graingrowth/pilot/n32_bdf1_h40.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  seed = 20261009
[]

[Mechanics]
  out_dir      = ../../results/cp_study/checks/cfix_h40
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
