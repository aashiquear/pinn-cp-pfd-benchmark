# GRAIN-INSTANCE CHECK, instance B: the calibrated h40 set (cal_center2.i) on a
# second grain structure -- a new Voronoi seeding (rng seed 20261012, deck
# decks/graingrowth/instances/n32_bdf1_h40_sB.i) grown with the same
# phase-field settings, and a new orientation draw (seed 20261012). Everything
# else (mesh h40, loading, parameters) is the calibrated run's. The curve is
# compared with the calibrated structure's to show the parameter set is not
# tied to one grain realization.
[Growth]
  run_dir = ../../results/graingrowth_instance/sB
  deck    = ../graingrowth/instances/n32_bdf1_h40_sB.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  seed = 20261012
[]

[Mechanics]
  out_dir      = ../../results/cp_study/grain_instance/sB
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
  snapshots_every = 4
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
