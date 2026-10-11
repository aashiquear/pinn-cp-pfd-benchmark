# TEXTURE / GEOMETRY SWAP: the calibrated grain geometry (30 grains) with the
# first 30 orientations of grain instance A (results/cp_study/grain_instance/sA/
# grains.csv), everything else the calibrated run's. Instance A is 17.5%
# stronger than the calibrated structure while its area-weighted Taylor
# estimate is only 4.5% stronger; this run says how much of the offset the
# orientation draw alone carries on a fixed geometry.
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
  # Texture/geometry swap: instance A's orientations on the calibrated geometry
  file = ../../results/cp_study/grain_instance/sA/grains.csv
[]

[Mechanics]
  out_dir      = ../../results/cp_study/grain_instance/swap_calgeom_sAori
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
