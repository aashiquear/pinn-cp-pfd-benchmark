# Case 'tau_c0_x0p8' of study sens_h40.i: {'tau_c0': 4.23318}
# CP calibration, FE Gauss-Newton CENTER, iteration 1.
#
# Plate 100 x 50 x 1 um loaded along x, mesh level h40 = 40 x 20 x 1 HEX8
# (2.5 um elements), the 30-grain structure grown for 2.5 interface times
# (results/graingrowth_pilot/n32_bdf1_h40), orientation seed 20261009.
# Loading to 20% nominal strain at 1e-3 /s (t_end = 200 s) in 256 steps
# (7.8e-4 strain per step, under the G20 limit of 8e-4). The field series is
# written every 4 steps (snapshots_every = 4) for field comparisons.
#
# STRAIN RATE is an assumption: Vaughan (2024) does not report it; 1e-3 /s is
# the quasi-static value. With n = 8.33 a factor 3 in rate moves the stress 14%.
#
# Fixed by decision: C11 C12 C44 b mu q alpha0 (and gamma0, n unless a
# convergence issue forces it). Free, in study order: rho0, k1 (k_m), y_c,
# tau_c0. Starting values: rho0 41.15 um^-2, k1 0.1677, y_c 0.0032 um,
# tau_c0 5.35 MPa. Iteration 2 center: the Levenberg-Marquardt step of iteration 1.

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
  out_dir      = /home/aashiquear/projects/pinn-cp-pfd-benchmark/results/cp_study/sens/h40/tau_c0_x0p8
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
  tau_c0  = 4.23318
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 47.9872
  k1      = 0.198833
  y_c     = 0.00343758
[]
