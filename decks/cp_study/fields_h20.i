# FIELD SERIES rerun of the h20 mesh with the calibrated set (same deck as
# mesh_h20.i except out_dir and snapshots_every = 4: one field snapshot every 4
# load steps = every 0.3125% strain, 64 frames to 20%), for the equivalent
# plastic strain animation of the mesh study (author's request 2026-10-11).
# CP calibration, FE Gauss-Newton FINAL h40 set on mesh h20 (mesh sensitivity): the Taylor-surrogate
# optimum of results/cp_study/taylor_screen/fit.json (RMSE 13.2 MPa at the
# material point). The FE polycrystal is softer than the isostrain aggregate,
# so this is the starting point of the FE correction, not the answer.
#
# Plate 100 x 50 x 1 um loaded along x (decided 2026-10-09), mesh level h40
# = 40 x 20 x 1 HEX8 (2.5 um elements), the 30-grain structure grown for 2.5
# interface times (results/graingrowth_pilot/n32_bdf1_h40), orientation seed
# 20261009. Loading to 20% nominal strain at 1e-3 /s (t_end = 200 s) in 256
# steps (7.8e-4 strain per step, under the G20 limit of 8e-4).
#
# STRAIN RATE is an assumption: Vaughan (2024) does not report it in the
# extracted record; 1e-3 /s is the quasi-static value the earlier Taylor
# calibration used. With n = 8.33 a factor 3 in rate moves the stress 14%.
#
# Parameters: Rezwan 2025 Tables 1-2 with the paper's y_c = 0.00268 (not the
# MOOSE test deck's 0.00286). Fixed by decision: C11 C12 C44 b mu q alpha0
# (and gamma0, n unless a convergence issue forces it). Free, in study order:
# rho0, k1 (k_m), y_c, tau_c0.

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
  out_dir      = ../../results/cp_study/fields/h20
  load_axis    = x
  strain       = 0.20
  t_end        = 200.0
  n_steps      = 256
  nx           = 20
  ny           = 10
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
