# Ni single crystal, specimen L.1: CALIBRATED CP set on the production mesh h020 (Gauss-Newton
# iteration 2, RMSE 4.1 MPa against Fig. 3; docs/12 §4). The run is the centre of
# decks/ni/cp/h020/cal_gn2.i (results/ni_study/cal_h020/gn2/base); every later Ni stage starts here.
# The h040 calibration (decks/ni/cp/calibrated_L1_h040.i) was 14-20% stiff against
# h020 (mesh study, docs/12 §3); the author chose to recalibrate on h020 (2026-10-17).
# Starts from the h040-calibrated set; grip strain 30% in 300 steps (1e-3 per step,
# 5e-2/s on the bar; the curves end at the first cutback, beyond 50% gauge strain).
# PINNACLE implicit state update (ADR 0033); torch CPU, 2 threads, memory cap 16 GB.
[Specimen]
  out_dir          = ../../../../results/ni_study/cal_h020/gn2/base
  mesh             = ../../../../results/ni_study/mesh/h020/L.1/mesh.e
  orientations     = ../../../../data/extracted/ni_table1_orientations.csv
  group            = L
  thickness        = ../../../../data/extracted/ni_table2_thickness.csv
  specimen         = L.1
  gage_width       = 1.85
  length           = 25.4
  strain           = 0.30
  t_end            = 6.0          # bar strain rate 5e-2/s: strain / 0.05
  n_steps          = 300
  snapshots_every  = 20
  backend          = torch
  use_gpu          = false
  threads          = 2
  newton_tol       = 1.0e-3
  newton_rel_tol   = 0
  max_cutbacks     = 10
  linear_method    = direct
  sparse_threshold = 500
  state_update     = implicit
  state_iterations = 2
  C11     = 253000
  C12     = 152000
  C44     = 124000
  mu      = 26000
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 0.24397
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 2.5e-7
  rho0    = 725576
  k1      = 0.0114677
  y_c     = 1.54209e-06
  predictor = none
[]
