# BASE specimen deck of the last stage, route 'fem' (docs/12 §8; author 2026-10-25): the carried
# CP+PFD set (chi 1, W_c 0, G_c 30.95 MPa mm, l0 0.6 mm: the W_c = 0 calibration of docs/12 §6.2) with the calibrated h020 CP set,
# predictor = none. The 15 specimens are the [Cases] of fem_15_h020.i (scripts/ni_study.py).
[Specimen]
  out_dir          = ../../../results/ni_study/final/fem_wc0/L.1
  mesh             = ../../../results/ni_study/mesh/h020/L.1/mesh.e
  orientations     = ../../../data/extracted/ni_table1_orientations.csv
  group            = L
  thickness        = ../../../data/extracted/ni_table2_thickness.csv
  specimen         = L.1
  gage_width       = 1.85
  length           = 25.4
  strain           = 0.40          # 2026-10-26: the 30% path ends near 60% gauge strain, below the 61.6 / 62.6% failures of 9.3 / 9.4 and the soft synthetic specimens
  t_end            = 8.0          # bar strain rate 5e-2/s: strain / 0.05
  n_steps          = 400          # same 1e-3 grip strain per step
  snapshots_every  = 2
  backend          = torch
  use_gpu          = false
  threads          = 2
  newton_tol       = 5.0e-2   # 2026-10-26: Newton floor of the implicit stress's FD tangent on multi-slip orientations (S01 1e-2 N); 5e-4 of the reaction; reactions unchanged at 5e-3 to 1e-7
  newton_rel_tol   = 1.0e-6
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
  memory_cap_gb    = 16
  Gc               = 30.95   # 2026-10-28: W_c = 0 pair of l0 = 0.6 (author: W_c = 75 calibrated on L transfers late on the other orientations)
  l0               = 0.6
  zeta             = 1          # plastic-work driving force prefactor
  W_c              = 0
  eta              = 0
  chi              = 1   # as in the W_c = 0 calibration (chi had its default 1 there)
  damage_model     = at2
  stop_damage      = 0
  checkpoint_every = 8
[]
