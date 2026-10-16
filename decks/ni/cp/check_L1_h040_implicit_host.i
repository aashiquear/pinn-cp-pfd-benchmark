# CP-only check of the Ni dog-bone pipeline (docs/11 §5 stage 2): specimen L.1
# on the coarse h040 mesh, the paper's parameter set (Table 3), to a
# grip strain of 5% in 50 steps (1e-3 grip increment), the implicit update (ADR 0033);
# compared with the true-uniaxial material point (docs/12 §1.1).
# Local run, numpy, 4 threads. Not a result.
[Specimen]
  out_dir         = ../../../results/ni_study/cp/check_L1_h040_implicit_host
  mesh            = ../../../results/ni_study/mesh/h040/L.1/mesh.e
  orientations    = ../../../data/extracted/ni_table1_orientations.csv
  group           = L
  thickness       = ../../../data/extracted/ni_table2_thickness.csv
  specimen        = L.1
  gage_width      = 1.85
  length          = 25.4
  strain          = 0.05
  t_end           = 12.7          # grip displacement 0.1 mm/s (author, 2026-10-15): strain x 25.4 mm / 0.1 mm/s
  n_steps         = 50
  snapshots_every = 10
  backend         = torch
  use_gpu         = false
  threads         = 2
  newton_tol      = 1.0e-3
  newton_rel_tol  = 0
  max_cutbacks    = 10
  linear_method   = direct
  sparse_threshold = 500
  C11     = 253000
  C12     = 152000
  C44     = 124000
  mu      = 26000
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 0.26
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 2.5e-7
  rho0    = 2.35e6
  k1      = 0.0295
  y_c     = 0.79e-6
  predictor = none
  state_update    = implicit
[]
