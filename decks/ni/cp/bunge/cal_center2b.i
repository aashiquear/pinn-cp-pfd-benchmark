# Ni L.1 h020, total-Lagrangian, Bunge convention (ADR 0036): candidate center of Gauss-Newton
# iteration 2. Damped step from the paper's set on the iteration-1 Jacobian (rho0, k_m, y_c in log
# space, Levenberg-Marquardt 0.001 of the mean diagonal, tau_c0 held): rho0 x1.154, k_m x0.711, y_c x0.419.
# Rate 1e-3/s.
[Specimen]
  out_dir = ../../../../results/ni_study/bunge/cal/gn2/center_b
  mesh             = ../../../../results/ni_study/mesh/h020/L.1/mesh.e
  orientations     = ../../../../data/extracted/ni_table1_orientations.csv
  group            = L
  thickness        = ../../../../data/extracted/ni_table2_thickness.csv
  specimen         = L.1
  gage_width       = 1.85
  length           = 25.4
  strain           = 0.30
  t_end = 600          # gauge strain rate about 1e-3/s (assumed quasi-static rate; the paper gives none)
  n_steps          = 300
  snapshots_every  = 20
  backend          = torch
  use_gpu          = false
  threads = 4
  newton_tol = 0.05
  newton_rel_tol   = 0
  max_cutbacks     = 10
  linear_method    = direct
  sparse_threshold = 500
  kinematics       = total_lagrangian
  state_update     = implicit
  state_iterations = 2
  C11     = 253000
  C12     = 152000
  C44     = 124000
  mu      = 26000
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0 = 0.26
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 2.5e-7
  rho0 = 2.7119e+06
  k1 = 0.020974
  y_c = 3.3101e-07
  predictor = none
[]
