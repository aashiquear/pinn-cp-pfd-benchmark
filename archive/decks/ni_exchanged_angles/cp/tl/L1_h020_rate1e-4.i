# Ni single crystal, specimen L.1, h020: the paper's Table 3 set UNCHANGED, total-Lagrangian
# kinematics (PINNACLE ADR 0035), strain rate 1e-4/s on the gauge. Restart of Part III
# (2026-11-15): the strain rate is the one quantity the paper does not give and is set here
# against Fig. 3. CP only (no damage); grip strain 30 % in 300 steps.
[Specimen]
  out_dir = ../../../../results/ni_study/tl/cp_L1_h020_rate1e-4
  mesh             = ../../../../results/ni_study/mesh/h020/L.1/mesh.e
  orientations     = ../../../../data/extracted/ni_table1_orientations.csv
  group            = L
  thickness        = ../../../../data/extracted/ni_table2_thickness.csv
  specimen         = L.1
  gage_width       = 1.85
  length           = 25.4
  strain           = 0.30
  t_end = 6000          # gauge strain rate about 9.8e-05/s (gauge strain 0.59 at 30 % bar strain)
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
  rho0 = 2.35e6
  k1 = 0.0295
  y_c = 0.79e-6
  predictor = none
[]
