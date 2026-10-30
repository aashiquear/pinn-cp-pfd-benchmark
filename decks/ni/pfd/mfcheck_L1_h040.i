# Multi-fidelity check (docs/13 §4.7 P7.4; author 2026-10-29): can the surrogate's training
# runs use h040 instead of h020? Specimen L.1 with the carried CP+PFD set (the only orientation
# with its own W_c: chi 0, W_c 75 MPa, G_c 0.1 MPa mm, l0 0.6 mm), identical to the h020
# reference run results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1 (failure 0.4639, peak
# 241.5 MPa) except for the mesh. Note l0/h = 1.5 on h040 (3 on h020). One run.
[Specimen]
  out_dir          = ../../../results/ni_study/pfd/mfcheck_L1_h040
  mesh             = ../../../results/ni_study/mesh/h040/L.1/mesh.e   # the ONLY change against the h020 reference
  orientations     = ../../../data/extracted/ni_table1_orientations.csv
  group            = L
  thickness        = ../../../data/extracted/ni_table2_thickness.csv
  specimen         = L.1
  gage_width       = 1.85
  length           = 25.4
  strain           = 0.30
  t_end            = 6.0          # bar strain rate 5e-2/s: strain / 0.05
  n_steps          = 300
  snapshots_every  = 2
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
  memory_cap_gb    = 16
  Gc               = 0.1
  l0               = 0.6
  zeta             = 1          # plastic-work driving force prefactor
  W_c              = 75         # plastic-work threshold (MPa), the L calibration
  eta              = 0
  chi              = 0          # purely ductile driving force
  damage_model     = at2
  stop_damage      = 0
  checkpoint_every = 8
[]
