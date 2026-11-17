# Ni single crystal, specimen L.1, CP+PFD BASE deck of the G_c / l0 calibration (docs/12 §6):
# the calibrated h020 CP set with the Al2219 CP+PFD phase field (AT2, zeta 1, W_c 0, eta 0,
# three staggered alternations, tolerance 1e-3). Gc and l0 are set per run by
# scripts/ni_pfd_calibrate.py (decks/ni/pfd/cal_L1_h020.i). Path: grip strain 30% in 300
# steps (bar rate 5e-2/s); the run stops at through-section failure (PINNACLE
# stop_separation: d > 0.85 cuts the grips apart); snapshots every 2 steps, checkpoints every 8.
[Specimen]
  out_dir          = ../../../results/ni_study/pfd/base_L1_h020
  mesh             = ../../../results/ni_study/mesh/h020/L.1/mesh.e
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
  Gc               = 125
  l0               = 0.4
  zeta             = 1          # plastic-work driving force prefactor
  W_c              = 0          # plastic-work threshold (MPa): the W_c = 0 condition; the W_c stage sets it per deck
  eta              = 0
  chi              = 1          # elastic prefactor on psi+ (PINNACLE elastic_prefactor); 0 drops the elastic energy
  damage_model     = at2
  stop_damage      = 0
  checkpoint_every = 8
[]
