# Ni single crystal, specimen L.1, CP+PFD CHECK run (Part III stage 6, docs/12 §6): the
# calibrated h020 CP set (decks/ni/cp/h020/calibrated_L1_h020.i) with the phase field of the
# Al2219 CP+PFD plate (AT2, plastic-work driving force zeta = 1, W_c = 0, eta = 0, up to three
# staggered alternations per load step, tolerance 1e-3 on the damage increment). G_c = 125
# MPa mm is the paper's value (rezwan 2026 predicting, Sec. 3.5; its l0 is not stated);
# l0 = 0.4 mm = 2h on the h020 gage cell, the Al2219 production ratio. Purpose: the run cost,
# where and how L.1 fails at the paper's G_c, and the end-to-end check of the failure
# criterion (scripts/ni_specimen.py fracture: damage > 0.85 cutting the specimen in two)
# before the G_c / l0 calibration to the Fig. 4(b) failure strain (L.1: 0.463).
# Path: grip strain 30% in 300 load steps (6 s, bar rate 5e-2/s) as the CP runs; the run
# stops once a node reaches d = 0.995; snapshots every 2 steps, checkpoints every 8.
[Specimen]
  out_dir          = ../../../results/ni_study/pfd/check_L1_h020
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
  stop_damage      = 0.995
  checkpoint_every = 8
[]
