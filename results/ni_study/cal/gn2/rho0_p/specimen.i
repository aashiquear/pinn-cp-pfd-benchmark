# Ni single crystal, specimen L.1: calibration CENTRE (iteration 2, from gn1's step).
# calibration geometry (docs/11 §2, §5 stage 2): CP only, PINNACLE implicit state
# update (ADR 0033), grip strain 36% in 360 steps (1e-3 per step; the gauge strain
# runs about 1.5x the grip strain, to beyond the 0.51 of Fig. 3), bar strain rate 5e-2/s (grip speed 1.27 mm/s, about 7.5e-2/s on the gage).
# Parameters: Table 3 of the paper (data/extracted/ni_table3_cp_parameters.csv);
# rho0, k1 (k_m), y_c and tau_c0 are the calibration targets.
# Production arm: torch CPU, 2 threads, 7 side by side on the host (docs/11 §5).
[Specimen]
  out_dir          = .
  mesh             = ../../../mesh/h040/L.1/mesh.e
  orientations     = ../../../../../data/extracted/ni_table1_orientations.csv
  group            = L
  thickness        = ../../../../../data/extracted/ni_table2_thickness.csv
  specimen         = L.1
  gage_width       = 1.85
  length           = 25.4
  strain           = 0.36
  t_end            = 7.2          # bar strain rate 5e-2/s (author, 2026-10-15): strain / 0.05
  n_steps          = 360
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
  tau_c0  = 0.254803
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 2.5e-7
  rho0    = 1400370.0
  k1      = 0.0146493
  y_c     = 7.38183e-07
  predictor = none
[]
