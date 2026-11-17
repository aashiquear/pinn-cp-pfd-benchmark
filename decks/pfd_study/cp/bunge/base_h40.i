# CP+PFD BASE case on the h40 perforated polycrystal plate (Part II, docs/09
# §5.1 step 5): the Part I calibrated crystal-plasticity set and 30-grain
# point file, AT2 damage with the plastic-work driving force (zeta = 1, no
# threshold, eta = 0), 256 staggered steps to 20% at 1e-3 /s, up to 3
# alternations per step. Gc and l0 are the calibration targets (bisection by
# scripts/pfd_calibrate.py); the values here are the starting point.
# Production arm: torch CPU, 2 threads.
# Recheck of 2026-11-16: recalibrated Part I set (decks/cp_study/bunge/cal_center2.i), Bunge
# angles, total-Lagrangian balance (PINNACLE ADR 0035), implicit state update (ADR 0033).
# Gc is the recalibrated value at l0 = 5 (cal_h40.i in this folder: fracture at 17.11 %; the
# bisection itself started from the exchanged-angle value 1106.12).
[Plate]
  out_dir          = ../../../../results/pfd_study/cp/bunge/base_h40
  mechanics        = cp
  points           = ../../../../results/cp_study/bunge/cal/gn2/base/orientations_points.csv
  mesh             = ../../../../results/pfd_study/mesh/h40_hole5/mesh.e
  width            = 100
  height           = 50
  thickness        = 1
  strain           = 0.20
  t_end            = 200.0
  n_steps          = 256
  substeps         = 1
  max_iterations   = 3
  tolerance        = 1.0e-3
  snapshots_every  = 4
  stop_damage      = 0.95
  backend          = torch
  use_gpu          = false
  threads          = 2
  kinematics       = total_lagrangian
  state_update     = implicit
  reaction_set     = x+
  damage_set       = grips
  newton_tol       = 1.0e-3
  newton_rel_tol   = 0
  max_iter         = 50
  max_cutbacks     = 10
  line_search      = backtracking
  sparse_threshold = 500
  damage_threshold = 0.8
  C11              = 108.0e3
  C12              = 61.3e3
  C44              = 28.5e3
  mu               = 25.7e3
  gamma0           = 1.0e-5
  nexp             = 8.333333333333334
  tau_c0           = 5.38956
  q_lat            = 1.4
  alpha0           = 0.4
  burgers          = 0.000286
  rho0             = 53.9306
  k1               = 0.176771
  y_c              = 0.00271423
  Gc               = 1122.66        
  l0               = 5            
  eta              = 0.0
  zeta             = 1.0
  W_c              = 0.0
  damage_model     = at2
[]
