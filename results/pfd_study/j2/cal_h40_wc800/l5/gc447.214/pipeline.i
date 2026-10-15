# J2+PFD base at plastic-work threshold W_c = 800 MPa (author 2026-10-14: push the
# softening past the hardening range, refit Gc lower); derived from base_h40.i.
# J2+PFD BASE case on the h40 perforated plate (Part II, docs/09 §5.1 step 2):
# the calibrated J2 piecewise law (results/pfd_study/j2_fit), AT2 damage with
# the plastic-work driving force (zeta = 1, W_c = 20 MPa chosen by the author on
# 2026-10-14 after the W_c = 0 sensitivity, eta = 0: quasi-static), 256 staggered steps to 20% at 1e-3 /s as in
# Part I, up to 3 alternations per step. Gc and l0 are the calibration
# targets; the values here are the bisection's starting point and are
# overwritten by scripts/pfd_calibrate.py and scripts/pfd_study.py cases.
# Production arm: torch CPU, 2 threads, 7 cases side by side (docs/09 §2).
[Plate]
  out_dir          = .
  mechanics        = j2
  mesh             = ../../../../mesh/h40_hole5/mesh.e
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
  reaction_set     = x+
  damage_set       = grips
  newton_tol       = 1.0e-3
  newton_rel_tol   = 0
  max_iter         = 50
  max_cutbacks     = 6
  line_search      = backtracking
  sparse_threshold = 500
  damage_threshold = 0.8
  j2_fit           = ../../../../j2_fit/j2_fit.json
  j2_law           = piecewise
  finite_strain    = true
  Gc               = 447.214
  l0               = 5
  eta              = 0.0
  zeta             = 1.0
  W_c              = 800.0
  damage_model     = at2
[]
