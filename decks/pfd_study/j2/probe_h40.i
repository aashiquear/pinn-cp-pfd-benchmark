# TIMING PROBE: the base deck truncated to its first 16 steps (same strain per step), host torch CPU 2 threads.
# J2+PFD BASE case on the h40 perforated plate (Part II, docs/09 §5.1 step 2):
# the calibrated J2 piecewise law (results/pfd_study/j2_fit), AT2 damage with
# the plastic-work driving force (zeta = 1, no threshold, eta = 0: the
# author's quasi-static choice), 256 staggered steps to 20% at 1e-3 /s as in
# Part I, up to 3 alternations per step. Gc and l0 are the calibration
# targets; the values here are the bisection's starting point and are
# overwritten by scripts/pfd_calibrate.py and scripts/pfd_study.py cases.
# Production arm: torch CPU, 2 threads, 7 cases side by side (docs/09 §2).
[Plate]
  out_dir          = ../../../results/pfd_study/j2/probe_h40
  mechanics        = j2
  mesh             = ../../../results/pfd_study/mesh/h40_hole5/mesh.e
  width            = 100
  height           = 50
  thickness        = 1
  strain           = 0.0125
  t_end            = 12.5
  n_steps          = 16
  substeps         = 1
  max_iterations   = 3
  tolerance        = 1.0e-4
  snapshots_every  = 0
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
  j2_fit           = ../../../results/pfd_study/j2_fit/j2_fit.json
  j2_law           = piecewise
  finite_strain    = true
  Gc               = 300.0
  l0               = 5.0
  eta              = 0.0
  zeta             = 1.0
  W_c              = 0.0
  damage_model     = at2
[]
