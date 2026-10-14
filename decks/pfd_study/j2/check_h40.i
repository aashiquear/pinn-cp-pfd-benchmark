# J2+PFD pipeline CHECK on the h40 perforated plate: 8 staggered steps to 1%
# strain, local machine, to validate the template, the reaction output and the
# damage bookkeeping before any campaign (docs/09 §5.1 step 2 prerequisite).
[Plate]
  out_dir         = ../../../results/pfd_study/j2/check_h40
  mechanics       = j2
  mesh            = ../../../results/pfd_study/mesh/h40_hole5/mesh.e
  width           = 100
  height          = 50
  thickness       = 1
  strain          = 0.01
  t_end           = 10.0
  n_steps         = 8
  substeps        = 1
  max_iterations  = 3
  tolerance       = 1.0e-4
  snapshots_every = 4
  backend         = numpy
  use_gpu         = false
  threads         = 4
  reaction_set    = x+
  damage_set      = grips
  newton_tol      = 1.0e-3
  newton_rel_tol  = 0
  max_iter        = 50
  max_cutbacks    = 6
  line_search     = backtracking
  sparse_threshold = 500
  damage_threshold = 0.8
  j2_fit          = ../../../results/pfd_study/j2_fit/j2_fit.json
  j2_law          = piecewise
  finite_strain   = true
  Gc   = 20.0
  l0   = 5.0
  eta  = 0.0
  zeta = 1.0
  W_c  = 0.0
  damage_model = at2
[]
