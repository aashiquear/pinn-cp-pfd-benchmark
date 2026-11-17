# CP+PFD pipeline CHECK on the h40 perforated polycrystal plate: 4 staggered
# steps to 0.5% strain on the local machine, to validate the CP template
# (cp_quasi_static inside staggered with the Part I point file, degradation,
# reaction output, damage bookkeeping) before any campaign. Brittle Gc on
# purpose so damage moves within the short path.
[Plate]
  out_dir          = ../../../results/pfd_study/cp/check_h40_linear
  mechanics        = cp
  mesh             = ../../../results/pfd_study/mesh/h40_hole5/mesh.e
  points           = ../../../results/cp_study/cal/gn2/base/orientations_points.csv
  width            = 100
  height           = 50
  thickness        = 1
  strain           = 0.005
  t_end            = 5.0
  n_steps          = 4
  substeps         = 2
  max_iterations   = 3
  tolerance        = 1.0e-3
  snapshots_every  = 2
  stop_damage      = 0.95
  backend          = numpy
  use_gpu          = false
  threads          = 3
  reaction_set     = x+
  damage_set       = grips
  newton_tol       = 1.0e-3
  newton_rel_tol   = 0
  max_iter         = 50
  max_cutbacks     = 6
  line_search      = backtracking
  sparse_threshold = 500
  damage_threshold = 0.8
  C11     = 108.0e3
  C12     = 61.3e3
  C44     = 28.5e3
  mu      = 25.7e3
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 5.29147
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 47.9872
  k1      = 0.198833
  y_c     = 0.00343758
  Gc           = 50.0
  predictor    = linear
  l0           = 5.0
  eta          = 0.0
  zeta         = 1.0
  W_c          = 0.0
  damage_model = at2
[]
