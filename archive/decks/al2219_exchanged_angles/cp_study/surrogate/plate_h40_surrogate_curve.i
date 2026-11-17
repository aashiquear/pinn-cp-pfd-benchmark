# PLATE SURROGATE, reduced target (ADR 0030): a network sigma_eng(rho0, k1, y_c,
# tau_c0, strain) of the engineering stress-strain curve of the calibrated h40
# plate over the four calibrated parameters (x0.5 .. x2, log), trained on the
# 24 FE samples of decks/cp_study/surrogate/design_h40.i (data source; crystal
# plasticity is history dependent, so the physics source is refused). Run from
# the repository root after `surrogate_design.py ... all`:
#   pinnacle run decks/cp_study/surrogate/plate_h40_surrogate_curve.i
# The network trains on the torch GPU; `evaluate` predicts 1,000 random
# parameter rows (scripts/surrogate_queries.py) and logs the time per query.
[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 0
  output_dir = results/cp_study/surrogate/h40_curve
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]

[Mesh]
  type      = box_mesh
  domain    = '0 0 0 100 50 1'
  divisions = '40 20 1'
[]

[Variables]
  [./u]
    family     = fem
    components = 3
  []
[]

[Functions]
  [./pull]
    # Top-face displacement: strain x height x t / t_end.
    kind       = parsed
    expression = 0.2*100.0*t/200
  []
[]

[Materials]
  [./cp]
    type                  = dislocation_density_cp
    C11                   = 108000
    C12                   = 61300
    C44                   = 28500
    mu                    = 25700
    gamma0                = 1e-05
    rate_exponent_n       = 8.33333333333
    tau_c0                = 5.29147
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 0.000286
    rho0                  = 47.9872
    k1                    = 0.198833
    annihilation_diameter = 0.00343758
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'
    []
  []
[]

[Kernels]
  [./sd_cp]
    type          = stress_divergence_cp
    variable      = u
    material      = cp
    dim           = 3
    finite_strain = true
    gnd           = false
  []
[]

[OrientationField]
  [./source]
    type          = point_file
    path          = results/cp_study/cal/gn2/base/orientations_points.csv
    interpolation = nearest
  []
[]

[BCs]
  [./x_left]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./y_bottom]
    variable   = u
    kind       = dirichlet
    node_set   = ymin
    components = [1]
    values     = [0.0]
  []
  [./z_back]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
  [./x_pull]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
    function   = pull
  []
[]

[Parameters]
  [./rho0]
    path  = materials/cp/rho0
    range = '23.9936 95.9744'
    scale = log
  []
  [./k1]
    path  = materials/cp/k1
    range = '0.0994165 0.397666'
    scale = log
  []
  [./y_c]
    path  = materials/cp/annihilation_diameter
    range = '0.00171879 0.00687516'
    scale = log
  []
  [./tau_c0]
    path  = materials/cp/tau_c0
    range = '2.64574 10.5829'
    scale = log
  []
[]
[Executioner]
  type          = parametric_surrogate
  sources       = data
  target        = reduced
  columns       = eng_stress
  abscissa      = nominal_strain
  data_dir      = results/cp_study/surrogate/h40_data/samples
  samples_file  = results/cp_study/surrogate/h40_data/design.csv
  holdout       = 0.25
  seed          = 0
  hidden_layers = 3
  neurons       = 64
  activation    = tanh
  epochs_adam   = 4000
  epochs_lbfgs  = 50
  lr            = 2.0e-3
  certify       = false
  evaluate      = results/cp_study/surrogate/queries_1000.csv
  [./generator]
    type             = cp_quasi_static
    n_steps          = 256
    t_end            = 200
    newton_tol       = 0.001
    # Relative test off (0): a warm start from a predictor must be held to the
    # SAME absolute criterion as a cold start; the relative-to-first-residual
    # test penalises a good initial guess (docs/07 §10.0 of the study).
    newton_rel_tol   = 0
    max_iter         = 50
    max_cutbacks     = 10
    # direct (sparse LU) by default: conjugate gradients stalls once the
    # finite-strain CP tangent loses positive definiteness (7.8-9.8% strain on
    # this plate); gmres is the iterative alternative for larger meshes.
    linear_method    = direct
    sparse_threshold = 0
    # Per-step field snapshots (0 = final state only); N writes
    # fields_step{NNNN}.vtk + fields.pvd every N accepted steps.
    snapshots_every  = 4
    # newton (default) corrects every step; none accepts a [Utilities]
    # predictor's proposal as the solution (PINN as solver, ADR 0028).
    corrector        = newton
  []
[]


[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement stress strain equivalent_plastic_strain slip_activity mean_rho'
[]
