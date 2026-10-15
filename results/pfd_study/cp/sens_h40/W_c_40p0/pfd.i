# PINNACLE deck TEMPLATE: crystal plasticity + phase-field ductile fracture on
# the perforated polycrystal plate (hexaform mesh, Part I orientation point
# file). Rendered by scripts/pfd_plate.py from a pipeline deck; every @KEY@ is
# substituted, nothing is edited by hand.
#
# UNITS: micrometer, megapascal, second (the CP study's system).
[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = solver
[]
[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]
[Mesh]
  type   = file
  path   = ../../../mesh/h40_hole5/mesh.e
  format = exodus
  scale  = 1.0
  verify = true
[]
[Variables]
  [./u]
    family     = fem
    components = 3
    order      = 2        # 2x2x2 Gauss points: order 1 (one point) hourglasses on hex8
  []
  [./d]
    family     = fem
    components = 1
    order      = 2
  []
[]
[Functions]
  [./pull]
    # End-face displacement: strain x width x t / t_end.
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
  [./frac]
    kind         = pfd_fracture
    Gc           = 1106.12
    length_scale = 5
    damage_model = at2
    anisotropy   = isotropic
    viscosity    = 0
  []
  [./drive]
    kind                   = pfd_driving_force
    initial                = 0.0
    plastic_work_prefactor = 1
    plastic_work_threshold = 40
  []
  [./degradation]
    kind = pfd_degradation
    form = quadratic
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
    degradation   = degradation
  []
  [./pfd]
    type            = pfd_evolution
    variable        = d
    dim             = 3
    fracture        = frac
    driving_force   = drive
    boundary_damage = 0.0
  []
[]
[OrientationField]
  [./source]
    type          = point_file
    path          = ../../../../cp_study/cal/gn2/base/orientations_points.csv
    interpolation = nearest
  []
[]
[BCs]
  [./x_left]
    variable   = u
    kind       = dirichlet
    node_set   = x-
    components = [0]
    values     = [0.0]
  []
  [./x_pull]
    variable   = u
    kind       = dirichlet
    node_set   = x+
    components = [0]
    values     = [1.0]
    function   = pull
  []
  [./y_bottom]
    variable   = u
    kind       = dirichlet
    node_set   = y-
    components = [1]
    values     = [0.0]
  []
  [./z_back]
    variable   = u
    kind       = dirichlet
    node_set   = z-
    components = [2]
    values     = [0.0]
  []
[]
[Executioner]
  type            = staggered
  n_steps         = 256
  max_iterations  = 3
  tolerance       = 0.001
  snapshots_every = 4
  stop_damage     = 0.95    # end the path once max d reaches this (0 = full path)
  stop_on_mechanics_failure = true
  [./mechanics]
    variable          = u
    type              = cp_quasi_static
    n_steps           = 1
    t_end             = 200
    newton_tol        = 0.001
    newton_rel_tol    = 0
    max_iter          = 50
    max_cutbacks      = 10
    linear_method     = direct
    sparse_threshold  = 0
    line_search       = backtracking
    reaction_node_set = x+
    corrector         = newton       # newton; with a predictor: auto (verified PINN-as-solver route)
    corrector_tol     = 0.05
    predictor_reject_ratio = 10
  []
  [./damage]
    variable          = d
    type              = pfd_steady
    n_load_steps      = 1
    boundary_node_set = grips
    boundary_damage   = 0.0
    sparse_threshold  = 500
    linear_method     = direct
    damage_threshold  = 0.8
    dt                = 0.78125
  []
[]

[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement damage stress strain equivalent_plastic_strain'
[]
