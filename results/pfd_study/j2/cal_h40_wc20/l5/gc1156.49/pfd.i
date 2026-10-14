# PINNACLE deck TEMPLATE: J2 plasticity + phase-field ductile fracture on the
# perforated plate (hexaform mesh). Rendered by scripts/pfd_plate.py from a
# pipeline deck; every @KEY@ is substituted, nothing is edited by hand.
#
# UNITS: micrometer, megapascal, second (the CP study's system); Gc in
# MPa um = J/m^2 x 1e-6 ... see docs/workflow_pfd_plate.md for the conversion.
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
  path   = ../../../../mesh/h40_hole5/mesh.e
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
[Materials]
  [./j2]
    type           = j2_isotropic_hardening
    E              = 71160.8776844
    nu             = 0.345704948646
    sigma_y0       = 348.973838801
    K              = 0.0
    n              = 1.0
    hardening_mode = piecewise
    plastic_strain_points = [0, 0.001103404113, 0.001457760442, 0.001890566117, 0.002236465332, 0.002829409103, 0.00380594512, 0.004578460728, 0.006311751205, 0.007282276745, 0.009562935224, 0.01258409712, 0.01521540543, 0.02045979418, 0.02571600259, 0.03268534209, 0.03996487847, 0.05163012401, 0.06409474713, 0.08333478431, 0.1051847222, 0.133863871, 0.1487245395]
    yield_stress_points   = [348.9738388, 348.9738388, 350.2251486, 350.2251486, 367.5473601, 371.6000299, 381.3202262, 385.6983504, 394.0693825, 397.3442748, 405.7411151, 418.061362, 421.9701515, 438.7878043, 452.6626513, 470.6982324, 488.3785644, 517.1306282, 546.9541556, 589.0047449, 632.0621951, 680.9188223, 694.0109532]
    finite_strain  = true
    split_smoothing = 0   # MPa: rounds the tension-compression split (Newton robustness)
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = 1156.49
    length_scale = 5
    damage_model = at2
    anisotropy   = isotropic
    viscosity    = 0
  []
  [./drive]
    kind                   = pfd_driving_force
    initial                = 0.0
    plastic_work_prefactor = 1
    plastic_work_threshold = 20
  []
  [./degradation]
    kind = pfd_degradation
    form = quadratic
  []
[]
[Kernels]
  [./sd_j2]
    type        = stress_divergence_j2
    variable    = u
    dim         = 3
    material    = j2
    degradation = degradation
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
    values     = [20]      # strain x width, ramped by the load factor
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
    type              = fem_quasi_static
    n_steps           = 1
    ramp_dirichlet    = true
    max_iter          = 50
    tol               = 0.001
    rel_tol           = 0
    linear_method     = direct
    sparse_threshold  = 0          # sparse LU always: the plate has thousands of DOFs
    line_search       = backtracking   # backtracking keeps Newton inside the basin on degraded states
    max_cutbacks      = 6
    reaction_node_set = x+
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
  fields          = 'displacement damage stress strain'
[]
