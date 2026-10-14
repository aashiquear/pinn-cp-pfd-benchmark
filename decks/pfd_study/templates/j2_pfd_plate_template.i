# PINNACLE deck TEMPLATE: J2 plasticity + phase-field ductile fracture on the
# perforated plate (hexaform mesh). Rendered by scripts/pfd_plate.py from a
# pipeline deck; every @KEY@ is substituted, nothing is edited by hand.
#
# UNITS: micrometer, megapascal, second (the CP study's system); Gc in
# MPa um = J/m^2 x 1e-6 ... see docs/workflow_pfd_plate.md for the conversion.
[Global]
  backend    = @BACKEND@
  use_gpu    = @USEGPU@
  precision  = float64
  seed       = 0
  output_dir = @OUTDIR@
[]
[Problem]
  type          = composed
  dimension     = 3
  finite_strain = @FINITESTRAIN@
[]
[Mesh]
  type   = file
  path   = @MESH@
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
    E              = @E@
    nu             = @NU@
    sigma_y0       = @SIGMAY0@
    K              = @KHARD@
    n              = @NHARD@
    hardening_mode = @HARDENING@
    plastic_strain_points = @ALPHAPOINTS@
    yield_stress_points   = @SIGMAPOINTS@
    finite_strain  = @FINITESTRAIN@
    split_smoothing = @SPLITSMOOTHING@   # MPa: rounds the tension-compression split (Newton robustness)
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = @GC@
    length_scale = @L0@
    damage_model = @DAMAGEMODEL@
    anisotropy   = isotropic
    viscosity    = @ETA@
  []
  [./drive]
    kind                   = pfd_driving_force
    initial                = 0.0
    plastic_work_prefactor = @ZETA@
    plastic_work_threshold = @WC@
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
    values     = [@PULL@]      # strain x width, ramped by the load factor
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
  n_steps         = @NSTEPS@
  max_iterations  = @MAXALT@
  tolerance       = @ALTTOL@
  snapshots_every = @SNAPSHOTS@
  stop_damage     = @STOPDAMAGE@    # end the path once max d reaches this (0 = full path)
  stop_on_mechanics_failure = true
  [./mechanics]
    variable          = u
    type              = fem_quasi_static
    n_steps           = @SUBSTEPS@
    ramp_dirichlet    = true
    max_iter          = @MAXITER@
    tol               = @NEWTONTOL@
    rel_tol           = @NEWTONRELTOL@
    linear_method     = direct
    sparse_threshold  = 0          # sparse LU always: the plate has thousands of DOFs
    line_search       = @LINESEARCH@   # backtracking keeps Newton inside the basin on degraded states
    max_cutbacks      = @MAXCUTBACKS@
    reaction_node_set = @REACTIONSET@
  []
  [./damage]
    variable          = d
    type              = pfd_steady
    n_load_steps      = 1
    boundary_node_set = @DAMAGESET@
    boundary_damage   = 0.0
    sparse_threshold  = @SPARSETHRESHOLD@
    linear_method     = direct
    damage_threshold  = @DTHRESHOLD@
    dt                = @DT@
  []
[]
[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement damage stress strain'
[]
