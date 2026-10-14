# PINNACLE deck TEMPLATE: crystal plasticity + phase-field ductile fracture on
# the perforated polycrystal plate (hexaform mesh, Part I orientation point
# file). Rendered by scripts/pfd_plate.py from a pipeline deck; every @KEY@ is
# substituted, nothing is edited by hand.
#
# UNITS: micrometer, megapascal, second (the CP study's system).
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
  finite_strain = true
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
[Functions]
  [./pull]
    # End-face displacement: strain x width x t / t_end.
    kind       = parsed
    expression = @STRAIN@*100.0*t/@TEND@
  []
[]
[Materials]
  [./cp]
    type                  = dislocation_density_cp
    C11                   = @C11@
    C12                   = @C12@
    C44                   = @C44@
    mu                    = @MU@
    gamma0                = @GAMMA0@
    rate_exponent_n       = @NEXP@
    tau_c0                = @TAUC0@
    q_lat                 = @QLAT@
    alpha0                = @ALPHA0@
    burgers               = @BURGERS@
    rho0                  = @RHO0@
    k1                    = @K1@
    annihilation_diameter = @YC@
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'
    []
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
    path          = @POINTFILE@
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
  n_steps         = @NSTEPS@
  max_iterations  = @MAXALT@
  tolerance       = @ALTTOL@
  snapshots_every = @SNAPSHOTS@
  stop_damage     = @STOPDAMAGE@    # end the path once max d reaches this (0 = full path)
  stop_on_mechanics_failure = true
  [./mechanics]
    variable          = u
    type              = cp_quasi_static
    n_steps           = @SUBSTEPS@
    t_end             = @TEND@
    newton_tol        = @NEWTONTOL@
    newton_rel_tol    = @NEWTONRELTOL@
    max_iter          = 50
    max_cutbacks      = @MAXCUTBACKS@
    linear_method     = direct
    sparse_threshold  = 0
    line_search       = @LINESEARCH@
    reaction_node_set = @REACTIONSET@
    corrector         = @CORRECTOR@       # newton; with a predictor: auto (verified PINN-as-solver route)
    corrector_tol     = @CORRECTORTOL@
    predictor_reject_ratio = 10
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
@UTILITIES@
[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement damage stress strain equivalent_plastic_strain'
[]
