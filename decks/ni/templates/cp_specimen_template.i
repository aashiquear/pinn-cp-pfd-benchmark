# PINNACLE deck TEMPLATE: crystal plasticity on a single-crystal Ni dog-bone
# specimen (hexaform `dogbone` mesh, Part III). Rendered by
# scripts/ni_specimen.py from a specimen deck; every @KEY@ is substituted,
# nothing is edited by hand.
#
# UNITS: millimetre, megapascal, second (the paper's system, docs/11 §4).
# Loading: the x- grip face clamped (u = 0), the x+ grip face clamped in y
# and z and pulled in x by the ramp `pull` to @STRAIN@ x the bar length at
# t_end; the gauge strain between the first and last thickness station is
# recorded per step (gauge_node_sets), the load as the reaction on x+.
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
[]
[Functions]
  [./pull]
    # End-face displacement: grip strain x bar length x t / t_end.
    kind       = parsed
    expression = @STRAIN@*@LENGTH@*t/@TEND@
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
      euler_angles = '@EULER@'
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
[BCs]
  [./clamp_x]
    variable   = u
    kind       = dirichlet
    node_set   = x-
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
  [./grip_yz]
    variable   = u
    kind       = dirichlet
    node_set   = x+
    components = [1, 2]
    values     = [0.0, 0.0]
  []
  [./pull_x]
    variable   = u
    kind       = dirichlet
    node_set   = x+
    components = [0]
    values     = [1.0]
    function   = pull
  []
[]
[Executioner]
  type              = cp_quasi_static
  n_steps           = @NSTEPS@
  t_end             = @TEND@
  newton_tol        = @NEWTONTOL@
  newton_rel_tol    = @NEWTONRELTOL@
  max_iter          = 50
  max_cutbacks      = @MAXCUTBACKS@
  linear_method     = @LINEAR@
  sparse_threshold  = @SPARSETHRESHOLD@
  line_search       = backtracking
  snapshots_every   = @SNAPSHOTS@
  reaction_node_set = x+
  gauge_node_sets   = 'gage- gage+'
  corrector         = @CORRECTOR@
  corrector_tol     = @CORRECTORTOL@
  predictor_reject_ratio = @REJECTRATIO@
[]
@UTILITIES@
[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement stress strain equivalent_plastic_strain slip_activity mean_rho'
[]
