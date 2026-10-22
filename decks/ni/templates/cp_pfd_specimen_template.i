# PINNACLE deck TEMPLATE: crystal plasticity + phase-field ductile fracture on a
# single-crystal Ni dog-bone specimen (hexaform `dogbone` mesh, Part III). Rendered by
# scripts/ni_specimen.py from a specimen deck that sets Gc (and l0); every @KEY@ is
# substituted, nothing is edited by hand.
#
# UNITS: millimeter, megapascal, second (the paper's system, docs/11 §4).
# Mechanics as decks/ni/templates/cp_specimen_template.i (x- clamped, x+ clamped in y
# and z and pulled in x; gauge strain from gage-/gage+, reaction on x+); fracture as the
# Al2219 CP+PFD plate (decks/pfd_study/templates/cp_pfd_plate_template.i): AT2 phase
# field driven by the plastic work (zeta, threshold W_c), quadratic degradation,
# staggered mechanics/damage alternations per load step, no damage on the grip faces.
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
    state_update     = @STATEUPDATE@   # implicit = backward-Euler slip update (PINNACLE ADR 0033)
    state_iterations = @STATEITER@
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
  type            = staggered
  n_steps         = @NSTEPS@
  max_iterations  = @MAXALT@
  tolerance       = @ALTTOL@
  snapshots_every = @SNAPSHOTS@
  stop_damage     = @STOPDAMAGE@    # end the path once max d reaches this (0 = full path)
  stop_separation = 'x- x+'         # end the path at through-section failure: d > @SEPDAMAGE@ cuts the grips apart (PINNACLE 14fa4c6)
  stop_separation_damage = @SEPDAMAGE@
  stop_on_mechanics_failure = true
  checkpoint_every = @CHECKPOINT@   # restartable state every N load steps (0 = none); a run that dies keeps its curve
  checkpoint_path  = solver/restart
  [./mechanics]
    variable          = u
    type              = cp_quasi_static
    n_steps           = @SUBSTEPS@
    t_end             = @TEND@
    newton_tol        = @NEWTONTOL@
    newton_rel_tol    = @NEWTONRELTOL@
    max_iter          = 50
    max_cutbacks      = @MAXCUTBACKS@
    linear_method     = @LINEAR@
    sparse_threshold  = @SPARSETHRESHOLD@
    line_search       = backtracking
    reaction_node_set = x+
    gauge_node_sets   = 'gage- gage+'
    corrector         = @CORRECTOR@
    corrector_tol     = @CORRECTORTOL@
    predictor_reject_ratio = @REJECTRATIO@
  []
  [./damage]
    variable          = d
    type              = pfd_steady
    n_load_steps      = 1
    boundary_node_set = grips
    boundary_damage   = 0.0
    sparse_threshold  = 0
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
