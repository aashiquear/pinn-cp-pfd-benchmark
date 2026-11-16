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
  path   = ../../../mesh/h020/L.3/mesh.e
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
    expression = 0.3*25.4*t/6
  []
[]
[Materials]
  [./cp]
    type                  = dislocation_density_cp
    C11                   = 253000
    C12                   = 152000
    C44                   = 124000
    mu                    = 26000
    gamma0                = 1e-05
    rate_exponent_n       = 8.33333333333
    tau_c0                = 0.24397
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 2.5e-07
    rho0                  = 725576
    k1                    = 0.0114677
    annihilation_diameter = 1.54209e-06
    [./crystal]
      structure    = fcc
      euler_angles = '60 23.5 80.2'
    []
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = 0.1
    length_scale = 0.6
    damage_model = at2
    anisotropy   = isotropic
    viscosity    = 0
  []
  [./drive]
    kind                   = pfd_driving_force
    initial                = 0.0
    plastic_work_prefactor = 1
    plastic_work_threshold = 75
    elastic_prefactor      = 0   # chi on psi+ (1 = the form before the key; 0 = purely ductile driving force)
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
    state_update     = implicit   # implicit = backward-Euler slip update (PINNACLE ADR 0033)
    state_iterations = 2
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
  n_steps         = 300
  max_iterations  = 3
  tolerance       = 0.001
  snapshots_every = 2
  stop_damage     = 0    # end the path once max d reaches this (0 = full path)
  stop_separation = 'x- x+'         # end the path at through-section failure: d > 0.85 cuts the grips apart (PINNACLE 14fa4c6)
  stop_separation_damage = 0.85
  stop_on_mechanics_failure = true
  checkpoint_every = 8   # restartable state every N load steps (0 = none); a run that dies keeps its curve
  checkpoint_path  = solver/restart
  [./mechanics]
    variable          = u
    type              = cp_quasi_static
    n_steps           = 1
    t_end             = 6
    newton_tol        = 0.001
    newton_rel_tol    = 0
    max_iter          = 50
    max_cutbacks      = 10
    linear_method     = direct
    sparse_threshold  = 500
    line_search       = backtracking
    reaction_node_set = x+
    gauge_node_sets   = 'gage- gage+'
    corrector         = newton
    corrector_tol     = 0.05
    predictor_reject_ratio = 10
  []
  [./damage]
    variable          = d
    type              = pfd_steady
    n_load_steps      = 1
    boundary_node_set = grips
    boundary_damage   = 0.0
    sparse_threshold  = 0
    linear_method     = direct
    damage_threshold  = 0.8
    dt                = 0.02
  []
[]

[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement damage stress strain equivalent_plastic_strain'
[]
