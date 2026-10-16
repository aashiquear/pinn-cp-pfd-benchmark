# PINNACLE deck TEMPLATE: crystal plasticity on a single-crystal Ni dog-bone
# specimen (hexaform `dogbone` mesh, Part III). Rendered by
# scripts/ni_specimen.py from a specimen deck; every @KEY@ is substituted,
# nothing is edited by hand.
#
# UNITS: millimetre, megapascal, second (the paper's system, docs/11 §4).
# Loading: the x- grip face clamped (u = 0), the x+ grip face clamped in y
# and z and pulled in x by the ramp `pull` to 0.36 x the bar length at
# t_end; the gauge strain between the first and last thickness station is
# recorded per step (gauge_node_sets), the load as the reaction on x+.
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
  path   = ../../../mesh/h040/L.1/mesh.e
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
    expression = 0.36*25.4*t/7.2
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
    tau_c0                = 0.249155
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 2.5e-07
    rho0                  = 763856
    k1                    = 0.0129672
    annihilation_diameter = 1.47916e-06
    [./crystal]
      structure    = fcc
      euler_angles = '60 23.5 80.2'
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
    state_update     = implicit   # explicit | iterated | implicit (PINNACLE ADR 0033; implicit = backward-Euler slip update)
    state_iterations = 2     # passes of the iterated update per load step
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
  n_steps           = 360
  t_end             = 7.2
  newton_tol        = 0.001
  newton_rel_tol    = 0
  max_iter          = 50
  max_cutbacks      = 10
  linear_method     = direct
  sparse_threshold  = 500
  line_search       = backtracking
  snapshots_every   = 20
  reaction_node_set = x+
  gauge_node_sets   = 'gage- gage+'
  corrector         = newton
  corrector_tol     = 0.05
  predictor_reject_ratio = 10
[]

[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement stress strain equivalent_plastic_strain slip_activity mean_rho'
[]
