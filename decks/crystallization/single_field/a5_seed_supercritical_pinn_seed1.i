# Part IV, case A5 (PINN route): a crystal seed 1.1 times the critical size in undercooled Ni,
# Delta T / T_f = 0.20, spherical symmetry.
#
# PHYSICS. Allen-Cahn dynamics of the single-field model on a radial line with the spherical
# volume element:
#
#     dphi/dt = M [ eps^2 (phi_rr + 2/r phi_r) - w g'(phi) - p'(phi) Delta_f ].
#
# The critical nucleus of case A3 is a saddle point of the free energy: a seed with the nucleus
# profile stretched to 1.1 of its radius lies on the growing side of it.
# The seed is the tanh droplet fitted to the nucleus (amplitude 0.9127, radius 0.9825 nm, width
# 0.21791 nm; solid volume 5.39 nm^3 against 5.29 of the nucleus) with the radius multiplied by 1.1.
#
# HOW IT IS SOLVED (PINN route, `pf_pinn`). One network phi_theta(r, t) (4 x 48, tanh)
# is trained on the strong residual of the equation above at 2000 collocation points over
# [0, 6] nm x [0, 4], in 4 time slabs of one time unit (1000 Adam epochs and 30
# L-BFGS iterations each). The initial condition and the range 0 <= phi <= 1 hold by
# construction of the ansatz; the boundary values are penalised.
# On the spherical mesh the network reads r^2 (the field is even in r) and the strong form
# carries the term (2/r) phi_r.
#
# EXPECTED RESULT. Compared with a5_seed_supercritical_window_fem.i (same window, FE): values in the report. The
# solid amount at the slab ends is in summary.json (`invariants`, `mass`).
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa; mobility M = 1, so time is
# in units of 1/(M GPa).
#
# PARAMETERS at T = 0.8 T_f: eps^2 = 0.329871 aJ/nm, w = 12.740387 GPa, Delta_f = -0.4923638 GPa.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 1
  output_dir = results/crystallization/single_field/a5_seed_supercritical_pinn_seed1
[]

[Problem]
  type          = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type              = line_mesh
  domain            = '0.0 6.0'
  divisions         = 600             # h = 0.01 nm
  coordinate_system = spherical
[]

[Variables]
  [./phi]
    family     = pinn
    components = 1
    network    = mlp
  []
[]

[Functions]
  [./seed]
    kind       = parsed
    expression = 0.9127*0.5*(1-tanh((x-1.0807)/(2*0.21791)))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = seed
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'
    [./coefficients]
      w  = 12.740387
      df = -0.4923638
    []
  []
[]

[Kernels]
  [./dphi_dt]
    type        = time_derivative
    variable    = phi
    coefficient = 1.0
  []
  [./bulk]
    type                 = ac_bulk
    variable             = phi
    material             = f_bulk
    mobility             = 1.0
    gradient_coefficient = 0.329871
    spatial_dim          = 1
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.329871    # eps^2(T)
    spatial_dim          = 1
  []
[]

[BCs]
  [./liquid]
    variable   = phi
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [0.0]
  []
[]

[Executioner]
  type                 = pf_pinn
  residual_form        = strong
  n_steps              = 200
  dt                   = 0.02        # horizon t = 4; the network is continuous in time
  slabs                = 4           # time-slab curriculum, one shared network
  continuity           = data
  slab_epochs          = 1000        # Adam epochs per slab
  adam_lr              = 2.0e-3
  lbfgs_iters          = 30          # then L-BFGS per slab
  causal               = false
  weighting            = fixed
  hidden_layers        = 4
  neurons              = 48
  transform_bounded    = true        # phi in [0, 1] by construction
  bound_lower          = 0.0
  bound_upper          = 1.0
  transform_ic_exact   = true        # the initial condition holds exactly
  collocation_interior = 2000
  n_grid               = 201
  seed                 = 1
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
