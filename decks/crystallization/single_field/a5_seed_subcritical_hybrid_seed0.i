# Part IV, case A5 (PINN+FEM route): a crystal seed 0.9 times the critical size in undercooled Ni,
# Delta T / T_f = 0.20, spherical symmetry.
#
# PHYSICS. Allen-Cahn dynamics of the single-field model on a radial line with the spherical
# volume element:
#
#     dphi/dt = M [ eps^2 (phi_rr + 2/r phi_r) - w g'(phi) - p'(phi) Delta_f ].
#
# The critical nucleus of case A3 is a saddle point of the free energy: a seed with the nucleus
# profile stretched to 0.9 of its radius lies on the melting side of it.
# The seed is the tanh droplet fitted to the nucleus (amplitude 0.9127, radius 0.9825 nm, width
# 0.21791 nm; solid volume 5.39 nm^3 against 5.29 of the nucleus) with the radius multiplied by 0.9.
#
# HOW IT IS SOLVED (PINN+FEM route, `pf_pinn_fem_hybrid`). Predictor: the network of the
# PINN deck, trained over the horizon. Corrector: the FE march of a5_seed_subcritical_window_fem.i (BDF2,
# dt = 0.02); the Newton iteration of every step starts from the network's field at the end
# time of that step instead of from the previous state. Every state is the corrector's.
# solver_history.csv records per step `newton_iters`, `predictor_rms` (distance of the
# prediction from the converged state) and `predictor_mass` (the solid amount of the
# predicted field).
#
# EXPECTED RESULT. The fields of a5_seed_subcritical_window_fem.i to the Newton tolerance. What the network
# changes is the number of Newton iterations per step; values in the report.
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa; mobility M = 1, so time is
# in units of 1/(M GPa).
#
# PARAMETERS at T = 0.8 T_f: eps^2 = 0.329871 aJ/nm, w = 12.740387 GPa, Delta_f = -0.4923638 GPa.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a5_seed_subcritical_hybrid_seed0
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
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  [./seed]
    kind       = parsed
    expression = 0.9127*0.5*(1-tanh((x-0.8842)/(2*0.21791)))
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
  type = pf_pinn_fem_hybrid
  [./pinn]
    type                 = pf_pinn
    residual_form        = strong
    n_steps              = 150
    dt                   = 0.02        # horizon t = 3; the network is continuous in time
    slabs                = 3           # time-slab curriculum, one shared network
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
    seed                 = 0
  []
  [./fem]
    type            = pf_fem_transient
    scheme          = bdf2
    n_steps         = 150
    dt              = 0.02       # to t = 3
    newton_tol      = 1.0e-10
    snapshots_every = 5
  []
[]

[Postprocessors]
  [./free_energy]
    type          = total_free_energy
    bulk_field    = field_phi
    free_energy   = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.329871
      []
    []
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
