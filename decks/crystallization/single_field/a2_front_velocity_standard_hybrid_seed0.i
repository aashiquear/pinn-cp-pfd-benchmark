# Part IV, case A2 (PINN+FEM route): steady crystallization front in undercooled Ni,
# Delta T / T_f = 0.20, interpolation function p = phi^3 (10 - 15 phi + 6 phi^2).
#
# PHYSICS. Allen-Cahn dynamics of the single-field model (review Eq. 16, no noise):
#
#     dphi/dt = M [ eps^2 phi_xx - w g'(phi) - p'(phi) Delta_f ].
#
# A planar interface between crystal (phi = 1, left) and liquid (phi = 0, right) moves into the
# liquid at a constant velocity once the profile has adjusted:
#
#     no closed form; reference from the travelling-wave solution of scripts/cryst_reference.py
#     (decks/crystallization/reference/ni_front.i): v = 0.644490 nm per time unit
#
# The velocity is read from the solid length: `mass` of solver_history.csv is int phi dx, and
# v = d(mass)/dt (fit over t > 5, after the start-up transient).
#
# HOW IT IS SOLVED (PINN+FEM route, `pf_pinn_fem_hybrid`). Predictor: the network of the
# PINN deck, trained over the horizon. Corrector: the FE march of a2_front_velocity_standard_window_fem.i (BDF2,
# dt = 0.02); the Newton iteration of every step starts from the network's field at the end
# time of that step instead of from the previous state. Every state is the corrector's.
# solver_history.csv records per step `newton_iters`, `predictor_rms` (distance of the
# prediction from the converged state) and `predictor_mass` (the solid amount of the
# predicted field).
#
# EXPECTED RESULT. The fields of a2_front_velocity_standard_window_fem.i to the Newton tolerance. What the network
# changes is the number of Newton iterations per step; values in the report.
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa; mobility M = 1, so time is
# in units of 1/(M GPa).
#
# PARAMETERS at T = 0.8 T_f: eps^2 = 0.329871 aJ/nm, w = 12.740387 GPa, Delta_f = -0.4923638 GPa,
# delta = 0.227560 nm.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a2_front_velocity_standard_hybrid_seed0
[]

[Problem]
  type          = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type      = line_mesh
  domain    = '0.0 10.0'       # nm
  divisions = 1000             # h = 0.01 nm
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[ICs]
  [./phi0]
    variable  = phi
    kind      = tanh_interface   # crystal on the left; equilibrium width 4 delta = 0.91 nm
    axis      = 0
    position  = 2.5
    width     = 0.91024
    amplitude = 1.0
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

[Executioner]
  type = pf_pinn_fem_hybrid
  [./pinn]
    type                 = pf_pinn
    residual_form        = strong
    n_steps              = 300
    dt                   = 0.02        # horizon t = 6; the network is continuous in time
    slabs                = 6           # time-slab curriculum, one shared network
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
    n_steps         = 300
    dt              = 0.02       # to t = 6
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
