# Part IV, case A2 (FE route): steady crystallization front in undercooled Ni,
# Delta T / T_f = 0.20, interpolation function p = phi^2 (3 - 2 phi).
#
# PHYSICS. Allen-Cahn dynamics of the single-field model (review Eq. 16, no noise):
#
#     dphi/dt = M [ eps^2 phi_xx - w g'(phi) - p'(phi) Delta_f ].
#
# A planar interface between crystal (phi = 1, left) and liquid (phi = 0, right) moves into the
# liquid at a constant velocity once the profile has adjusted:
#
#     v = -6 M Delta_f delta = 6 x 0.4923638 x 0.227560 = 0.672254 nm per time unit  (closed form:
#     with this p the travelling profile is the equilibrium logistic, review Eq. 16 ff.)
#
# The velocity is read from the solid length: `mass` of solver_history.csv is int phi dx, and
# v = d(mass)/dt (fit over t > 5, after the start-up transient).
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa; mobility M = 1, so time is
# in units of 1/(M GPa).
#
# PARAMETERS at T = 0.8 T_f: eps^2 = 0.329871 aJ/nm, w = 12.740387 GPa, Delta_f = -0.4923638 GPa,
# delta = 0.227560 nm.
#
# EXPECTED RESULT. mass(t) is a straight line for t > 5; the free energy falls linearly (rate
# v Delta_f); the profile keeps its shape. With bdf2 and dt = 0.02 the velocity error is below 1 %.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a2_front_velocity_tanh_fem
[]

[Problem]
  type          = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type      = line_mesh
  domain    = '0.0 40.0'       # nm
  divisions = 4000             # h = 0.01 nm = delta / 23
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
    position  = 5.0
    width     = 0.91024
    amplitude = 1.0
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**2*(3-2*phi)*df'
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
  type            = pf_fem_transient
  scheme          = bdf2
  n_steps         = 1750
  dt              = 0.02       # to t = 35: the front travels about 23.5 nm
  newton_tol      = 1.0e-10
  snapshots_every = 25
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
