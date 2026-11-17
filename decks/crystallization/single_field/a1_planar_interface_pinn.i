# Part IV, case A1 (PINN route): planar crystal-liquid interface of Ni at the melting point.
#
# PHYSICS. Single-field model of (granasy 2019 nucleation), Eq. 12 and 16, with Delta_f = 0:
#
#     dphi/dt = M [ eps^2 phi_xx - w g'(phi) ],   g = phi^2 (1 - phi)^2 / 4.
#
# A diffuse interface that starts too wide (tanh width twice the equilibrium one) relaxes to
# the stationary profile. First integral eps^2/2 phi_x^2 = w g gives the closed forms
#
#     phi(x) = 1 / (1 + exp(x / delta)),   delta = sqrt(2 eps^2 / w) = 0.227560 nm
#     gamma  = eps^2 / (6 delta)                                      = 0.302 J/m^2
#     d(10-90 %) = 2 ln(9) delta                                      = 1.000 nm
#
# which are the two material inputs (review Table 2) the model parameters were derived from:
# the run checks that PINNACLE returns them.
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa; mobility M = 1, so time is
# in units of 1/(M GPa).
#
# PARAMETERS. eps^2 = 0.412338 aJ/nm, w = 15.925484 GPa (derived from gamma and the width).
#
# HOW IT IS SOLVED. `pf_pinn_stationary`: one network phi_theta(x) is fitted to the too-wide
# initial profile and trained on the squared stationary residual eps^2 phi_xx - w g'(phi).
# The bulk phases are pinned at the two ends (Dirichlet), so the uniform states are excluded
# and no deflation is needed. The position of the interface is not fixed by the equation on
# an infinite line (translation invariance); here it stays where the guess put it, at x = 0.
#
# EXPECTED RESULT. `free_energy` of postprocessors.csv = gamma = 0.3020 aJ/nm^2 within 1e-3;
# the profile is the logistic within 1e-3.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a1_planar_interface_pinn
[]

[Problem]
  type          = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type      = line_mesh
  domain    = '-3.0 3.0'       # nm; phi is within 2e-6 of the bulk values at the ends
  divisions = 600              # h = 0.01 nm = delta / 23
[]

[Variables]
  [./phi]
    family     = pinn
    components = 1
    order      = 2
  []
[]

[ICs]
  [./phi0]
    variable  = phi
    kind      = tanh_interface   # 0.5 A (1 - tanh(2 (x - position)/width)); equilibrium width = 4 delta = 0.91
    axis      = 0
    position  = 0.0
    width     = 1.82
    amplitude = 1.0
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4'
    [./coefficients]
      w = 15.925484
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
    gradient_coefficient = 0.412338
    spatial_dim          = 1
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.412338    # eps^2(T)
    spatial_dim          = 1
  []
[]

[BCs]
  # bulk crystal on the left, bulk liquid on the right
  [./crystal]
    variable   = phi
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [1.0]
  []
  [./liquid]
    variable   = phi
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [0.0]
  []
[]

[Executioner]
  type                 = pf_pinn_stationary
  hidden_layers        = 3
  neurons              = 32
  collocation_interior = 600
  pretrain_epochs      = 2000
  lbfgs_iters          = 150
  w_dirichlet          = 10.0
[]

[Postprocessors]
  [./free_energy]
    type          = total_free_energy
    bulk_field    = field_phi
    free_energy   = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.412338
      []
    []
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
