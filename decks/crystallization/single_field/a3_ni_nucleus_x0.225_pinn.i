# Part IV, case A3 (PINN route): critical nucleus of crystalline Ni in its undercooled melt,
# reduced undercooling Delta T / T_f = 0.225 (T = 1339.43 K).
#
# PHYSICS. Single-field phase-field model of (granasy 2019 nucleation), Eq. 12-15;
# phi = 0 liquid, phi = 1 crystal:
#
#     F[phi] = int { eps^2/2 |grad phi|^2 + w g(phi) + p(phi) Delta_f } dV
#     g = phi^2 (1 - phi)^2 / 4,   p = phi^3 (10 - 15 phi + 6 phi^2)
#
# The critical nucleus is the stationary point (a saddle) of F with phi -> 0 far away:
#
#     eps^2 (phi'' + 2/r phi') = w g'(phi) + p'(phi) Delta_f,   phi'(0) = 0.
#
# Its excess free energy is the nucleation barrier W* = F[phi*] (the liquid has F = 0).
#
# HOW IT IS SOLVED. One network phi_theta(r) (no time input) is fitted to the initial
# condition (the classical droplet, R = 2 gamma / |Delta_f|) and then trained on the squared residual of the
# stationary equation above (`pf_pinn_stationary`): a saddle point of F is a zero, hence a
# minimum, of that loss, so the training can hold the state no time march can. The network
# reads r^2, so phi'(0) = 0 holds for any weights. The trivial liquid phi = 0 is a zero of the
# residual too: the deflation (`deflation_power`) removes it from the minima of the loss, and
# `collapse_tol` stops a run that fell onto it anyway.
# The mesh is used for the output and the energy integral only (spherical volume element).
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa, mobility 1.
#
# PARAMETERS (provenance in docs/14_phase_field_crystallization_plan.md, section 5).
#   Ni, review Table 2 (reported_table): gamma = 0.302 J/m^2, 10-90 % width 1 nm,
#   H_f = 17.29 kJ/mol, T_f = 1728.3 K, rho = 8357 kg/m^3, molar mass 58.6934 g/mol.
#   derived at T_f: eps^2 = 0.412338 aJ/nm, w = 15.92548 GPa, delta = sqrt(2 eps^2/w) = 0.22756 nm
#   derived at T = 0.775 T_f (eps^2 and w proportional to T):
#     eps^2 = 0.319562 aJ/nm,  w = 12.342250 GPa
#     Delta_f = -H_f (T_f - T)/(T_f v_m) = -0.5539093 GPa
#
# EXPECTED RESULT (reference: scripts/cryst_reference.py, shooting solution).
#   W* = 0.64302 aJ = 34.771 kT (kT = 0.0184929 aJ); phi(0) = 0.85505;
#   equimolar radius 1.0071 nm.
#   In results: `w_star_aj` of postprocessors.csv is W* in aJ. Seed 0: see the report for the
#   spread over seeds.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a3_ni_nucleus_x0.225_pinn
[]

[Problem]
  type      = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type              = line_mesh
  domain            = '0.0 6.0'       # nm; five nucleus radii
  divisions         = 1200            # output and quadrature mesh, h = 0.005 nm
  coordinate_system = spherical       # dV = 4 pi r^2 dr; r = 0 is the natural (symmetry) end
[]

[Variables]
  [./phi]
    family     = pinn
    components = 1
    order      = 2
  []
[]

[Functions]
  # rough guess: classical droplet, R = 2 gamma(T) / |Delta_f| = 0.8451 nm, amplitude 1
  [./droplet]
    kind       = parsed
    expression = 0.5*(1-tanh((x-0.8451)/(2*0.22756)))   # unquoted, no blanks
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = droplet
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'
    [./coefficients]
      w  = 12.342250
      df = -0.5539093
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
    gradient_coefficient = 0.319562
    spatial_dim          = 1
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.319562    # eps^2(T)
    spatial_dim          = 1
  []
[]

[BCs]
  # far field: undercooled liquid
  [./liquid]
    variable = phi
    kind     = dirichlet
    node_set = xmax
    components = [0]
    values     = [0.0]
  []
[]

[Executioner]
  type                 = pf_pinn_stationary
  hidden_layers        = 3
  neurons              = 32
  activation           = tanh
  collocation_interior = 600        # even grid in r
  pretrain_epochs      = 3000       # regression on the droplet
  pretrain_lr          = 5.0e-3
  lbfgs_iters          = 150        # x 20 inner iterations on the residual loss
  w_physics            = 1.0
  w_dirichlet          = 10.0
  deflation_power      = 1.0        # the liquid phi = 0 is not a minimum of the loss
  collapse_tol         = 0.1
[]

[Postprocessors]
  [./W_star_aJ]
    type          = total_free_energy
    bulk_field    = field_phi
    free_energy   = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.319562
      []
    []
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
