# Part IV, case B1 (PINN route): heterogeneous crystal nucleus on a foreign wall, Model I,
# contact angle Psi = 61.2 deg, undercooling 40 K (review Fig. 6, upper row and solid line).
#
# PHYSICS. Single-field model with a wall energy (granasy 2019 nucleation, Eq. 17-19;
# granasy 2007 heterogeneous), phi = 0 liquid, phi = 1 crystal, two dimensions:
#
#     F = int { eps^2/2 |grad phi|^2 + w g(phi) + p(phi) Delta_f } dA + int_wall Z(phi) ds
#     g = phi^2 (1 - phi)^2 / 4,   p = phi^3 (10 - 15 phi + 6 phi^2)
#
# Model I: the wall does not disturb the planar interface profile, which gives the natural
# condition (outward normal n)
#
#     eps^2 n.grad(phi) = -Z'(phi) = 6 gamma cos(Psi) phi (1 - phi),
#     Z = -gamma cos(Psi) (3 phi^2 - 2 phi^3),   gamma = eps^2 / (6 delta),
#
# i.e. Young's equation gamma_WL - gamma_WS = gamma cos(Psi) with sharp wall-liquid and
# wall-crystal interfaces. The critical nucleus is the stationary state (a saddle) with
# liquid far away; its excess free energy W = F is the barrier. The sharp-interface
# counterpart is the circular cap: W_cap = (gamma^2 / |Delta_f|) (Psi - sin Psi cos Psi).
#
# UNITS AND PARAMETERS. Scaled: lengths in delta = sqrt(2 eps^2 / w), energies per unit
# length in w(T) delta^2; eps^2 = 0.5, w = 1, gamma = 1/12. Ni data of the original paper
# (granasy 2007 heterogeneous, arXiv cond-mat/0608375, footnote 16; `reported_text`):
# gamma_SL = 0.364 J/m^2, interface thickness 2 nm, L = 2.666e9 J/m^3, T_m = 1728 K, eps^2 and
# w multiplied by T. Derived: delta = 2 nm / (2 ln 9) = 0.4551 nm, w(T_m) = 12 gamma/delta
# = 9.598 GPa, and
#     D = Delta_f / w(T) = -(L / w(T_m)) Delta T / (T_m - Delta T) = -0.277781 Delta T / (1728 K - Delta T)
#       = -0.0065825 at 40 K.
# Check: Model II's layer loses stability at D = -1/64, i.e. at 92.0 K, the value the papers state.
#
# HOW IT IS SOLVED (PINN route, `pf_pinn_stationary residual_form = fem`). A network
# phi_theta(x, y) (4 x 48, tanh) is fitted to the circular cap of the classical radius and then
# trained (L-BFGS) on the discrete L2 norm of the weak residual the FE arm assembles for its
# nodal values on the mesh of the FE deck: loss = sum_i R_i^2 / m_i / sum_i m_i (m_i: nodal
# measure), with the exact adjoint of the consistent Jacobian. The wall condition, the symmetry
# plane and the Dirichlet values are inside R; no strong form and no boundary penalty is used.
# The nucleus is a saddle of F but a zero, hence a minimum, of the loss; the state without a
# nucleus is deflated (the liquid phi = 0).
#
# EXPECTED RESULT. The free energy of the FE deck (b1_wall_model1_dT040_fem.i) within 1e-3 and the field
# within about 0.02 of it; values and seed scatter in the report.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/walls/b1_wall_model1_dT040_pinn_seed0
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 40 25'
  divisions = '160 100'       # h = 0.25 delta
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # circular cap of the classical radius (unquoted, no blanks)
  [./cap]
    kind       = parsed
    expression = 0.5*(1-tanh((sqrt(x**2+(y+6.0989)**2)-12.6598)/2))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = cap
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'
    [./coefficients]
      w  = 1.0
      df = -0.0065825
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
    gradient_coefficient = 0.5
    spatial_dim          = 2
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.5
    spatial_dim          = 2
  []
[]

[BCs]
  # Model I wall: eps^2 n.grad(phi) = 6 gamma cos(Psi) phi (1 - phi)
  [./wall]
    kind                 = pf_contact_angle
    variable             = phi
    side_set             = ymin
    contact_angle        = 61.2
    surface_tension      = 0.0833333333333   # gamma = 1/12
    gradient_coefficient = 0.5
    mobility             = 1.0
    phase_values         = '0 1'
    flux_form            = young
  []
  [./liquid_right]
    variable   = phi
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [0.0]
  []
  [./liquid_top]
    variable   = phi
    kind       = dirichlet
    node_set   = ymax
    components = [0]
    values     = [0.0]
  []
[]

[Executioner]
  type               = pf_pinn_stationary
  residual_form      = fem         # loss: discrete L2 norm of the assembled weak residual
  hidden_layers      = 4
  neurons            = 48
  pretrain_epochs    = 3000        # regression on the circular cap (the initial condition)
  lbfgs_iters        = 150
  deflation_power    = 1.0         # the state without a nucleus is not a minimum of the loss
  collapse_tol       = 0.1
  seed               = 0
[]

[Postprocessors]
  [./free_energy]
    type        = total_free_energy
    bulk_field  = field_phi
    free_energy = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.5
      []
    []
  []
  [./wall_energy]
    type = wall_free_energy
    bc   = wall
  []
  [./angle]
    type     = avg_contact_angle
    variable = phi
    side_set = ymin
    band     = 0.4
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
