# Part IV, case B1 (FE route): heterogeneous crystal nucleus on a foreign wall, Model I,
# contact angle Psi = 61.2 deg, undercooling 70 K (review Fig. 6, upper row and solid line).
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
#       = -0.0117278 at 70 K.
# Check: Model II's layer loses stability at D = -1/64, i.e. at 92.0 K, the value the papers state.
#
# HOW IT IS SOLVED. Half of the symmetric nucleus on [0, 25] x [0, 20]: wall at y = 0,
# symmetry (natural) at x = 0, liquid (phi = 0) on the two far sides. The nucleus is a
# saddle point whose Newton basin is narrow, so `pf_fem_saddle` is used: Newton at a fixed
# amount of crystal (int phi dA, where the nucleus is a minimum), the amount moved until the
# constraint multiplier vanishes, then a plain Newton polish. It starts from the circular cap
# of the classical radius R = gamma / |Delta_f| = 7.1056 with its centre 3.4232 below the wall.
#
# EXPECTED RESULT. W = 2 (free_energy + wall_energy) of postprocessors.csv (the factor 2 is
# the mirror half). Classical cap: W_cap = 0.38251. Model I follows the cap closely at small
# undercooling and falls below it as the nucleus shrinks toward the interface width
# (review Fig. 6, bottom panel, solid line); the measured ratio is in the report.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/walls/b1_wall_model1_dT070_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 25 20'
  divisions = '100 80'       # h = 0.25 delta
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
    expression = 0.5*(1-tanh((sqrt(x**2+(y+3.4232)**2)-7.1056)/2))
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
      df = -0.0117278
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
  type            = pf_fem_saddle
  scheme          = bdf1
  n_steps         = 1
  dt              = 1.0e12
  max_newton_iter = 60
  newton_tol      = 1.0e-9
  max_cutbacks    = 0
  snapshots_every = 1
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
