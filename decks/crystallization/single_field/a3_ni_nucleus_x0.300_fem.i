# Part IV, case A3 (FE route): critical nucleus of crystalline Ni in its undercooled melt,
# reduced undercooling Delta T / T_f = 0.3 (T = 1209.81 K).
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
# HOW IT IS SOLVED. One backward-Euler step of the Allen-Cahn flow with a very large time
# step is Newton's method on the stationary equation, started from the initial condition
# (a tanh droplet at the equimolar radius of the reference solution). The radial line carries
# the spherical volume element 4 pi r^2 ([Mesh] coordinate_system = spherical), so the weak
# form is the spherically symmetric one and the energy monitor returns W* directly.
# The nucleus is unstable under the flow itself: a time march from this droplet either
# melts or grows (case A5), so the stationary solve is the only way to hold it.
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa, mobility 1.
#
# PARAMETERS (provenance in docs/14_phase_field_crystallization_plan.md, section 5).
#   Ni, review Table 2 (reported_table): gamma = 0.302 J/m^2, 10-90 % width 1 nm,
#   H_f = 17.29 kJ/mol, T_f = 1728.3 K, rho = 8357 kg/m^3, molar mass 58.6934 g/mol.
#   derived at T_f: eps^2 = 0.412338 aJ/nm, w = 15.92548 GPa, delta = sqrt(2 eps^2/w) = 0.22756 nm
#   derived at T = 0.7 T_f (eps^2 and w proportional to T; gamma and the
#   10-90 % width matched for this well):
#     eps^2 = 0.288637 aJ/nm,  w = 11.147839 GPa
#     Delta_f = -H_f (T_f - T)/(T_f v_m) = -0.7385457 GPa
#
# EXPECTED RESULT (reference: scripts/cryst_reference.py, shooting solution).
#   W* = 0.24880 aJ = 14.895 kT (kT = 0.0167032 aJ); phi(0) = 0.70132;
#   equimolar radius 0.7991 nm.
#   In results: `energy` of solver_history.csv is W* in aJ.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a3_ni_nucleus_x0.300_fem
[]

[Problem]
  type      = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type              = line_mesh
  domain            = '0.0 12.0'      # nm; ten nucleus radii
  divisions         = 2400            # h = 0.005 nm = delta / 45
  coordinate_system = spherical       # dV = 4 pi r^2 dr; r = 0 is the natural (symmetry) end
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # initial guess of the Newton solve: a droplet close to the nucleus (a saddle point needs a
  # near start: amplitude, radius and width within about 10 % of the solution)
  [./droplet]
    kind       = parsed
    expression = 0.7380*0.5*(1-tanh((x-0.6273)/(2*0.19460)))   # unquoted, no blanks
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
      w  = 11.147839
      df = -0.7385457
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
    gradient_coefficient = 0.288637
    spatial_dim          = 1
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.288637    # eps^2(T)
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
  type            = pf_fem_transient
  scheme          = bdf1
  n_steps         = 1
  dt              = 1.0e12     # stationary solve: the rate term is 1e-12 of the others
  max_newton_iter = 50
  newton_tol      = 1.0e-10
  max_cutbacks    = 0
  snapshots_every = 1
[]

[Postprocessors]
  [./W_star_aJ]
    type          = total_free_energy
    bulk_field    = field_phi
    free_energy   = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.288637
      []
    []
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
