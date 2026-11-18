# Part IV, case F2b (FE route): growth of a eutectic grain of Ag-Cu at T = 900 K
# (undercooling 155.5 K) from a two-phase (lamellar) seed, liquid composition c0 = 0.345,
# with an orientation field: the two solid solutions keep a fixed orientation relation,
# after (lewis 2004 eutectic).
#
# PHYSICS. The binary regular-solution model of case F2 (equations below) plus the orientation
# theta in [0, 1) with the composition-dependent energy
#
#     f_ori = [1 - p(phi)] S E(c, u),   u = d |grad theta|
#     E = h(c) u + [1 - h(c)] {A + B |cos(2 pi m u)|},   h = (1/2){1 + cos[2 pi (c - c_alpha)/(c_beta - c_alpha)]}
#
# Inside either solid (h = 1) E is the grain-boundary energy, linear in the misorientation;
# across a boundary between the two solids (h -> 0) it is lowest at u = 1/(4 m) = 1/2:
# a step of half a period in theta. theta_t = chi(phi) div{dE-flux}, chi = chi_l p(phi).
# reported_text (lewis 2004 eutectic): A = 0.9, B = 0.1, m = 1/2, psi = 0, d = dx.
# assumed (not printed anywhere): S = 0.1, chi_l = 3, no noise, the seed.
#
#
# PHYSICS. Binary phase-field model with regular-solution thermodynamics, (lewis 2004 eutectic)
# and (ratkai 2017 nanoscale); phi = 1 liquid, phi = 0 solid, c = Cu fraction; energy in R T/v:
#
#     F = int { kappa_phi/2 |grad phi|^2 + kappa_c/2 |grad c|^2 + w(c) g(phi)
#               + [1 - p(phi)] f_S(c) + p(phi) f_L(c) } dV
#     f_L = omega_L c (1 - c) + c ln c + (1 - c) ln(1 - c)
#     f_S = (1 - c) f_A + c f_B + omega_S c (1 - c) + c ln c + (1 - c) ln(1 - c)
#     phi_t = -m_phi { w(c) g' + p' [f_L - f_S] - kappa_phi lap phi }
#     c_t   = div{ D(phi) c (1 - c) grad mu },   D = r + (1 - r) p(phi)
#     mu    = (w_B - w_A) g + [1 - p] f_S' + p f_L' - kappa_c lap c
#
# omega_S > 2: the solid has a miscibility gap (c_alpha = 0.0713, c_beta = 0.9287 at this
# temperature); the liquid solidifies into alternating lamellae of the two solid solutions,
# exchanging Cu sideways through the liquid ahead of the front. The gradient energy of c gives
# the lamella boundaries their energy; the composition equation is fourth order.
#
# UNITS. Length xi = 1 nm, time xi^2/D_L = 7.143e-10 s, energy density R T/v = 8.601e+08 J/m^3.
#
# PARAMETERS (provenance: docs/14_phase_field_crystallization_plan.md, sections 9.3, 9.4;
# computed by scripts/cryst_eutectic_thermo.py).
#   reported_text: L_Ag = 1097 J/cm^3, L_Cu = 1728 J/cm^3, v = 8.7 cm^3/mol, T_Ag = 1234 K,
#     T_Cu = 1358 K, Omega_L v = 16,000 - 4 T, Omega_S v = 26,000 - 4 T J/mol; gamma = 172 (Ag),
#     227 (Cu) mJ/m^2; 10-90 % interface thickness 1 nm; eps_c^2 = eps_phi^2/4; D_L = 1.4e-9 m^2/s;
#     M_phi = 0.25 m^3/(J s); D_S = D_L/100.
#   derived: f_A = -0.345208, f_B = -0.677567, omega_S = 2.99345, omega_L = 1.65709,
#     kappa_phi = 0.199134, kappa_c = 0.049783, w_A = 7.69102, w_B = 11.06139, m_phi = 0.153592.
#   assumed: the initial lamella pair (solid slab of 10 nm with the lever-rule fraction
#     0.000 of the Cu-rich phase), mesh size 0.25 nm, time step 0.004.
#
# HOW IT IS SOLVED. `pf_fem_explicit` (forward Euler, lumped mass, bilinear elements); mu is an
# auxiliary variable: the lumped projection of its kernel terms each step, so no fourth
# derivative is formed. Periodic across the lamellae, no flux at the bottom and top.
#
# EXPECTED RESULT. The lamellae of the seed grow outward side by side; new Cu-rich solid
# takes an orientation half a period from the Ag-rich solid next to it.
# Report: Part IV, eutectic section.

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/eutectic/grain/f2b_grain_T900_c0.345_S0.1_chi3_lam10_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '-32 -32 32 32'      # nm
  divisions = '256 256'                   # 0.25 nm
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
  [./c]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu]
    family     = fem
    components = 1
    order      = 2
  []
  [./theta]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # two-phase seed: a disc of radius 10 nm cut out of lamellae of spacing 10 nm along x
  # (Cu-rich fraction 0.319, the lever rule)
  [./seed_phi]
    kind       = parsed
    expression = 0.5*(1+tanh((sqrt(x**2+y**2)-10)/0.4551))
  []
  [./seed_c]
    kind       = parsed
    expression = (0.071334+0.857332*0.5*(1+tanh((cos(2*pi*x/10)-0.537930)/0.3)))*(1-0.5*(1+tanh((sqrt(x**2+y**2)-10)/0.4551)))+0.345*0.5*(1+tanh((sqrt(x**2+y**2)-10)/0.4551))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = seed_phi
  []
  [./c0]
    variable = c
    kind     = from_function
    function = seed_c
  []
  [./mu0]
    variable = mu
    kind     = constant
    value    = 0.0
  []
  # no orientational order in the liquid; the seed has theta = 0.25
  [./theta0]
    variable = theta
    kind     = uniform_random
    minimum  = 0.0
    maximum  = 1.0
    seed     = 20261118
    center   = '0 0'
    radius   = 10
    inside   = 0.25
  []
[]

[Kernels]
  [./dphi_dt]
    type        = time_derivative
    variable    = phi
    coefficient = 1.0
  []
  [./phi_gradient]
    type                 = ac_interface_anisotropic
    variable             = phi
    mobility             = 0.153592
    gradient_coefficient = 0.199134
    anisotropy_strength  = 0.0
  []
  [./phi_bulk]
    type          = rs_phase_bulk
    variable      = phi
    concentration = c
    mobility      = 0.153592
    solid_a       = -0.345208
    solid_b       = -0.677567
    omega_solid   = 2.99345
    omega_liquid  = 1.65709
    well_a        = 7.69102
    well_b        = 11.06139
  []
  [./dc_dt]
    type        = time_derivative
    variable    = c
    coefficient = 1.0
  []
  [./solute]
    type               = rs_solute_flux
    variable           = c
    phase              = phi
    chemical_potential = mu
    diffusivity        = 1.0
    diffusivity_ratio  = 0.01
  []
  [./potential]
    type                 = rs_chemical_potential
    variable             = mu
    phase                = phi
    concentration        = c
    gradient_coefficient = 0.049783
    solid_a              = -0.345208
    solid_b              = -0.677567
    omega_solid          = 2.99345
    omega_liquid         = 1.65709
    well_a               = 7.69102
    well_b               = 11.06139
  []
  [./phi_orientation]
    type          = eo_phase_coupling
    variable      = phi
    orientation   = theta
    concentration = c
    mobility      = 0.153592
    strength            = 0.1        # S (assumed)
    boundary_thickness  = 0.25       # d = dx
    solid_composition_a = 0.071334
    solid_composition_b = 0.928666
    offset              = 0.9
    amplitude           = 0.1
    misorientation_index = 0.5
    regularization      = 1.0e-3
  []
  [./mu_orientation]
    type          = eo_potential_coupling
    variable      = mu
    phase         = phi
    orientation   = theta
    concentration = c
    strength            = 0.1        # S (assumed)
    boundary_thickness  = 0.25       # d = dx
    solid_composition_a = 0.071334
    solid_composition_b = 0.928666
    offset              = 0.9
    amplitude           = 0.1
    misorientation_index = 0.5
    regularization      = 1.0e-3
  []
  [./dtheta_dt]
    type        = time_derivative
    variable    = theta
    coefficient = 1.0
  []
  [./theta_flux]
    type             = eo_grain_flux
    variable         = theta
    phase            = phi
    concentration    = c
    mobility_liquid  = 3       # chi_l (assumed)
    mobility_crystal = 0.0
    strength            = 0.1        # S (assumed)
    boundary_thickness  = 0.25       # d = dx
    solid_composition_a = 0.071334
    solid_composition_b = 0.928666
    offset              = 0.9
    amplitude           = 0.1
    misorientation_index = 0.5
    regularization      = 1.0e-3
  []
[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = 100000
  dt                  = 0.004
  periodic_axes       = 'x y'
  angular_variables   = 'theta'
  subcycles           = 'theta 10'
  auxiliary_variables = 'mu'
  snapshots_every     = 5000
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c theta'
[]
