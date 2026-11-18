# Part IV, case F2 (FE route): lamellar eutectic growth of Ag-Cu at T = 955.5 K
# (undercooling 100.0 K below the eutectic temperature 1055.5 K), liquid composition
# c0 = 0.345 (mole fraction of Cu), lamellar spacing lambda = 8 nm.
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
# omega_S > 2: the solid has a miscibility gap (c_alpha = 0.0938, c_beta = 0.9062 at this
# temperature); the liquid solidifies into alternating lamellae of the two solid solutions,
# exchanging Cu sideways through the liquid ahead of the front. The gradient energy of c gives
# the lamella boundaries their energy; the composition equation is fourth order.
#
# UNITS. Length xi = 1 nm, time xi^2/D_L = 7.143e-10 s, energy density R T/v = 9.132e+08 J/m^3.
#
# PARAMETERS (provenance: docs/14_phase_field_crystallization_plan.md, sections 9.3, 9.4;
# computed by scripts/cryst_eutectic_thermo.py).
#   reported_text: L_Ag = 1097 J/cm^3, L_Cu = 1728 J/cm^3, v = 8.7 cm^3/mol, T_Ag = 1234 K,
#     T_Cu = 1358 K, Omega_L v = 16,000 - 4 T, Omega_S v = 26,000 - 4 T J/mol; gamma = 172 (Ag),
#     227 (Cu) mJ/m^2; 10-90 % interface thickness 1 nm; eps_c^2 = eps_phi^2/4; D_L = 1.4e-9 m^2/s;
#     M_phi = 0.25 m^3/(J s); D_S = D_L/100.
#   derived: f_A = -0.271126, f_B = -0.560873, omega_S = 2.79163, omega_L = 1.53289,
#     kappa_phi = 0.199134, kappa_c = 0.049783, w_A = 7.69102, w_B = 11.06139, m_phi = 0.163064.
#   assumed: the initial lamella pair (solid slab of 10 nm with the lever-rule fraction
#     0.309 of the Cu-rich phase), mesh size 0.25 nm, time step 0.004.
#
# HOW IT IS SOLVED. `pf_fem_explicit` (forward Euler, lumped mass, bilinear elements); mu is an
# auxiliary variable: the lumped projection of its kernel terms each step, so no fourth
# derivative is formed. Periodic across the lamellae, no flux at the bottom and top.
#
# EXPECTED RESULT. The two-phase front advances at a steady speed that depends on lambda (a
# maximum at an intermediate spacing; too narrow a pair is overgrown by one phase).
# The front speed is the slope of the solid fraction against time times the box length:
# `phi_mean` of solver_history.csv. Report: Part IV, eutectic section.

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/eutectic/f2_lamellar_T955.5_c0.345_lam8_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 8 60'      # nm: one lamella pair across, growth along y
  divisions = '32 240'                   # 0.25 nm
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
[]

[Functions]
  # solid below y = 10 nm; the Cu-rich lamella occupies the central fraction of the width
  [./front]
    kind       = parsed
    expression = 0.5*(1+tanh((y-10.0)/0.4551))
  []
  [./composition]
    kind       = parsed
    expression = 0.345*0.5*(1+tanh((y-10.0)/0.4551))+0.5*(1-tanh((y-10.0)/0.4551))*(0.093829+0.812343*0.5*(1+tanh((cos(2*pi*(x/8-0.5))-(0.564176))/0.15)))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = front
  []
  [./c0]
    variable = c
    kind     = from_function
    function = composition
  []
  [./mu0]
    variable = mu
    kind     = constant
    value    = 0.0
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
    mobility             = 0.163064
    gradient_coefficient = 0.199134
    anisotropy_strength  = 0.0
  []
  [./phi_bulk]
    type          = rs_phase_bulk
    variable      = phi
    concentration = c
    mobility      = 0.163064
    solid_a       = -0.271126
    solid_b       = -0.560873
    omega_solid   = 2.79163
    omega_liquid  = 1.53289
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
    solid_a              = -0.271126
    solid_b              = -0.560873
    omega_solid          = 2.79163
    omega_liquid         = 1.53289
    well_a               = 7.69102
    well_b               = 11.06139
  []
[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = 300000
  dt                  = 0.004
  periodic_axes       = 'x'
  auxiliary_variables = 'mu'
  snapshots_every     = 15000
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c mu'
[]
