# Part IV, case F3 (FE route): ternary eutectic growth at a fixed driving force, two-dimensional
# section. Free energy of the pure solids f_s = -0.86 (liquid: 0), liquid composition
# (c1, c2, c3) = (0.455, 0.455, 0.09), initial lamellar spacing 16 cells.
#
# PHYSICS. Ternary phase-field model of (ratkai 2015 ternary), dimensionless; phi = 1 liquid,
# phi = 0 solid (the paper uses the opposite labels), c3 = 1 - c1 - c2:
#
#     F = int { kappa_phi/2 |grad phi|^2 + w g(phi) + p(phi) f_l + [1 - p(phi)] [f_s + kappa_c/2 sum_i |grad c_i|^2] } dV
#     f_l = sum_i c_i ln c_i,   f_s = sum_i c_i (f_s,i + ln c_i) + Omega_12 c1 c2 + Omega_23 c2 c3 + Omega_31 c3 c1
#     phi_t = -M_phi { w g' + p' [f_l - f_s - kappa_c/2 sum_i |grad c_i|^2] - kappa_phi lap phi }
#     c_i,t = sum_j div[ p(phi) M_ij grad mu_j ],   M = 1 on the diagonal, -1/2 off it
#     mu_j  = ln c_j + [1 - p] (f_s,j + sum_k Omega_jk c_k) - kappa_c div([1 - p] grad c_j)
#
# Components 1 and 2 separate into two solid solutions (Omega_12 = 3.05 > 2); component 3 is
# soluble in neither (Omega_23 = Omega_31 = 3.0) and is rejected into the liquid, where it
# piles up ahead of the front and, in three dimensions, destabilizes it into two-phase dendrites.
#
# PARAMETERS (reported_text, the paper's reference conditions): kappa_phi = 0.75, kappa_c = 0.4,
#   w = 0.0469, M_phi = 1, Omega_12 = 3.05, Omega_23 = Omega_31 = 3.0, dx = 1, dt = 0.0025,
#   liquid composition (0.455, 0.455, 0.09). The paper pulls the sample at a fixed velocity
#   through a gradient of f_s (-0.9640 at the cold end, -0.9140 at the tip) in a 96 x 96 x 612
#   box; here the driving force is uniform and the section two-dimensional (`assumed` set-up).
#   assumed: the initial lamellae (solid slab of 40 cells, compositions 0.8225 / 0.0875 / 0.09).
#
# HOW IT IS SOLVED. `pf_fem_explicit`; the three chemical potentials are auxiliary variables
# (lumped projections), so the fourth-order composition equations need no third derivative.
# Periodic across the lamellae, no flux at the bottom and top.
#
# EXPECTED RESULT. See the Part IV report (ternary eutectic section): a lamellar front that
# rejects component 3; with rising driving force the front speeds up and traps it.

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/ternary/f3_long_fs-0.86_c30.09_lam16_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 64 192'
  divisions = '64 192'
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
  [./c1]
    family     = fem
    components = 1
    order      = 2
  []
  [./c2]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu1]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu2]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu3]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # solid below y = 40; interface profile of width sqrt(2 kappa_phi / w) = 5.66
  [./front]
    kind       = parsed
    expression = 0.5*(1+tanh((y-40.0)/11.31))
  []
  [./comp1]
    kind       = parsed
    expression = 0.455+0.5*(1-tanh((y-40.0)/11.31))*0.3675*tanh(cos(2*pi*x/16)/0.3)
  []
  [./comp2]
    kind       = parsed
    expression = 0.455-0.5*(1-tanh((y-40.0)/11.31))*0.3675*tanh(cos(2*pi*x/16)/0.3)
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = front
  []
  [./c10]
    variable = c1
    kind     = from_function
    function = comp1
  []
  [./c20]
    variable = c2
    kind     = from_function
    function = comp2
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
    mobility             = 1.0
    gradient_coefficient = 0.75
    anisotropy_strength  = 0.0
  []
  [./phi_bulk]
    type                             = ts_phase_bulk
    variable                         = phi
    well                             = 0.0469
    mobility                         = 1.0
    composition_gradient_coefficient = 0.4
    compositions                     = 'c1 c2'
    solid_energies                   = '-0.86 -0.86 -0.86'
    interactions                     = '3.05 3.0 3.0'
  []
  [./dc1_dt]
    type        = time_derivative
    variable    = c1
    coefficient = 1.0
  []
  [./c1_flux]
    type                = ts_solute_flux
    variable            = c1
    phase               = phi
    chemical_potentials = 'mu1 mu2 mu3'
    mobilities          = '1.0 -0.5 -0.5'
  []
  [./dc2_dt]
    type        = time_derivative
    variable    = c2
    coefficient = 1.0
  []
  [./c2_flux]
    type                = ts_solute_flux
    variable            = c2
    phase               = phi
    chemical_potentials = 'mu1 mu2 mu3'
    mobilities          = '-0.5 1.0 -0.5'
  []
  [./potential1]
    type                 = ts_chemical_potential
    variable             = mu1
    phase                = phi
    component            = 1
    gradient_coefficient = 0.4
    compositions         = 'c1 c2'
    solid_energies       = '-0.86 -0.86 -0.86'
    interactions         = '3.05 3.0 3.0'
  []
  [./potential2]
    type                 = ts_chemical_potential
    variable             = mu2
    phase                = phi
    component            = 2
    gradient_coefficient = 0.4
    compositions         = 'c1 c2'
    solid_energies       = '-0.86 -0.86 -0.86'
    interactions         = '3.05 3.0 3.0'
  []
  [./potential3]
    type                 = ts_chemical_potential
    variable             = mu3
    phase                = phi
    component            = 3
    gradient_coefficient = 0.4
    compositions         = 'c1 c2'
    solid_energies       = '-0.86 -0.86 -0.86'
    interactions         = '3.05 3.0 3.0'
  []
[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = 600000
  dt                  = 0.0025
  periodic_axes       = 'x'
  auxiliary_variables = 'mu1 mu2 mu3'
  snapshots_every     = 30000
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c1 c2'
[]
