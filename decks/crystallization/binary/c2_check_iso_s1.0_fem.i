# Part IV, case C2 (FE route), small check: growth of a single crystal in the binary Cu-Ni
# model at S = 1.0, isotropic, no noise, 200 x 200 cells; compared field by field with the
# independent finite-difference reference (decks/crystallization/reference/c2_check_iso_s1.0.i).
#
# PHYSICS. Binary ideal-solution phase-field model of (granasy 2005 spherulites), Eq. 2.3-2.4,
# orientation uniform; phi = 1 liquid, phi = 0 crystal, c = Cu fraction:
#
#     phi_t = m_phi { lap(phi) - a_w(c) g'(phi) - p'(phi) b(c) }
#     c_t   = div{ p(phi) [ grad c - c (1 - c) (q_f p'(phi) + q_w g'(phi)) grad phi ] }
#
# UNITS. Lengths in xi = 2.1e-4 cm, time in xi^2/D_l (dimensionless form of the paper).
#
# PARAMETERS (provenance: docs/14_phase_field_crystallization_plan.md, section 7.2).
#   a_w = xi^2 w / eps^2: 10204.6 (Ni), 10181.9 (Cu);  b = xi^2 (f_L - f_S)/(eps^2 T):
#   109551 (Ni), -143771 (Cu);  q_f = 0.273838, q_w = 2.46e-5;  m_phi = 1, D_s = 0.
#   T = 1574 K: c_L = 0.466219, c_S = 0.399112; S = 1.0: c_inf = c_S = 0.399112.
#
# HOW IT IS SOLVED. `pf_fem_explicit`: forward Euler with lumped mass on bilinear elements,
# dx = 6.25e-3, dt = 4.75e-6 (the paper's grid and step), periodic in x and y. Seed: a disc of
# 12 dx at the solidus composition with the equilibrium interface profile.
#
# EXPECTED RESULT. A circular crystal growing at a nearly constant velocity with almost
# complete solute trapping; the solute integral is conserved to round-off. Reference
# (finite differences): equivalent radius 93.0 dx after 500 steps.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/binary/c2_check_iso_s1.0_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '-0.625 -0.625 0.625 0.625'
  divisions = '200 200'
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
[]

[Functions]
  # phi = (1 + tanh((r - R0)/(2 sqrt(2) delta)))/2, R0 = 12 dx = 0.075, delta = 9.90476e-3
  [./seed_phi]
    kind       = parsed
    expression = 0.5*(1+tanh((sqrt(x**2+y**2)-0.075)/0.0280149))
  []
  # c = c_S + (c_inf - c_S) phi: uniform at S = 1.0
  [./seed_c]
    kind       = parsed
    expression = 0.399112+0*x
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
    gradient_coefficient = 1.0
    anisotropy_strength  = 0.0
    anisotropy_mode      = 4
  []
  [./phi_bulk]
    type          = wb_phase_bulk
    variable      = phi
    concentration = c
    mobility      = 1.0
    well_a        = 10204.6      # Ni
    well_b        = 10181.9      # Cu
    driving_a     = 109551.0     # Ni: f_L - f_S > 0 below its melting point
    driving_b     = -143771.0    # Cu: above its melting point
  []
  [./dc_dt]
    type        = time_derivative
    variable    = c
    coefficient = 1.0
  []
  [./solute]
    type              = wb_solute_flux
    variable          = c
    phase             = phi
    diffusivity       = 1.0
    diffusivity_ratio = 0.0
    partition_driving = 0.273838
    partition_well    = 2.46e-5
  []
[]

[Executioner]
  type            = pf_fem_explicit
  n_steps         = 500
  dt              = 4.75e-6
  periodic_axes   = 'x y'
  snapshots_every = 100
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c'
[]
