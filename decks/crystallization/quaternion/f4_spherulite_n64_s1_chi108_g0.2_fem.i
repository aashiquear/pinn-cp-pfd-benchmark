# Part IV, case F4 (FE route): a polycrystalline spherulite in three dimensions by growth front nucleation,
# after (pusztai 2005 threedimensions). 64 x 64 x 64 cells, supersaturation S = 1.
#
# PHYSICS. Binary ideal-solution model of the two-dimensional cases (phi = 1 liquid, phi = 0
# crystal, c = Cu fraction) in three dimensions; the orientation is a unit quaternion
# q = (q0, q1, q2, q3), q0 = cos(delta/2), q_i = c_i sin(delta/2) (rotation by delta about c):
#
#     f_ori = [1 - p(phi)] h F(u),   u = xi0 |grad q| k / pi,   F = |sin(pi u)| (u < 1/2), 1 beyond
#     phi_t = M(n') { lap phi - a_w(c) g'(phi) - p'(phi) [ b(c) - h F(u) ] } + zeta_phi
#     c_t   = div{ p(phi) [ grad c - c (1 - c) (q_f p' + q_w g') grad phi ] }
#     q_i,t = chi_l p(phi) { div J_i - q_i sum_k q_k div J_k },  J_i = [1 - p] (pi/2) cos(pi u) xi0 grad q_i / u
#     M(n') = 1 - 3 gamma + 4 gamma (n'_1^4 + n'_2^4 + n'_3^4),  n' = R(q)^T grad phi / |grad phi|
#
# 2 xi0 |grad q| is the misorientation angle across one correlation length; for a rotation
# about a fixed axis the model is the two-dimensional one with theta in units of 2 pi / k.
# The orientation is random in the liquid, orders in the interface and is frozen in the
# crystal. The cubic mobility is taken in the frame of the local crystal, so every grain
# grows along its own <100> directions. Crystal symmetries are not applied to the
# orientational energy (triclinic, as in the paper's runs).
#
# UNITS. Lengths in xi = 2.1e-4 cm, time in xi^2/D_l (as the two-dimensional Cu-Ni cases).
#
# PARAMETERS. Binary model: as case C2 (reported_table / derived, plan section 7.2).
#   reported_text (the paper): S = 1.0 for spherulites; M_q reduced by 3 to 3.5 for the spherulite.
#   assumed: h = 1900 and chi_l as in the two-dimensional cases (the paper prints mobilities in
#   other units and no H), k = 4, gamma = 0.2 (the paper prints no magnitude; 0.2 is quoted
#   for Cu-Ni by the same authors in a later proceedings paper), the seeds, the noise
#   (standard deviation 0.0025 of the nodal increment per step, liquid only), the box size
#   (paper: 300^3; here 64^3).
#
# HOW IT IS SOLVED. `pf_fem_explicit` on trilinear elements, dx = 6.25e-3, dt = 4.75e-6,
# periodic. q is a unit-norm variable (tangent projection of the rate, normalisation, q and
# -q identified) advanced in 20 sub-steps per step. 450 steps, field series every 45.
#
# EXPECTED RESULT. With the orientational mobility reduced, the front cannot order the orientation: new grains form at the front and the crystal is a compact, nearly spherical polycrystal (compare Fig. 4).

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/quaternion/f4_spherulite_n64_s1_chi108_g0.2_fem
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '-0.2 -0.2 -0.2 0.2 0.2 0.2'
  divisions = '64 64 64'
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
  [./q]
    family     = fem
    components = 4
    order      = 2
  []
[]

[Functions]
  # seeds of radius 8 dx: the product of the tanh profiles of the single seeds (unquoted, no blanks)
  [./seed_phi]
    kind       = parsed
    expression = 0.5*(1+tanh((sqrt((x-(0.00000))**2+(y-(0.00000))**2+(z-(0.00000))**2)-0.05)/0.0280149))
  []
  [./seed_c]
    kind       = parsed
    expression = 0.399112+0.000000*(0.5*(1+tanh((sqrt((x-(0.00000))**2+(y-(0.00000))**2+(z-(0.00000))**2)-0.05)/0.0280149)))
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
  # random orientations in the liquid; every seed has one orientation (random, from seed + 1)
  [./q0]
    variable = q
    kind     = random_unit_quaternion
    seed     = 20261118
    centers  = '0.00000 0.00000 0.00000'
    radius   = 0.05
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
  []
  [./phi_mobility]
    type        = cubic_kinetic_anisotropy
    variable    = phi
    strength    = 0.2
    orientation = q
  []
  [./phi_orientation]
    type               = quaternion_phase_coupling
    variable           = phi
    orientation        = q
    correlation_length = 6.25e-3
    period_index       = 4
    strength           = 1900.0
    mobility           = 1.0
  []
  [./phi_bulk]
    type          = wb_phase_bulk
    variable      = phi
    concentration = c
    mobility      = 1.0
    well_a        = 10204.6
    well_b        = 10181.9
    driving_a     = 109551.0
    driving_b     = -143771.0
  []
  [./phi_noise]
    type         = langevin_noise
    variable     = phi
    amplitude    = 2983.6       # 0.0025 / (dt x 0.1764): measured factor of the lumped trilinear load
    distribution = normal
    state_factor = phi**3*(10-15*phi+6*phi**2)
    seed         = 20261118
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
  [./dq_dt]
    type        = time_derivative
    variable    = q
    coefficient = 1.0
  []
  [./q_flux]
    type               = quaternion_grain_flux
    variable           = q
    phase              = phi
    correlation_length = 6.25e-3
    period_index       = 4
    mobility_liquid    = 108
    mobility_crystal   = 0.0
    regularization     = 1.0e-3
  []
[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = 450
  dt                  = 4.75e-6
  periodic_axes       = 'x y z'
  unit_norm_variables = 'q'
  subcycles           = 'q 20'
  snapshots_every     = 45
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c q'
[]
