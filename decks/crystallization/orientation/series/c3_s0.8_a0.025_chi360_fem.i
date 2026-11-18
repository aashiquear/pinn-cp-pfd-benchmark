# Part IV, case C3 (FE route), parameter series of the orientation-field model:
# supersaturation S = 0.8, anisotropy of the interface free energy s0 = 0.025 (four-fold),
# orientational mobility chi_l = 360 (360 = value of (granasy 2005 spherulites)).
# Member of the series: anisotropy. Conditions otherwise those of Fig. 4 of the paper.
#
# PHYSICS. Model C of (granasy 2004 modelling), appendix A.3, with x = 0 (single-minimum
# orientational energy); phi = 1 liquid, phi = 0 crystal, c = Cu fraction, theta in [0, 1) the
# orientation in units of 2 pi / 4:
#
#     phi_t   = div[s^2 grad phi + s s' grad_perp phi] - a_w(c) g'(phi) - p'(phi) [ b(c) - h F(u) ] + zeta_phi
#     c_t     = div{ p(phi) [ grad c - c (1 - c) (q_f p'(phi) + q_w g'(phi)) grad phi ] }
#     theta_t = chi_l p(phi) div{ [1 - p(phi)] (pi/2) cos(pi u) grad(theta)/|grad(theta)| } + zeta_theta
#     s = 1 + s0 cos(4 vartheta - 2 pi theta),  F(u) = |sin(pi u)| (u < 1/2),  u = xi0 |grad theta|
#
# In the liquid the orientation fluctuates freely; in the interface it orders against the
# front; in the crystal it is frozen (chi = 0). If ordering is slow against the front
# (small chi_l), orientational defects are trapped: new grains form at the growth front.
#
# UNITS. Lengths in xi = 2.1e-4 cm, time in xi^2/D_l (the papers' dimensionless form).
#
# PARAMETERS (values and provenance: docs/14_phase_field_crystallization_plan.md, sections 7.2, 8.1).
#   Binary model as in case C2; c_inf = 0.412533; xi0 = dx.
#   ASSUMED, not given in any source read: h = 1900 (from H T / 2 = gamma_SL: the largest
#   grain-boundary energy equal to the crystal-liquid interface energy); the regularization
#   of |grad theta| (1e-3 on u); the seed; the step count; the noise amplitudes 0.0025 (phi)
#   and 0.0375 (theta) taken as standard deviations of the nodal increment per step.
#
# HOW IT IS SOLVED. `pf_fem_explicit` on the GPU: forward Euler with lumped mass on bilinear
# elements, dx = 6.25e-3, dt = 4.75e-6, periodic; theta is an angular variable (differences to
# the nearest equivalent value) advanced in 20 sub-steps per step. 2000 steps, field series
# every 100 steps.
#
# EXPECTED RESULT. See the series table of the Part IV report (orientation-field model, "What to expect when the parameters change").
# One realization of a noise-driven process: grain shapes differ with the seeds, the form does not.

[Global]
  backend    = torch
  use_gpu    = true        # explicit residual on the GPU (rocm-server); false runs on the CPU
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/orientation/series/c3_s0.8_a0.025_chi360_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '-1.5625 -1.5625 1.5625 1.5625'   # 500 dx = 3.125 xi = 6.6 um
  divisions = '500 500'
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
  # orientation: the crystal angle in units of 2 pi / k, a periodic variable in [0, 1)
  [./theta]
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
  # c = c_S + (c_inf - c_S) phi
  [./seed_c]
    kind       = parsed
    expression = 0.399112+0.013421*0.5*(1+tanh((sqrt(x**2+y**2)-0.075)/0.0280149))
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
  # no orientational order in the liquid: independent random values; the seed crystal has theta = 0
  [./theta0]
    variable = theta
    kind     = uniform_random
    minimum  = 0.0
    maximum  = 1.0
    seed     = 20261117
    center   = '0 0'
    radius   = 0.075
    inside   = 0.0
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
    anisotropy_strength  = 0.025
    anisotropy_mode      = 4
    orientation          = theta      # s = 1 + s0 cos(4 vartheta - 2 pi theta): each grain has its frame
  []
  # orientational free energy seen by the phase field: + p'(phi) h F(xi0 |grad theta|)
  [./phi_orientation]
    type               = orientation_phase_coupling
    variable           = phi
    orientation        = theta
    correlation_length = 6.25e-3     # xi0 = dx (paper)
    branch_fraction    = 0.0         # x = 0: single-minimum energy (no fixed-angle branching)
    branch_n           = 0.5
    strength           = 1900.0      # h = xi^2 H / (2 xi0 eps^2): ASSUMED (see header)
    mobility           = 1.0
    liquid_value       = 1.0
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
  # Thermal noise of the phase field in the liquid (phi = 1) only. The paper gives the
  # amplitude 0.0025 without saying what it multiplies (`unavailable`); here it is the standard
  # deviation of the nodal increment per step (`assumed`, as in the reference run):
  # A = 0.0025 / (dt x 0.3334), 0.3334 being the averaging factor of the lumped bilinear load.
  [./phi_noise]
    type         = langevin_noise
    variable     = phi
    amplitude    = 1578.7
    distribution = normal
    state_factor = phi**3*(10-15*phi+6*phi**2)
    seed         = 20261117
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
  [./dtheta_dt]
    type        = time_derivative
    variable    = theta
    coefficient = 1.0
  []
  # theta_t = chi(phi) div{ [1 - p(phi)] (pi/2) cos(pi u) grad(theta)/|grad(theta)| },  u = xi0 |grad theta|
  [./theta_flux]
    type               = orientation_grain_flux
    variable           = theta
    phase              = phi
    correlation_length = 6.25e-3
    branch_fraction    = 0.0
    branch_n           = 0.5
    mobility_liquid    = 360       # chi_l (360 in the paper)
    mobility_crystal   = 0.0         # the grain structure is frozen in the crystal
    regularization     = 1.0e-3
    liquid_value       = 1.0
  []
  # orientational noise in the liquid: per-step nodal increment of standard deviation 0.0375 p(phi)
  [./theta_noise]
    type                  = langevin_noise
    variable              = theta
    amplitude             = 23679.0     # 0.0375 / (dt x 0.3334)
    distribution          = normal
    state_factor          = phi**3*(10-15*phi+6*phi**2)
    state_factor_variable = phi
    seed                  = 20261118
  []
[]

[Executioner]
  type            = pf_fem_explicit
  n_steps         = 2000
  dt              = 4.75e-6
  periodic_axes   = 'x y'
  angular_variables = 'theta'
  subcycles       = 'theta 20'    # the orientation equation needs about 1/20 of the step (2004 review)
  snapshots_every = 100
  energy_every    = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c theta'
[]
