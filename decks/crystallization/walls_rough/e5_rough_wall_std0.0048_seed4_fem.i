# Part IV, case E5 (FE route): heterogeneous nucleation on a rough wall, after Fig. 3(a) of
# (granasy 2006 review). A liquid layer of the binary Cu-Ni model with the orientation field
# (conditions of Fig. 4(f) of (granasy 2005 spherulites): S = 1.0, s0 = 0.1, chi_l = 54) lies on a
# sawtooth wall: grooves 30 cells deep, 80 cells apart. Every boundary is a no-flux boundary,
# which for the phase field is a wall with a contact angle of 90 degrees (as in the paper).
# No seed: crystals form from the noise. In a groove a nucleus touches two walls and needs a
# fraction of the homogeneous nucleus, so with a noise level just below the one that nucleates in
# the bulk the crystals appear in the grooves first.
# MESH. Written by `python3 scripts/cryst_rough_wall_mesh.py` (not shipped: 5 MB); mapped grid
# of 400 x 200 bilinear elements, 6.25e-3 wide.
# assumed: the wall shape and size, the layer thickness, the noise: standard deviation 0.0048
# of the nodal increment per step (0.77 of the level used for nucleation in the bulk, case C5).
# correlated. The paper prints 0.0025 without saying what it multiplies; its time constant is 0.0106.
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
#   Binary model as in case C2; c_inf = 0.399112; xi0 = dx.
#   ASSUMED, not given in any source read: h = 1900 (from H T / 2 = gamma_SL: the largest
#   grain-boundary energy equal to the crystal-liquid interface energy); the regularization
#   of |grad theta| (1e-3 on u); the seed; the step count; the noise amplitudes 0.0025 (phi)
#   and 0.0375 (theta) taken as standard deviations of the nodal increment per step.
#
# HOW IT IS SOLVED. `pf_fem_explicit` on the GPU: forward Euler with lumped mass on bilinear
# elements, dx = 6.25e-3, dt = 4.75e-6, no-flux boundaries; theta is an angular variable (differences to
# the nearest equivalent value) advanced in 20 sub-steps per step. 3000 steps, field series
# every 150 steps.
#
# EXPECTED RESULT. The first crystals sit in the grooves of the wall; they grow out as spherulites and cover the wall before any crystal forms in the bulk of the liquid.
# Realization 4 of the noise (seeds of the two noise terms shifted by 4).
# One realization of a noise-driven process: grain shapes differ with the seeds, the form does not.

[Global]
  backend    = torch
  use_gpu    = true        # explicit residual on the GPU (rocm-server); false runs on the CPU
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/walls_rough/e5_rough_wall_std0.0048_seed4_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type   = file
  path   = decks/crystallization/walls_rough/e5_rough_wall.msh
  format = gmsh
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

[ICs]
  # homogeneous liquid of the solidus composition (S = 1), no orientational order
  [./phi0]
    variable = phi
    kind     = constant
    value    = 1.0
  []
  [./c0]
    variable = c
    kind     = constant
    value    = 0.399112
  []
  [./theta0]
    variable = theta
    kind     = uniform_random
    minimum  = 0.0
    maximum  = 1.0
    seed     = 20261118
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
    anisotropy_strength  = 0.1
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
  # A = 0.0048 / (dt x 0.3334), 0.3334 being the averaging factor of the lumped bilinear load.
  [./phi_noise]
    type         = langevin_noise
    variable     = phi
    amplitude    = 3031.0
    distribution = normal
    state_factor = phi**3*(10-15*phi+6*phi**2)
    seed         = 20261115
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
    mobility_liquid    = 54        # chi_l (360 in the paper; 54 = reduced by the factor 0.15)
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
    seed                  = 20261122
  []
[]

[Executioner]
  type            = pf_fem_explicit
  n_steps         = 3000
  dt              = 4.75e-6
  angular_variables = 'theta'
  subcycles       = 'theta 20'    # the orientation equation needs about 1/20 of the step (2004 review)
  snapshots_every = 150
  energy_every    = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c theta'
[]
