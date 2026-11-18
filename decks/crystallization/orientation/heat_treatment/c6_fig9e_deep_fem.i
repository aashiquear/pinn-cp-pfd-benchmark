# Part IV, case C6 (FE route): heat treatments before spherulitic solidification, after
# Fig. 9(d), (e) of (granasy 2005 spherulites): two thermal histories with the same final
# temperature. Conditions of its Fig. 4(f): S = 1.0, s0 = 0.1, chi_l = 54; no seed, the crystals
# nucleate from the noise; 500 x 500 as in the paper.
# This deck: (e) deep quench: 3000 steps at 1595 K, 100 steps at 1520 K, then 1574 K.
# The temperature enters through the two driving terms of the phase-field equation,
# b_i(T) = K_i (T_i/T - 1) (`derived` from the values of the model at 1574 K and 1595 K:
# K_Ni = 1.11761e6, T_Ni = 1728.3 K, K_Cu = 1.04905e6, T_Cu = 1358.3 K), given as functions of
# time. The other coefficients are kept at their 1574 K values (`assumed`).
# `assumed`: the paper quenches to 1350 K; there the driving term is three times larger and the
# explicit step of this deck is not stable, so 1520 K is used, for 100 steps.
# The phase-field noise is `calibrated` (standard deviation 0.0062 of the nodal increment per
# step, as in case C5). Paper: the durations of the stages are not printed.
# EXPECTED RESULT. Both histories end in the same growth form (space-filling spherulites); the
# deep quench leaves more crystals, because nucleation is faster at the lower temperature.
#
[Global]
  backend    = torch
  use_gpu    = true        # explicit residual on the GPU (rocm-server); false runs on the CPU
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/orientation/heat_treatment/c6_fig9e_deep_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '-1.5625 -1.5625 1.5625 1.5625'   # 500 dx
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
  # driving forces b_Ni, b_Cu against time: 1595 K for 3000 steps (t < 0.01425), 1520 K for
  # 100 steps (t < 0.014725), then 1574 K
  [./b_ni]
    kind       = parsed
    expression = 93394.0+59753.6*(t>=0.01425)-43596.6*(t>=0.014725)
  []
  [./b_cu]
    kind       = parsed
    expression = -155690.0+44080.1*(t>=0.01425)-32161.1*(t>=0.014725)
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
    driving_a_function    = b_ni
    driving_b_function    = b_cu
    bounded_interpolation = true     # the liquid is the stable phase during the heat treatment
  []
  # Thermal noise of the phase field in the liquid (phi = 1) only. The paper gives the
  # amplitude 0.0025 without saying what it multiplies (`unavailable`); here it is the standard
  # deviation of the nodal increment per step (`assumed`, as in the reference run):
  # A = 0.0062 / (dt x 0.3334), 0.3334 being the averaging factor of the lumped bilinear load.
  [./phi_noise]
    type         = langevin_noise
    variable     = phi
    amplitude    = 3915.0
    distribution = normal
    state_factor = phi**3*(10-15*phi+6*phi**2)
    seed         = 20261111
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
    seed                  = 20261118
  []
[]

[Executioner]
  type            = pf_fem_explicit
  n_steps         = 4300
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
