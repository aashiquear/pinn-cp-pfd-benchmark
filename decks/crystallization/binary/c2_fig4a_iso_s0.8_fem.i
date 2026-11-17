# Part IV, case C2 (FE route): single-crystal growth form of (granasy 2005 spherulites),
# Fig. 4(a): supersaturation S = 0.8, interface-energy anisotropy s0 = 0.0 (four-fold),
# isotropic mobility; binary Cu-Ni ideal solution at 1574 K, 500 x 500 grid (6.6 um).
#
# PHYSICS. Binary phase-field model, Eq. 2.3-2.4 of the paper, with a uniform orientation
# (one crystal); phi = 1 liquid, phi = 0 crystal, c = Cu fraction:
#
#     phi_t = m_phi { div[s^2 grad phi + s s' grad_perp phi] - a_w(c) g'(phi) - p'(phi) b(c) } + zeta
#     c_t   = div{ p(phi) [ grad c - c (1 - c) (q_f p'(phi) + q_w g'(phi)) grad phi ] }
#     s = 1 + s0 cos(4 vartheta),  vartheta = atan2(phi_y, phi_x)
#
# UNITS. Lengths in xi = 2.1e-4 cm, time in xi^2/D_l (the paper's dimensionless form).
#
# PARAMETERS (values and provenance: docs/14_phase_field_crystallization_plan.md, section 7.2).
#   a_w: 10204.6 (Ni), 10181.9 (Cu); b: 109551 (Ni), -143771 (Cu); q_f = 0.273838,
#   q_w = 2.46e-5; m_phi = 1; D_s = 0. c_L = 0.466219, c_S = 0.399112,
#   c_inf = c_L - S (c_L - c_S) = 0.412533.
#   Not given in the paper and assumed here: the seed (disc of 12 dx at the solidus composition
#   with the equilibrium interface profile), the number of steps, the meaning of the noise
#   amplitude. The composition noise of the paper (0.00125, as a random flux) is not applied.
#
# HOW IT IS SOLVED. `pf_fem_explicit`: forward Euler with lumped mass on bilinear elements,
# dx = 6.25e-3, dt = 4.75e-6 (the paper's grid and step), periodic in x and y, 2400 steps.
# The field series is written every 120 steps.
#
# EXPECTED RESULT. A seaweed-like form with radial liquid channels: without anisotropy the front is morphologically unstable (Mullins-Sekerka) and solute is rejected into the channels
# (compare Fig. 4(a) of the paper). The pattern is one realization of a noise-driven
# process: details differ with the seed, the growth form does not. The solute integral
# (`c_mean`) is constant to round-off.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/binary/c2_fig4a_iso_s0.8_fem
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
[]

[Executioner]
  type            = pf_fem_explicit
  n_steps         = 2400
  dt              = 4.75e-6
  periodic_axes   = 'x y'
  snapshots_every = 120
  energy_every    = 120
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c'
[]
