# Part IV, case B2 (FE route), realization 3 of the seed study (noise seed 20261120): noise-induced crystal nucleation, growth and impingement in
# two dimensions; transformed fraction against the Johnson-Mehl-Avrami-Kolmogorov law.
#
# PHYSICS. Single-field model with a Langevin noise (granasy 2019 nucleation, Eq. 16):
#
#     dphi/dt = M [ eps^2 lap(phi) - w g'(phi) - p'(phi) Delta_f ] + zeta,
#     <zeta(r,t) zeta(r',t')> = 2 k_B T M delta(r - r') delta(t - t'),
#     g = phi^2 (1 - phi)^2 / 4,   p = phi^3 (10 - 15 phi + 6 phi^2).
#
# The noise produces crystal-like fluctuations in the undercooled liquid (phi = 0); those
# beyond the critical size grow at the front velocity v until they impinge. For a constant
# nucleation rate J and a constant growth rate in d dimensions the transformed fraction is
#
#     X(t) = 1 - exp{ -[(t - t0)/tau]^p },   p = 1 + d = 3 in two dimensions
#
# (review Eq. 5-6 and section 3.1.1.3; phase-field runs of this kind give p = 2.9 to 3.0).
#
# UNITS. Scaled: lengths in delta = sqrt(2 eps^2 / w), energies in w delta^2 per unit
# thickness, time in 1/(M w): eps^2 = 0.5, w = 1, M = 1.
#
# PARAMETERS (`assumed`: a generic single-field case, not a material).
#   Delta_f / w = -0.20, a strongly driven liquid: front velocity 0.8208 delta per time unit
#   (travelling-wave reference), critical nucleus of the 2-D model W* = 0.0251 w delta^2 with
#   centre value phi(0) = 0.187 (shooting solution), smaller than the interface width.
#   Noise strength A^2 = 2 M k_B T = 0.0064 (k_B T = 0.0032 w delta^2 per unit thickness,
#   W*/k_B T = 7.8), drawn as discrete white noise (`scaling = white`: variance
#   A^2 / (cell volume x step)) and applied in the liquid only (`state_factor = 1 - p(phi)`):
#   inside the crystal the standard p makes the free energy unbounded for phi > 1, and a
#   noise acting there drives phi away without limit (measured: phi = 16 before failure).
#
# HOW IT IS SOLVED. 160 x 160 delta, 320 x 320 bilinear elements (h = delta / 2), natural
# (no-flux) boundaries, backward Euler with Newton, dt = 0.1. The field series is written
# every 5 steps for the animation.
#
# EXPECTED RESULT. After an incubation period nuclei appear at random places and times,
# grow as discs and impinge. `mass` of solver_history.csv divided by the area 25600 is the
# transformed fraction X(t) (the figure script uses the fraction of nodes with phi > 0.5 from
# the field series); the Avrami plot ln[-ln(1 - X)] against ln(t - t0) is a straight
# line for 0.01 < X < 0.95 with slope p close to 3 (the fit and its uncertainty are in the
# report; the pattern is one realization of a random process: another seed gives other
# nuclei and the same law).

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/noise/b2_jmak_2d_seed3_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 160 160'
  divisions = '320 320'
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[ICs]
  [./liquid]
    variable = phi
    kind     = constant
    value    = 0.0
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'
    [./coefficients]
      w  = 1.0
      df = -0.20
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
    gradient_coefficient = 0.5
    spatial_dim          = 2
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.5
    spatial_dim          = 2
  []
  [./thermal_noise]
    type         = langevin_noise
    variable     = phi
    amplitude    = 0.08          # A = sqrt(2 M k_B T)
    scaling      = white
    distribution = normal
    # noise in the liquid only: zeta = zeta_l [1 - p(phi)] (no fluctuations in the crystal)
    state_factor = 1-phi**3*(10-15*phi+6*phi**2)     # unquoted, no blanks
    seed         = 20261120
  []
[]

[Executioner]
  type            = pf_fem_transient
  scheme          = bdf1
  n_steps         = 350
  dt              = 0.1
  newton_tol      = 1.0e-8
  max_newton_iter = 20
  snapshots_every = 5
[]

[Postprocessors]
  [./transformed_fraction]
    type      = feature_volume_fraction
    field     = field_phi
    threshold = 0.5
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
