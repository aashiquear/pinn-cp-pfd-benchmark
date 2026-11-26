# Part V, 316L case (secondary benchmark): stress left by one laser track, crystal
# plasticity on the representative volume of (grilli 2022 elimination), coarse mesh.
#
# PROBLEM: 200 x 160 x 108 um of polycrystalline 316L. The temperature history of
#   track_thermal_h8.i (run that deck first) loads the crystal through `thermal_load`:
#   thermal expansion, a slip resistance and a stiffness that fall with temperature,
#   molten material that carries (almost) no stress and solidifies stress free. The run
#   follows the scan (457 us) and the cooling to 8 ms; what remains is the residual
#   stress of the track.
#
# FROM THE PAPER (its Table 2, `reported_table`): slip law gamma_dot0 = 1e-7 /s, m = 0.1;
#   hardening h0 = 3,839 MPa, exponent 2.5, initial and saturation resistance 213 and
#   302 MPa; C11, C12, C44 = 204.6, 137.7, 126.2 GPa at 303 K; resistance factor
#   k_A + k_B exp(-k_C (T - T0)) with 0.53, 0.47, 0.008; melting 1648.15 K; residual
#   stiffness 0.01; volumetric expansion 44.73e-6 /K with slope 0.01011e-6 /K2 (a third
#   of each is the linear value used here); boundary conditions: no normal displacement
#   on x = 0, x = 200 um, y = 0, y = 160 um and on the bottom; top free; mesh 25 x 20 x 14.
#
# THE HARDENING LAW IS A SUBSTITUTE. The paper uses a phenomenological law,
#   d tau_c = sum h0 (1 - tau_c / tau_sat)^a |d gamma|. PINNACLE's material of that form
#   (`kalidindi`) has no batched field evaluation, so this deck uses the dislocation-density
#   law (`dislocation_density_cp`, backward-Euler update, consistent tangent) with
#   parameters derived to give the paper's initial resistance (213 MPa), saturation
#   resistance (302 MPa) and initial hardening rate (h0 (1 - 213/302)^2.5 = 181 MPa):
#   tau_c = tau_c0 + alpha0 mu b sqrt(sum q rho), d rho = (k1 sqrt(rho) - k2 rho) |d gamma| / b.
#   Derivation: both laws harden system a by (rate) x sum_b q_ab |d gamma_b|; the rate of
#   this law at rho0 is alpha0 mu k1 / (2 sqrt(S)), S = 1 + 11 q_lat, which gives
#   k1 = 0.0652; the saturation density (k1 / k2)^2 gives k2 = 2 y_c. (A first version
#   of this deck had k1 = 4.0e-3, y_c = 5.2e-10: the same saturation, but a hardening
#   rate lower by the factor S = 16.4; its results are superseded.)
#   The flow rule (power law, gamma_dot0, m) is the paper's.
#
# NOT FROM THE PAPER: the temperature field (conduction substitute); the grains (844
#   Voronoi grains of 16 um, the mean size the paper states, with uniform random
#   orientations; the paper mapped phase-field grains and EBSD orientations); the latent
#   hardening ratio 1.4 (assumed; the paper does not print it); one stiffness factor for
#   all elastic constants, slope -3.93e-4 /K (the mean relative slope of the three
#   constants of the paper); steps of 45 us (the paper: 1 us).
#
# LOADING: the temperature only (no displacement or traction is prescribed).
#
# TWO STAGES, as in the paper: (1) scan and cooling to 8.01 ms with the four side walls
#   and the bottom fixed; (2) relaxation at the final temperature with the walls x = 200 um
#   and y = 160 um removed (`release_time`), so that the stress which remains is carried
#   by the volume and not by the walls. The paper does not print the length of its
#   relaxation stage; 22 steps of 45 us (0.99 ms) are assumed here. history.csv holds
#   the mean stress of every step, so the state before the release (step 178) is kept.
#
# UNITS: SI (m, Pa, s, K).

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ../../../results/lpbf/ss316l/mechanics_h8_seed2
[]

[Mesh]
  type      = box_mesh
  domain    = '0.0 0.0 0.0 200.0e-6 160.0e-6 108.0e-6'
  divisions = '25 20 14'
[]

[Variables]
  [./u]
    family     = fem
    components = 3
    order      = 2
  []
[]

[Materials]
  [./cp]
    kind                  = dislocation_density_cp
    C11                   = 204.6e9
    C12                   = 137.7e9
    C44                   = 126.2e9
    mu                    = 75.0e9        # assumed (shear modulus of the Taylor term)
    gamma0                = 1.0e-7
    rate_exponent_n       = 10.0          # 1 / m, m = 0.1
    tau_c0                = 210.7e6       # derived: 213 MPa less the Taylor term at rho0
    q_lat                 = 1.4           # assumed
    alpha0                = 0.3           # assumed
    burgers               = 2.54e-10
    rho0                  = 1.0e10        # 1/m2, chosen
    k1                    = 6.5155e-02    # derived: initial hardening rate of the paper's law (181 MPa)
    annihilation_diameter = 8.2581e-09    # derived: saturation resistance 302 MPa
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'       # replaced per point by [OrientationField]
    []
  []
[]

[OrientationField]
  [./source]
    type          = point_file
    path          = grains_seed2.csv
    interpolation = nearest
  []
[]

[Functions]
  [./temperature]
    kind = field_series
    file = ../../../results/lpbf/ss316l/thermal_h8/temperature_series.npz
  []
[]

[Utilities]
  [./heat]
    kind                  = thermal_load
    temperature_function  = temperature
    reference_temperature = 303.0
    expansion_coefficient = 14.91e-6     # 44.73e-6 / 3
    expansion_slope       = 3.37e-9      # 0.01011e-6 / 3
    resistance_ka         = 0.53
    resistance_kb         = 0.47
    resistance_kc         = 0.008
    rate_sensitivity      = 0.1
    stiffness_slope       = -3.93e-4
    melting_temperature   = 1648.15
    residual_stiffness    = 0.01
  []
[]

[Kernels]
  [./sd_cp]
    type          = stress_divergence_cp
    variable      = u
    material      = cp
    dim           = 3
    finite_strain = true
    kinematics    = total_lagrangian
    state_update  = implicit
    thermal_load  = heat
  []
[]

[BCs]
  [./xmin]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./xmax]
    variable     = u
    kind         = dirichlet
    node_set     = xmax
    components   = [0]
    values       = [0.0]
    release_time = 8.01e-3          # relaxation stage: this wall is removed
  []
  [./ymin]
    variable   = u
    kind       = dirichlet
    node_set   = ymin
    components = [1]
    values     = [0.0]
  []
  [./ymax]
    variable     = u
    kind         = dirichlet
    node_set     = ymax
    components   = [1]
    values       = [0.0]
    release_time = 8.01e-3          # relaxation stage: this wall is removed
  []
  [./zmin]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
[]

[Executioner]
  type            = cp_quasi_static
  n_steps         = 200              # 178 temperature frames (45 us each), then 22 relaxation steps
  t_end           = 9.0e-3
  newton_tol      = 1.0e-7           # N; forces are of order 1e-2 N
  newton_rel_tol  = 1.0e-6
  max_iter        = 30
  max_cutbacks    = 200
  linear_method   = direct
  line_search     = backtracking
  snapshots_every = 4
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
