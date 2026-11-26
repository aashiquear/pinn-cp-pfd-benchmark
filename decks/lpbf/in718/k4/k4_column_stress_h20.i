# Part V, case K4: stress left by the process in the build column of alloy 718 (coarse).
#
# PROBLEM: the block of case K2 (600 x 600 um: 200 um of plate, four layers of 40 um,
#   scanned patch 440 x 440 um) as a crystal-plasticity body. The temperature history of
#   the column run k2_column_h10_r2.i (run that deck first) loads it through
#   `thermal_load`: thermal expansion, a slip resistance and a stiffness that fall with
#   temperature, molten material that carries (almost) no stress and solidifies stress
#   free. Material that does not exist yet (unbuilt layers, loose powder) is treated as
#   molten material: the consolidated share of the column run tells which points are
#   metal at which time (`consolidated_function`, PINNACLE ADR 0054). Powder that never
#   melts (outside the scanned patch) stays without stress to the end.
#
# GRAINS: those the column run produced (final grain map, file grains_k2_column_r2.csv
#   written by `scripts/lpbf_k3.py cube <run> <file> 0 0 0 600e-6`); every quadrature
#   point takes the orientation of the nearest node of the column mesh.
#
# MATERIAL: `dislocation_density_cp` with the set CALIBRATED in case K3 at room
#   temperature (k3_cube_n6_seed1.i: tau_c0 = 263.37 MPa, k1 = 0.052690,
#   annihilation_diameter = 1.4891 nm; elastic constants and rate constants as there).
#
# ASSUMED (no measurement of this report behind them; study deck):
#   - linear expansion 13.0e-6 /K at the plate temperature with slope 5.0e-9 /K2
#     (16.5e-6 /K at 1,050 K), stopped at the solidus;
#   - slip-resistance factor with temperature (table): 1.0 to 800 K, 0.85 at 923 K,
#     0.60 at 1,073 K and above (the shape of the yield stress of wrought alloy 718 up
#     to 1,073 K, not a fit);
#   - above 1,273 K the material carries no stress: it is treated as molten (a cut-off
#     temperature, in common use in process-stress models; the solidus is 1,528 K);
#     residual stiffness 0.01; the thermal expansion stops at the same temperature;
#   - one factor for all elastic constants, slope -3.5e-4 /K (floor 0.2).
#   WHY THE FLOOR OF 0.60: the slip law has the exponent 20. A resistance factor f
#   multiplies the slip rate by f^-20; a first version of this deck let f fall to 0.05
#   at the solidus (a factor 1e26) and did not converge in its first load step.
#
# LOADING: the temperature only. Stress-free state: the plate at 353.15 K.
#
# THREE STAGES: (1) the build, 0 to 22.6 ms, with the four side walls and the bottom held
#   (no normal displacement): the block is part of a larger leg; (2) cool-down to the
#   plate temperature, 22.6 to 24.6 ms (one frame of uniform 353.15 K appended to the
#   series by `scripts/lpbf_mechanics.py cooldown`; the thermal run stops 2.7 ms after
#   the last track); (3) relaxation with the walls x = 600 um and y = 600 um removed
#   (`release_time`), 1 ms, as in the 316L case. history.csv holds every step, so the
#   state before the release is kept.
#
# TIME: load steps of 100 us (the beam travels 96 um per step: the stress during the
#   scan is not resolved in time; the stress after freezing is set by the cooling).
#
# MESH: 20 um (30 x 30 x 18 elements, twice the element edge of the column run).
#
# NOT COMPARED WITH A MEASUREMENT: the residual-strain record of the build (NIST
#   mds2-2711) is at the scale of the bridge; this block has no part around it.
#
# UNITS: SI (m, Pa, s, K).

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ../../../../results/lpbf/in718/k4/column_stress_h20
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '0.0 0.0 0.0 600.0e-6 600.0e-6 360.0e-6'
  divisions = '30 30 18'
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
    C11                   = 250.0e9
    C12                   = 160.0e9
    C44                   = 120.0e9
    mu                    = 81.0e9
    gamma0                = 1.0e-3
    rate_exponent_n       = 20.0           # 1 / m, m = 0.05
    tau_c0                = 2.633653e+08   # calibrated (case K3)
    q_lat                 = 1.4
    alpha0                = 0.3
    burgers               = 2.54e-10
    rho0                  = 1.0e10
    k1                    = 5.268964e-02   # calibrated (case K3)
    annihilation_diameter = 1.489077e-09   # calibrated (case K3)
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'       # replaced per point by [OrientationField]
    []
  []
[]

[OrientationField]
  [./source]
    type          = point_file
    path          = grains_k2_column_r2.csv
    interpolation = nearest
  []
[]

[Functions]
  [./temperature]
    kind = field_series
    file = ../../../../results/lpbf/in718/k2/column_h10_r2/temperature_series_cooled.npz
  []
  [./built]
    kind = field_series
    file = ../../../../results/lpbf/in718/k2/column_h10_r2/consolidated_series.npz
  []
[]

[Utilities]
  [./heat]
    kind                  = thermal_load
    temperature_function  = temperature
    consolidated_function = built
    reference_temperature = 353.15
    expansion_coefficient = 13.0e-6      # assumed
    expansion_slope       = 5.0e-9       # assumed
    expansion_cutoff      = 1273.0
    resistance_table      = '353.15 1.0 800.0 1.0 923.0 0.85 1073.0 0.60'   # assumed; constant above
    rate_sensitivity      = 0.05
    stiffness_slope       = -3.5e-4      # assumed
    stiffness_floor       = 0.2
    melting_temperature   = 1273.0       # assumed cut-off: no stress above it
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
    release_time = 24.61e-3         # relaxation stage: this wall is removed
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
    release_time = 24.61e-3         # relaxation stage: this wall is removed
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
  n_steps         = 256              # 226 build steps, 20 cool-down steps, 10 relaxation steps of 100 us
  t_end           = 25.6e-3
  newton_tol      = 1.0e-7
  newton_rel_tol  = 1.0e-6
  max_iter        = 30
  max_cutbacks    = 400
  linear_method   = direct
  line_search     = backtracking
  snapshots_every = 4
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
