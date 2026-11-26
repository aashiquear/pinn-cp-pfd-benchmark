# Part V, case K2: build column of the AM Bench 2022 alloy 718 build, STANDALONE deck
# (5.0 um, float64).
#
# 650 x 650 um, 200 um of plate and the 17 layers (0 to 680 um) the 3-D EBSD record of
# (schwalbach 2024 sectioningdata) covers; 5.0 um elements: 130 x 130 x 176 = 2.97 million
# elements, 30 order parameters, float64.
# SIZE: the explicit path keeps about 27 kB per element in float64 (13 kB in float32) at
#   30 order parameters and 8 Gauss points: about 80 GB in all. Split it over several
#   GPUs with `devices` in [Global] (NVIDIA CUDA or AMD ROCm, ADR 0050), e.g.
#       devices = 'cuda:0 cuda:1 cuda:2 cuda:3 cuda:4 cuda:5 cuda:6 cuda:7'
#   This deck has NOT been run by the authors of the benchmark (one 64 GB GPU); the
#   10 um deck k2_column_h10.i is the checked configuration.
# STEPS: dt = 6.386e-07 s (0.8 of the thermal limit), about 161888 steps for 17 layers.
# GRAIN MODEL at this mesh: boundary width 20 um (4 elements); L(T_l) kappa = 2.94e-06 m2/s;
#   melting rate 7.8e+05 /s, regrowth rate 4.0e+05 /s (front speed 2 sqrt(D L kappa) = 2.17 m/s).
#   Site density and undercooling are the assumed values of the coarse check: they are
#   the parameters to calibrate against the measured grain size (median 11.5 um).
#
# One difference to the checked 10 um deck: the nucleation block names the process state, so
# that a nucleus at the build surface does not reach into powder or empty space (the
# 10 um runs of the report were made before that option existed).
# Everything else (material tables, kernels, loading, units SI) is as in k2_column_h10.i,
# whose header explains the blocks. Run:
#     python3 scripts/lpbf_column.py decks/lpbf/in718/k2/<this deck> run
# (PINNACLE on PYTHONPATH), then `grains` and `figure`; or `pinnacle run <this deck>`.

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 20261125
  output_dir = ../../../../results/lpbf/in718/k2/k2_column_full_h5
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '0.0 0.0 0.0 650.0e-6 650.0e-6 880.0e-6'
  divisions = '130 130 176'
[]

[Variables]
  [./temperature]
    family     = fem
    components = 1
    order      = 2
  []
  [./eta]
    family     = fem
    components = 30
    order      = 2
  []
[]

[Materials]
  [./in718]
    type     = phase_change_thermal
    solidus  = 1528.15
    liquidus = 1610.15
    latent_heat = 1.6122e9
    heat_capacity_table = '298.0 3.4808e6 400.0 3.6609e6 500.0 3.8329e6 600.0 4.0049e6 700.0 4.1769e6 800.0 4.3489e6 1100.0 5.2007e6 1200.0 5.2007e6 1300.0 5.2007e6 1400.0 5.1925e6 1528.15 5.3972e6 1610.15 6.3473e6'
    conductivity_solid_table = '298.0 9.94 400.0 11.59 500.0 13.24 600.0 14.91 700.0 16.61 800.0 18.34 1100.0 22.72 1200.0 23.61 1300.0 24.47 1400.0 25.32'
    conductivity_liquid_table = '1610.15 26.80 1673.15 27.70 1773.15 29.10 1873.15 30.45 1973.15 31.80 2073.15 33.13'
    conductivity_powder = 0.3        # assumed: loose powder bed in argon
    reference_temperature = 298.0
  []
[]

[Utilities]
  [./scan]
    kind              = scan_path
    power             = 285.0
    speed             = 0.96
    radius            = 38.5e-6      # D4sigma 77 um / 2
    track_length      = 490.0e-6
    tracks            = 5
    hatch             = 110.0e-6
    layers            = 17
    layer_thickness   = 40.0e-6
    layer_rotation    = 90.0
    first_layer_angle = 90.0         # odd layers along Y, even layers along X
    serpentine        = true
    center            = '325.0e-6 325.0e-6 240.0e-6'
    initial_surface   = 200.0e-6     # top of the build plate
    track_dwell       = 2.0e-4
    layer_dwell       = 2.5e-3
  []
  [./state]
    kind                    = process_state
    variable                = temperature
    material                = in718
    scan_path               = scan
    substrate_top           = 200.0e-6
    layer_reset_temperature = 353.15
    record_file             = ../../../../results/lpbf/in718/k2/k2_column_full_h5/process_state.npz
  []
  [./nuclei]
    kind                = undercooling_nucleation
    variable            = eta
    temperature         = temperature
    material            = in718
    site_density        = 3.7e13     # assumed: one site per (30 um)^3
    mean_undercooling   = 5.0        # assumed
    undercooling_spread = 2.0        # assumed
    seed_radius         = 7.5e-6
    clear_above         = 200.0e-6   # the initial grains fill the plate only
    process_state       = state      # a seed takes consolidated material only (PINNACLE after 037245e)
  []
[]

[Kernels]
  [./storage]
    type          = enthalpy_time_derivative
    variable      = temperature
    material      = in718
    process_state = state
  []
  [./conduction]
    type          = phase_change_conduction
    variable      = temperature
    material      = in718
    process_state = state
  []
  [./laser]
    type              = beam_heat_source
    variable          = temperature
    scan_path         = scan
    absorptivity      = 0.7724       # K1 calibration (bare plate; kept for the powder layers)
    penetration_depth = 45.01e-6
  []
  [./surface]
    type                     = surface_heat_loss
    variable                 = temperature
    scan_path                = scan
    band_thickness           = 5.0e-6
    ambient_temperature      = 353.15
    convection_coefficient   = 10.0
    emissivity               = 0.4
    boiling_temperature      = 3186.0
    vaporization_latent_heat = 6.4e6
    molar_mass               = 0.0588
    recondensation_fraction  = 0.18
  []
  [./sides]
    type            = far_field_heat_loss
    variable        = temperature
    faces           = 'xmin xmax ymin ymax'
    band_thickness  = 5.0e-6
    coefficient     = 2.0e4          # k / L with k = 20 W/(m K), L = 1 mm
    far_temperature = 353.15
  []
  [./eta_dot]
    type     = time_derivative
    variable = eta
  []
  [./grains]
    type                = thermal_grain_growth
    variable            = eta
    temperature         = temperature
    material            = in718
    boundary_energy     = 0.9
    boundary_width      = 20.0e-6
    mobility_prefactor  = 0.218      # L at the liquidus, m3/(J s)
    activation_energy   = 1.5e6      # numerical: the front mobility is gone 100 K below the liquidus
    mobility_reference_temperature = 1610.15
    melting_rate        = 7.83e+05
    solidification_rate = 4.00e+05
    front_range         = 2.0
    selection_exponent  = 4.0
    orientation_seed    = 20261125
    process_state       = state      # grains grow into consolidated material only
  []
[]

[ICs]
  [./plate]
    variable = temperature
    kind     = constant
    value    = 353.15
  []
  [./grains0]
    variable       = eta
    kind           = polycrystal_coloring
    n_grains       = 2974
    n_order_params = 30
    int_width      = 10.0e-6
    seed           = 7
  []
[]

[BCs]
  [./plate_bottom]
    variable   = temperature
    kind       = dirichlet
    node_set   = zmin
    components = [0]
    values     = [353.15]
  []
[]

[Executioner]
  type            = pf_fem_explicit
  dt              = 6.386e-07
  n_steps         = 1000000
  end_time        = 1.0339e-01
  energy_every    = 0
  snapshots_every = 3131
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
