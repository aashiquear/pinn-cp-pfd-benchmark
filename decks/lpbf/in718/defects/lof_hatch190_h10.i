# Part V, defect indicators: lack of fusion against the hatch spacing (conduction model).
#
# PROBLEM: a block 600 x 600 um of alloy 718, 200 um of plate and two powder layers of
#   40 um. Three tracks of 440 um per layer at the build parameters (285 W, 960 mm/s,
#   D4sigma 77 um), the same track positions in both layers (no rotation: the
#   unfavorable case), hatch spacing 190 um. The series lof_hatch110 / 150 / 190 / 230
#   changes the hatch only. The build of AM Bench 2022 used 110 um.
#
# WHAT IS READ: the `defect_indicators` report (PINNACLE ADR 0055): the share of the
#   deposited material between the outer track axes that never melted (voids, resolved
#   to the 10 um element), and per track the overlap index (H/W)^2 + (L/D)^2 of
#   (tang 2017 lackoffusion) from the pool width and depth of the same run. The two are
#   independent statements about lack of fusion from one run.
#
# THE BLOCK STORY: the thermal blocks of k2_column_h10.i (no grain field); [Utilities]
#   defects is listed after state.
#
# LOADING: the beam only. Bottom at the plate temperature; far-field exchange on the sides.
#
# UNITS: SI (m, s, K, W, J).

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 20261125
  output_dir = ../../../../results/lpbf/in718/defects/lof_hatch190_h10
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '0.0 0.0 0.0 600.0e-6 600.0e-6 280.0e-6'
  divisions = '60 60 28'
[]

[Variables]
  [./temperature]
    family     = fem
    components = 1
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
    track_length      = 440.0e-6
    tracks            = 3
    hatch             = 190.0e-6
    layers            = 2
    layer_thickness   = 40.0e-6
    layer_rotation    = 0.0            # the same tracks in both layers: the unfavorable case
    first_layer_angle = 90.0         # odd layers along Y, even layers along X
    serpentine        = true
    center            = '300.0e-6 300.0e-6 240.0e-6'
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
    record_file             = ../../../../results/lpbf/in718/defects/lof_hatch190_h10/process_state.npz
  []
  [./defects]
    kind                  = defect_indicators
    scan_path             = scan
    process_state         = state
    material              = in718
    absorptivity          = 0.7724
    substrate_temperature = 353.15
    boiling_temperature   = 3186.0
    hatch                 = 190.0e-6
    layer_thickness       = 40.0e-6
    # the strip between the outer track axes, away from the track ends
    void_region           = '110.0e-6 150.0e-6 200.0e-6 490.0e-6 450.0e-6 280.0e-6'
    report_file           = ../../../../results/lpbf/in718/defects/lof_hatch190_h10/defects.json
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
    band_thickness           = 10.0e-6
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
    band_thickness  = 10.0e-6
    coefficient     = 2.0e4          # k / L with k = 20 W/(m K), L = 1 mm
    far_temperature = 353.15
  []
[]

[ICs]
  [./plate]
    variable = temperature
    kind     = constant
    value    = 353.15
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
  dt              = 2.5e-6
  n_steps         = 1000000
  end_time        = 6.5e-3           # two layers of three tracks, one dwell, time to freeze
  energy_every    = 0
  snapshots_every = 0
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
