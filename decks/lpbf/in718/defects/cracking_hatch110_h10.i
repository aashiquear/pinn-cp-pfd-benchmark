# Part V, defect indicators: solidification-cracking indicators at the build condition.
#
# PROBLEM: the block of lof_hatch110_h10.i (600 x 600 um of alloy 718, 200 um of plate,
#   three tracks of 440 um, 285 W, 960 mm/s, D4sigma 77 um, hatch 110 um), conduction
#   model, with ONE powder layer of 40 um and a cool-down until all material is solid.
#   Two layers are not used here: the change of layer sets the block to the plate
#   temperature in one step, which freezes the pool of the first layer artificially,
#   and a passage through the mushy zone that is cut short carries no cracking
#   indicator (first pass of 2026-11-26, kept as cracking_hatch110_h10_first_void).
#
# WHAT IS READ: two `cracking_indicators` reports (PINNACLE ADR 0056) that differ in the
#   curve of the solid fraction only:
#     crack_linear   f_s linear between the solidus and liquidus of the thermal model
#     crack_scheil   Scheil relation of the Ni-Nb stand-in of alloy 718: partition
#                    coefficient 0.48, liquidus slope -10.5 K per wt pct, 5 wt pct Nb
#                    (melting point of the solvent = liquidus + 52.5 K), remaining liquid
#                    frozen at the gamma-Laves eutectic, about 1185 C = 1458 K
#   Per element that froze: the time spent between f_s 0.40 and 0.90 and between 0.90
#   and 0.99 (clyne 1981 time ratio), and the pressure drop over the mushy zone of
#   (rappaz 1999 hottearing) from the gradient and the cooling rate at freezing.
#
# ASSUMED (no measurement for this alloy at these rates): arm spacing 1 um, shrinkage
#   factor 0.03, cavitation depression 2 kPa (the value of the paper, from castings);
#   strain rate = thermal expansion 1.6e-5 /K x cooling rate. Viscosity as in the flow
#   decks. The reports are indicators, to be read relative to each other.
#
# THE BLOCK STORY: the blocks of lof_hatch110_h10.i; [Utilities] crack_linear and
#   crack_scheil are listed after state.
#
# LOADING: the beam only. Bottom at the plate temperature; far-field exchange on the sides.
#
# UNITS: SI (m, s, K, W, J).

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 20261125
  output_dir = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '0.0 0.0 0.0 600.0e-6 600.0e-6 240.0e-6'
  divisions = '60 60 24'
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
    hatch             = 110.0e-6
    layers            = 1
    layer_thickness   = 40.0e-6
    layer_rotation    = 0.0
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
    record_file             = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10/process_state.npz
  []
  [./defects]
    kind                  = defect_indicators
    scan_path             = scan
    process_state         = state
    material              = in718
    absorptivity          = 0.7724
    substrate_temperature = 353.15
    boiling_temperature   = 3186.0
    hatch                 = 110.0e-6
    layer_thickness       = 40.0e-6
    # the strip between the outer track axes, away from the track ends
    void_region           = '190.0e-6 150.0e-6 200.0e-6 410.0e-6 450.0e-6 240.0e-6'
    report_file           = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10/defects.json
  []
  [./crack_linear]
    kind                  = cracking_indicators
    process_state         = state
    material              = in718
    arm_spacing           = 1.0e-6       # assumed
    viscosity             = 5.3e-3
    shrinkage             = 0.03         # assumed
    cavitation_depression = 2.0e3        # value of the paper (castings)
    expansion             = 1.6e-5       # strain rate = expansion x cooling rate
    region                = '190.0e-6 150.0e-6 60.0e-6 410.0e-6 450.0e-6 240.0e-6'
    solid_fraction        = linear
    report_file           = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10/cracking_linear.json
    map_file              = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10/cracking_linear.npz
  []
  [./crack_scheil]
    kind                  = cracking_indicators
    process_state         = state
    material              = in718
    arm_spacing           = 1.0e-6       # assumed
    viscosity             = 5.3e-3
    shrinkage             = 0.03         # assumed
    cavitation_depression = 2.0e3        # value of the paper (castings)
    expansion             = 1.6e-5       # strain rate = expansion x cooling rate
    region                = '190.0e-6 150.0e-6 60.0e-6 410.0e-6 450.0e-6 240.0e-6'
    solid_fraction        = back_diffusion
    melting_point         = 1662.65      # liquidus + 10.5 K per wt pct x 5 wt pct Nb
    partition             = 0.48
    eutectic_temperature  = 1458.15
    report_file           = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10/cracking_scheil.json
    map_file              = ../../../../results/lpbf/in718/defects/cracking_hatch110_h10/cracking_scheil.npz
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
  end_time        = 8.0e-3           # three tracks (1.8 ms) and the time to freeze
  energy_every    = 0
  snapshots_every = 0
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
