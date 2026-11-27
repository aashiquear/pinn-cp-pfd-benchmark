# Part V, defect indicators: gas pores as tracers in the melt pool (flow model).
#
# PROBLEM: a block 600 x 600 um of alloy 718, 200 um of plate and ONE powder layer of
#   40 um; three tracks of 440 um at the build parameters (285 W, 960 mm/s, D4sigma
#   77 um, hatch 110 um), with the flow of the liquid under a flat surface (the blocks
#   of k2_column_flow_h10.i without the grain field).
#
# WHAT IS READ: the `pore_tracers` report (PINNACLE ADR 0059). When a powder element
#   first melts it releases pores whose volume is 1 % of the element on average
#   (ASSUMED gas content of the layer), 10 um in diameter with a log-normal spread of
#   0.3 (ASSUMED). A pore moves with the liquid and rises by buoyancy, is trapped where
#   the liquid freezes and lost when it reaches the surface. Read: how many pores are
#   released, lost at the surface and trapped, where the trapped ones sit, and the gas
#   porosity of the melted material. The share that survives is the result; the amount
#   released is the input.
#
# DWELL: at least 1 ms between tracks; the scan is held until the material is solid
#   (`hold_scan_until_solid`), and the run ends when all is solid.
#
# LOADING: the beam only. Thermal walls and bottom, velocity walls as in
#   k2_column_flow_h10.i.
#
# UNITS: SI (m, s, K, W, J, Pa).

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 20261125
  output_dir = ../../../../results/lpbf/in718/defects/pores_flow_h10
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
  [./velocity]
    family     = fem
    components = 3
    order      = 2
  []
  [./pressure]
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
    layer_rotation    = 90.0
    first_layer_angle = 90.0         # odd layers along Y, even layers along X
    serpentine        = true
    center            = '300.0e-6 300.0e-6 240.0e-6'
    initial_surface   = 200.0e-6     # top of the build plate
    track_dwell       = 1.0e-3
    layer_dwell       = 3.0e-3
    recoat            = dwell_end
  []
  [./state]
    kind                    = process_state
    variable                = temperature
    material                = in718
    scan_path               = scan
    substrate_top           = 200.0e-6
    layer_reset_temperature = 353.15
    record_file             = ../../../../results/lpbf/in718/defects/pores_flow_h10/process_state.npz
    series_file             = ../../../../results/lpbf/in718/defects/pores_flow_h10/temperature_series.npz
    state_series_file       = ../../../../results/lpbf/in718/defects/pores_flow_h10/consolidated_series.npz
    series_interval         = 1.0e-4
    hold_scan_until_solid   = true
    stop_when_solid         = true
  []
  [./pores]
    kind                 = pore_tracers
    process_state        = state
    scan_path            = scan
    velocity             = velocity
    powder_gas_fraction  = 0.01         # ASSUMED: gas held by the powder layer, by volume
    pore_diameter        = 10.0e-6      # ASSUMED median; log-normal spread 0.3
    pore_diameter_spread = 0.3
    density              = 7578.0
    viscosity            = 5.3e-3
    seed                 = 20261126
    region               = '135.0e-6 80.0e-6 100.0e-6 465.0e-6 520.0e-6 240.0e-6'
    report_file          = ../../../../results/lpbf/in718/defects/pores_flow_h10/pores.json
    pore_file            = ../../../../results/lpbf/in718/defects/pores_flow_h10/pores.npz
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
  [./carried]
    type          = flow_heat_advection
    variable      = temperature
    velocity      = velocity
    material      = in718
    upwind        = 1
    mesh_size     = 10.0e-6
    reference_temperature = 1610.15
    process_state = state
  []
  [./u_dot]
    type     = time_derivative
    variable = velocity
  []
  [./momentum]
    type                  = flow_momentum
    variable              = velocity
    pressure              = pressure
    temperature           = temperature
    material              = in718
    density               = 7578.0        # kg/m3, liquid
    viscosity             = 5.3e-3        # Pa s
    gravity               = '0.0 0.0 -9.81'
    expansion_coefficient = 5.0e-5        # 1/K, assumed
    reference_temperature = 1610.15
    drag_coefficient      = 1.0e8         # 1/s, assumed
    upwind                = 1
    mesh_size             = 10.0e-6
    process_state         = state
  []
  [./marangoni]
    type                     = thermocapillary_shear
    variable                 = velocity
    temperature              = temperature
    material                 = in718
    density                  = 7578.0
    surface_tension_gradient = -1.1e-4    # N/(m K), assumed
    band_thickness           = 10.0e-6
    scan_path                = scan
    recoil_coefficient       = 0.54
    boiling_temperature      = 3186.0
    vaporization_latent_heat = 6.4e6
    molar_mass               = 0.0588
  []
  [./p_dot]
    type     = time_derivative
    variable = pressure
  []
  [./mass]
    type                 = flow_continuity
    variable             = pressure
    velocity             = velocity
    density              = 7578.0
    sound_speed          = 25
    pressure_diffusivity = 1.25e-5
    process_state        = state
  []
[]

[ICs]
  [./plate]
    variable = temperature
    kind     = constant
    value    = 353.15
  []
  [./rest]
    variable = velocity
    kind     = constant
    value    = 0.0
  []
  [./p0]
    variable = pressure
    kind     = constant
    value    = 0.0
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
  [./rest_xmin]
    variable   = velocity
    kind       = dirichlet
    node_set   = xmin
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
  [./rest_xmax]
    variable   = velocity
    kind       = dirichlet
    node_set   = xmax
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
  [./rest_ymin]
    variable   = velocity
    kind       = dirichlet
    node_set   = ymin
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
  [./rest_ymax]
    variable   = velocity
    kind       = dirichlet
    node_set   = ymax
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
  [./rest_zmin]
    variable   = velocity
    kind       = dirichlet
    node_set   = zmin
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
[]

[Executioner]
  type            = pf_fem_explicit
  dt              = 5.0e-7
  n_steps         = 1000000
  end_time        = 2.0e-2           # upper limit; the run ends when the path is over and all is solid
  energy_every    = 0
  snapshots_every = 200             # 100 us
  staged_variables = 'velocity pressure'
  staged_substeps  = 4
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
