# Part V, case K2 WITH MELT-POOL FLOW: build column of the AM Bench 2022 alloy 718 build.
#
# PROBLEM: as k2_column_h10.i (block 600 x 600 um, 200 um of plate, four powder layers of
#   40 um, 285 W, 960 mm/s, hatch 110 um, 90 degrees between layers), with the flow of the
#   liquid of k1_flow_template.i: thermocapillary shear and recoil pressure on the flat
#   build surface, buoyancy, drag of the mushy zone, heat carried by the flow. Three
#   fields are advanced together:
#
#     temperature          enthalpy heat equation with the heat carried by the flow
#     velocity, pressure   weakly compressible flow of the liquid (PINNACLE ADR 0052)
#     eta (30)             grains: melting, epitaxial regrowth, nucleation
#
# THE BLOCK STORY: the blocks of k2_column_h10.i, and from k1_flow_template.i the kernels
#   carried, u_dot, momentum, marangoni, p_dot, mass. Every flow kernel names the process
#   state: the terms vanish in elements that hold no material yet, so the build surface
#   of layers one to three, which has unbuilt mesh above it, is a free surface, and the
#   rates at its nodes use the mass of the deposited elements (PINNACLE ADR 0054). The
#   surface forces follow the surface height of the scan path. Powder is at rest through
#   the drag until it melts.
#
# LOADING: the beam only; thermal walls and bottom as in k2_column_h10.i. The velocity is
#   zero on the four cut sides and on the bottom: with flow the pool reaches the sides of
#   this small block (the liquid of two neighboring tracks merges and spreads over the
#   powder), so the sides are walls for the liquid. A block wide enough to keep the pool
#   away from its sides is the standalone deck.
#
# TIME: explicit steps of 0.5 us (the conduction deck: 2.5 us), bound by the upwind
#   diffusion at the speed of the pool; velocity and pressure take four sub-steps per
#   step (sound speed 25 m/s). 45,200 steps to 22.6 ms. Frames every 100 us.
#
# FLOW PARAMETERS are those of the single-track study k1_flow_study_h10.i (recoil
#   coefficient 0.54; d gamma / dT, expansion and drag constants ASSUMED); nothing is
#   fitted here. Grain parameters: see k2_column_h10.i.
#
# UNITS: SI (m, s, K, W, J, Pa).

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 20261125
  output_dir = ../../../../results/lpbf/in718/k2/column_flow_h10
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '0.0 0.0 0.0 600.0e-6 600.0e-6 360.0e-6'
  divisions = '60 60 36'
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
    tracks            = 5
    hatch             = 110.0e-6
    layers            = 4
    layer_thickness   = 40.0e-6
    layer_rotation    = 90.0
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
    record_file             = ../../../../results/lpbf/in718/k2/column_flow_h10/process_state.npz
    series_file             = ../../../../results/lpbf/in718/k2/column_flow_h10/temperature_series.npz
    state_series_file       = ../../../../results/lpbf/in718/k2/column_flow_h10/consolidated_series.npz
    series_interval         = 1.0e-4
    record_every            = 2000         # the record file is refreshed every 1 ms: a sign of life
  []
  [./nuclei]
    kind                = undercooling_nucleation
    variable            = eta
    temperature         = temperature
    material            = in718
    site_density        = 3.7e13     # assumed: one site per (30 um)^3
    mean_undercooling   = 5.0        # assumed
    undercooling_spread = 2.0        # assumed
    seed_radius         = 15.0e-6
    clear_above         = 200.0e-6   # the initial grains fill the plate only
    process_state       = state      # a seed stays inside consolidated material
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
    boundary_width      = 40.0e-6
    mobility_prefactor  = 0.111      # L at the liquidus, m3/(J s)
    activation_energy   = 1.5e6      # numerical: the front mobility is gone 100 K below the liquidus
    mobility_reference_temperature = 1610.15
    melting_rate        = 2.0e5
    solidification_rate = 1.0e5
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
    n_grains       = 1000
    n_order_params = 30
    int_width      = 20.0e-6
    seed           = 7
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
  end_time        = 2.26e-2          # four layers: 4 x 3.09 ms scanning + 3 x 2.5 ms dwell + 2.7 ms to freeze
  energy_every    = 0
  snapshots_every = 200             # 100 us
  staged_variables = 'velocity pressure'
  staged_substeps  = 4
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
