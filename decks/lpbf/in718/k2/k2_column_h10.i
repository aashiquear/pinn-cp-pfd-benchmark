# Part V, case K2: build column of the AM Bench 2022 alloy 718 build, coarse check.
#
# PROBLEM: a block 600 x 600 um of a leg of the AMB2022-01 bridge: 200 um of build plate
#   and four powder layers of 40 um, scanned with the build parameters (285 W, 960 mm/s,
#   hatch 110 um, 90 degrees between layers, beam D4sigma 77 um). Temperature and grain
#   structure are advanced together:
#
#     temperature   enthalpy heat equation: conduction with latent heat, beam source,
#                   surface losses, far-field exchange on the four cut sides
#     eta (30)      one order parameter per crystal orientation: melting, epitaxial
#                   regrowth with selection along the temperature gradient, nucleation
#
# THIS DECK: the coarse check: 10 um elements (129,600), 4 of the 17 layers the 3-D EBSD
#   record covers, a scanned patch of 440 x 440 um (five tracks) inside the block. At this
#   mesh the grain boundary is 40 um wide: the run shows the mechanisms and the cost, not
#   the measured grain size (median 11.5 um), which needs the 2 um standalone deck.
#
# BLOCK STORY: [Materials] in718: property tables (docs/16_lpbf_notebook.md section 2).
#   [Utilities] scan: generated raster; state: deposition, melting, freezing records, the
#   temperature of the block at the start of every layer; nuclei: nucleation sites, and
#   the removal of the initial grains above the plate. [Kernels]: storage, conduction,
#   laser, surface, sides act on the temperature; eta_dot, grains on the order parameters.
#   [ICs]: plate at the preheat temperature; Voronoi grains of about 50 um painted onto
#   30 order parameters (the nucleation service then clears everything above the plate).
#
# LOADING: the beam only. The bottom face is held at the plate temperature; the four sides
#   exchange heat with the surrounding part (film coefficient k / L, L = 1 mm).
#
# TIME: explicit steps of 2.5 us. Between two layers the scan path waits 2.5 ms (the pool
#   freezes), then `layer_reset_temperature` puts the block at the plate temperature: the
#   recoating dwell of the machine (tens of seconds) is not stepped through.
#
# MODEL PARAMETERS of the grain stage are numerical at this mesh (boundary width 4 h,
#   L(T_l) kappa = 3e-6 m2/s, P dt = D dt / 2 = 0.5 ... see the notebook, section 6) or
#   assumed (site density, undercooling); they are calibrated on the fine mesh.
#
# UNITS: SI (m, s, K, W, J).

[Global]
  backend    = torch
  use_gpu    = true
  precision  = float64
  seed       = 20261125
  output_dir = ../../../../results/lpbf/in718/k2/column_h10
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
    record_file             = ../../../../results/lpbf/in718/k2/column_h10/process_state.npz
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
    mobility_prefactor  = 1.78e4     # L(T_l) = 0.111 m3/(J s) with Q below
    activation_energy   = 1.6e5
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
  end_time        = 2.26e-2          # four layers: 4 x 3.09 ms scanning + 3 x 2.5 ms dwell + 2.7 ms to freeze
  energy_every    = 0
  snapshots_every = 400
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
