# Part V, 316L case (secondary benchmark): thermal field of the single track of
# (grilli 2022 elimination), coarse mesh.
#
# PROBLEM: the representative volume of the paper, 200 x 160 x 108 um of 316L, scanned by
#   one laser track along x (150 W, 700 mm/s) across its top face. The paper took the
#   temperature from a powder-scale thermal-fluid simulation that is not public; this
#   deck computes a substitute by conduction with latent heat. Its temperature history is
#   written every 45 us (the sampling of the paper) for the mechanics deck
#   track_mechanics_h8.i.
#
# WHAT IS FROM THE PAPER AND WHAT IS NOT:
#   reported   domain, power, speed, melting temperature 1648.15 K, mesh (25 x 20 x 14,
#              7,000 elements: the paper's 8 um mesh)
#   assumed    beam radius 35 um, absorptivity 0.5, penetration 30 um (the paper cites its
#              thermal-fluid papers for the beam; chosen so that the pool passes through
#              the volume with a margin, as the paper states); 316L thermal properties
#              (typical handbook values, not yet sourced: liquidus, heat capacity,
#              conductivity, latent heat, density); no powder layer (solid block)
#   Every assumed value is listed in docs/16_lpbf_notebook.md, section 7.
#
# LOADING: the beam enters at x = -60 um and leaves at x = 260 um (457 us). The four sides
#   and the bottom exchange heat with the surrounding part (film coefficient k / L,
#   L = 200 um), so the volume returns toward the ambient 298 K.
#
# TIME: explicit steps of 1.5 us to 8 ms.
#
# UNITS: SI.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ../../../results/lpbf/ss316l/thermal_h8
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '0.0 0.0 0.0 200.0e-6 160.0e-6 108.0e-6'
  divisions = '25 20 14'
[]

[Variables]
  [./temperature]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Materials]
  [./ss316l]
    type     = phase_change_thermal
    solidus  = 1648.15           # melting temperature of the paper (its Table 2)
    liquidus = 1723.15           # assumed
    latent_heat = 2.13e9         # assumed: 7900 kg/m3 x 270 kJ/kg
    heat_capacity_table = '298.0 3.95e6 1600.0 5.53e6 1723.15 6.32e6'   # assumed
    conductivity_solid_table = '298.0 14.0 1600.0 30.0'                 # assumed
    conductivity_liquid = 30.0                                          # assumed
    reference_temperature = 298.0
  []
[]

[Utilities]
  [./scan]
    kind         = scan_path
    power        = 150.0
    speed        = 0.7
    radius       = 35.0e-6
    track_length = 320.0e-6
    tracks       = 1
    center       = '100.0e-6 80.0e-6 108.0e-6'
  []
  [./state]
    kind            = process_state
    variable        = temperature
    material        = ss316l
    scan_path       = scan
    melting_temperature = 1648.15
    record_file     = ../../../results/lpbf/ss316l/thermal_h8/process_state.npz
    series_file     = ../../../results/lpbf/ss316l/thermal_h8/temperature_series.npz
    series_interval = 45.0e-6
  []
[]

[Kernels]
  [./storage]
    type          = enthalpy_time_derivative
    variable      = temperature
    material      = ss316l
    process_state = state
  []
  [./conduction]
    type          = phase_change_conduction
    variable      = temperature
    material      = ss316l
    process_state = state
  []
  [./laser]
    type              = beam_heat_source
    variable          = temperature
    scan_path         = scan
    absorptivity      = 0.5
    penetration_depth = 30.0e-6
  []
  [./surface]
    type                   = surface_heat_loss
    variable               = temperature
    scan_path              = scan
    band_thickness         = 8.0e-6
    ambient_temperature    = 298.0
    convection_coefficient = 10.0
    emissivity             = 0.4
  []
  [./surroundings]
    type            = far_field_heat_loss
    variable        = temperature
    faces           = 'xmin xmax ymin ymax zmin'
    band_thickness  = 8.0e-6
    coefficient     = 1.0e5
    far_temperature = 298.0
  []
[]

[ICs]
  [./ambient]
    variable = temperature
    kind     = constant
    value    = 298.0
  []
[]

[Executioner]
  type         = pf_fem_explicit
  dt           = 1.5e-6
  n_steps      = 1000000
  end_time     = 8.0e-3
  energy_every = 0
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
