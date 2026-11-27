# Part V, case K1 with a free surface: the seven NIST single tracks at 10 um (study deck).
#
# As k1_flow_study_h10.i (same heat-source absorptivity, same flow constants; nothing
# refitted), with the template k1_free_surface_template.i: 50 um of gas above the plate
# and the surface free to move. The penetration depth of the conduction fit is not
# used: the beam deposits its heat on the surface.
# New keys: gas_layer_um, surface_tension (N/m at the liquidus, assumed), phase_mobility
# (m/s, numerical), gas_damping (1/s, numerical).

[Study]
  template        = k1_free_surface_template.i
  workbook        = /mnt/extended/projects/ambench-data/nist/DATA_AMB2022_03/mds2-2718/AMB2022-718-SH1-MeltPool_Cross-Section_Measurement_Results.xlsx
  targets         = ../../../../data/lpbf/in718/k1_single_tracks.csv
  case_dir        = cases
  out_dir         = ../../../../results/lpbf/in718/k1
  tag             = free_surface_h10
  mesh_um         = 10.0
  backend         = torch
  use_gpu         = true
  absorptivity    = 0.7724
  penetration_um  = 45.01
  liquid_conductivity_factor = 1.0
  fit_cases       = 'case0 case1.2'
  end_time        = 1.2e-3
  sound_speed     = 25.0
  flow_speed_bound = 4.0
  flow_upwind     = 1.0
  surface_tension_gradient = -1.1e-4
  recoil_coefficient = 0.54
  frame_interval  = 3.0e-5
  gas_layer_um    = 50.0
  surface_tension = 1.8
  phase_mobility  = 3.0
  gas_damping     = 1.0e6
[]

[Cases]
  [./case0]
    line = 0
  []
  [./case1.1]
    line = 1.1
  []
  [./case1.2]
    line = 1.2
  []
  [./case2.1]
    line = 2.1
  []
  [./case2.2]
    line = 2.2
  []
  [./case3.1]
    line = 3.1
  []
  [./case3.2]
    line = 3.2
  []
[]
