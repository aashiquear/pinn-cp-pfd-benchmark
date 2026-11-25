# Part V, case K1 with melt-pool flow: seven bare-plate tracks of NIST AMB2022-03 at 10 um.
#
# Study deck of scripts/lpbf_k1.py (targets | decks | run | collect | figure). It fills
# k1_flow_template.i: the conduction model of k1_study_h10.i with the flow of the liquid
# (PINNACLE ADR 0052). Heat source: the calibrated values of the conduction study
# (absorptivity 0.7724, penetration 45.01 um), NOT refitted with flow, so that the change
# of depth and width is the effect of the flow alone.
#
# Flow settings: artificial sound speed 25 m/s; the step of the heat equation allows for
# the upwind diffusion of a flow of up to 4 m/s (a run whose flow is faster stops with a
# non-finite field and must be rerun with a larger bound); first-order upwind.
# d gamma / dT = -1.1e-4 N/(m K) (assumed: the value in common use for alloy 718);
# recoil coefficient 0.54.
# A frame of all fields is kept every 30 us (animation stack).

[Study]
  template        = k1_flow_template.i
  workbook        = /mnt/extended/projects/ambench-data/nist/DATA_AMB2022_03/mds2-2718/AMB2022-718-SH1-MeltPool_Cross-Section_Measurement_Results.xlsx
  targets         = ../../../../data/lpbf/in718/k1_single_tracks.csv
  case_dir        = cases
  out_dir         = ../../../../results/lpbf/in718/k1
  tag             = flow_h10
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
