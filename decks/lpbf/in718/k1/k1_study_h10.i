# Part V, case K1: study deck of the single-track runs at 10 um (coarse check) with the
# CALIBRATED heat source: one Gauss-Newton step of (absorptivity, penetration depth) on the
# relative depth and width errors of case0 and case1.2 (`lpbf_k1.py k1_study_h10_start.i fit`,
# results/lpbf/in718/k1/h10_start/k1_fit_step.json). The step lowers the rms error of the
# four fitted numbers from 10.3 % to a predicted 9.2 %: a conduction model cannot fit the
# depth of both spot sizes (docs/16_lpbf_notebook.md section 5).
#
# Read by scripts/lpbf_k1.py (targets | decks | run | collect | figure | fit). It lists
# the seven AM Bench 2022 single-track cases with their laser settings, names the cases
# the heat source is fitted on (decision D5: two fitted, five held out), and holds the
# heat-source parameters of this pass.
#
# Laser settings per case: power (W), speed (m/s), D4sigma beam diameter (um), as listed
# in the NIST cross-section workbook (weaver 2022 trackdata); the measured depth and width
# are read from that workbook by the `targets` command, not typed here.

[Study]
  template        = k1_template.i
  workbook        = /mnt/extended/projects/ambench-data/nist/DATA_AMB2022_03/mds2-2718/AMB2022-718-SH1-MeltPool_Cross-Section_Measurement_Results.xlsx
  targets         = ../../../../data/lpbf/in718/k1_single_tracks.csv
  case_dir        = cases
  out_dir         = ../../../../results/lpbf/in718/k1
  tag             = h10
  mesh_um         = 10.0
  backend         = torch
  use_gpu         = true
  absorptivity    = 0.7724        # calibrated: one Gauss-Newton step from 0.70 on case0 and case1.2
  penetration_um  = 45.01         # calibrated: the same step, from 52.0
  liquid_conductivity_factor = 1.0   # 1: Wiedemann-Franz values; > 1 stands in for pool convection
  fit_cases       = 'case0 case1.2'
  end_time        = 1.2e-3
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
