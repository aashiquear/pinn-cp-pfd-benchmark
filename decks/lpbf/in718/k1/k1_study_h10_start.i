# Part V, case K1: study deck of the single-track runs at 10 um with the STARTING heat source
# (absorptivity 0.70, penetration 52 um: the values the author's toolkit fitted with its own code).
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
  tag             = h10_start
  mesh_um         = 10.0
  backend         = torch
  use_gpu         = true
  absorptivity    = 0.70          # starting value of the fit
  penetration_um  = 52.0          # starting value of the fit
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
