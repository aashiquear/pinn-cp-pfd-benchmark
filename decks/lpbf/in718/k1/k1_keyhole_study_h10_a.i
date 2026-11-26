# Part V, case K1: single tracks with the heat-source depth taken from the keyhole number.
#
# STUDY DECK for scripts/lpbf_k1.py (template k1_keyhole_template.i). The seven NIST
# AMB2022-03 single tracks at 10 um, conduction only. The decay length of the source is
# keyhole_depth_share x e(Ke) per track (PINNACLE ADR 0055) in place of the fitted
# penetration depth of k1_study_h10.i (45.01 um). The absorptivity of the source is kept
# at the value of that fit (0.7724); no parameter is fitted to a depth here.
# Variant A: the keyhole number uses the absorptivity of the source (0.7724, the value of the conduction fit).
#
# Reported beside the fitted-depth study, not in place of it.

[Study]
  template        = k1_keyhole_template.i
  workbook        = /mnt/extended/projects/ambench-data/nist/DATA_AMB2022_03/mds2-2718/AMB2022-718-SH1-MeltPool_Cross-Section_Measurement_Results.xlsx
  targets         = ../../../../data/lpbf/in718/k1_single_tracks.csv
  case_dir        = cases
  out_dir         = ../../../../results/lpbf/in718/k1
  tag             = keyhole_h10_a
  mesh_um         = 10.0
  backend         = torch
  use_gpu         = true
  absorptivity    = 0.7724        # calibrated: one Gauss-Newton step from 0.70 on case0 and case1.2
  penetration_um  = 10.0          # NOT fitted: one element, the smallest decay length
  keyhole_absorptivity = 0.7724
  keyhole_depth_share  = 0.5
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
