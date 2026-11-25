# Part V, case K1: study deck of the STANDALONE single-track decks at 5 um.
#
# `scripts/lpbf_k1.py k1_study_h5.i decks` writes the seven decks under cases/h5/; each
# is a complete PINNACLE input (588,000 elements, about 1,900 explicit steps):
#     pinnacle run decks/lpbf/in718/k1/cases/h5/case0.i
# and `run`, `collect`, `figure` compare them with the NIST cross sections. The heat
# source carries the 10 um calibration (absorptivity 0.7724, penetration 45.01 um); a fit
# at this mesh is `lpbf_k1.py k1_study_h5.i fit`. For several GPUs add
# `devices = 'cuda:0 cuda:1 ...'` (or `num_gpus = N`) to [Global] of the generated decks,
# or set them in the template.
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
  tag             = h5
  mesh_um         = 5.0
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
