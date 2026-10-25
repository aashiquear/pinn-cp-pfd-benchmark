# FE runs of the surrogate design's 30 synthetic specimens (written by scripts/ni_surrogate_design.py;
# author 2026-10-25): the carried CP+PFD set of the base deck, one case per sample. Driver scripts/ni_study.py.
[Study]
  base         = ../../../../../decks/ni/final/fem_base_h020.i
  out_root     = ./runs
  concurrent   = 6
  include_base = false
[]
[Cases]
  [./S01]
    group        = S01
    specimen     = S01
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S01/mesh.e
  []
  [./S02]
    group        = S02
    specimen     = S02
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S02/mesh.e
  []
  [./S03]
    group        = S03
    specimen     = S03
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S03/mesh.e
  []
  [./S04]
    group        = S04
    specimen     = S04
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S04/mesh.e
  []
  [./S05]
    group        = S05
    specimen     = S05
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S05/mesh.e
  []
  [./S06]
    group        = S06
    specimen     = S06
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S06/mesh.e
  []
  [./S07]
    group        = S07
    specimen     = S07
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S07/mesh.e
  []
  [./S08]
    group        = S08
    specimen     = S08
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S08/mesh.e
  []
  [./S09]
    group        = S09
    specimen     = S09
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S09/mesh.e
  []
  [./S10]
    group        = S10
    specimen     = S10
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S10/mesh.e
  []
  [./S11]
    group        = S11
    specimen     = S11
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S11/mesh.e
  []
  [./S12]
    group        = S12
    specimen     = S12
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S12/mesh.e
  []
  [./S13]
    group        = S13
    specimen     = S13
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S13/mesh.e
  []
  [./S14]
    group        = S14
    specimen     = S14
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S14/mesh.e
  []
  [./S15]
    group        = S15
    specimen     = S15
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S15/mesh.e
  []
  [./S16]
    group        = S16
    specimen     = S16
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S16/mesh.e
  []
  [./S17]
    group        = S17
    specimen     = S17
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S17/mesh.e
  []
  [./S18]
    group        = S18
    specimen     = S18
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S18/mesh.e
  []
  [./S19]
    group        = S19
    specimen     = S19
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S19/mesh.e
  []
  [./S20]
    group        = S20
    specimen     = S20
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S20/mesh.e
  []
  [./S21]
    group        = S21
    specimen     = S21
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S21/mesh.e
  []
  [./S22]
    group        = S22
    specimen     = S22
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S22/mesh.e
  []
  [./S23]
    group        = S23
    specimen     = S23
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S23/mesh.e
  []
  [./S24]
    group        = S24
    specimen     = S24
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S24/mesh.e
  []
  [./S25]
    group        = S25
    specimen     = S25
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S25/mesh.e
  []
  [./S26]
    group        = S26
    specimen     = S26
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S26/mesh.e
  []
  [./S27]
    group        = S27
    specimen     = S27
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S27/mesh.e
  []
  [./S28]
    group        = S28
    specimen     = S28
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S28/mesh.e
  []
  [./S29]
    group        = S29
    specimen     = S29
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S29/mesh.e
  []
  [./S30]
    group        = S30
    specimen     = S30
    orientations = ./orientations.csv
    thickness    = ./thickness.csv
    mesh         = ./mesh/S30/mesh.e
  []
[]
