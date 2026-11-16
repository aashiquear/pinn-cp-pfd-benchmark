# Per-orientation W_c (docs/12 §8.7; author 2026-10-29): the other thickness variants of group L with
# L's calibrated W_c (75 MPa on L.1, docs/12 §6.4), campaign settings of fem_base_h020.i.
# Driver scripts/ni_study.py prepare run; then ni_specimen.py <case>/specimen.i curve fracture.
[Study]
  base         = ../final/fem_base_h020.i
  out_root     = ../../../results/ni_study/wc/variants
  concurrent   = 2
  include_base = false
[]
[Cases]
  [./L.3]
    group    = L
    specimen = L.3
    mesh     = ../../../results/ni_study/mesh/h020/L.3/mesh.e
    W_c      = 75
  []
  [./L.4]
    group    = L
    specimen = L.4
    mesh     = ../../../results/ni_study/mesh/h020/L.4/mesh.e
    W_c      = 75
  []
[]
