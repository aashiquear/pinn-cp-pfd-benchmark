# Per-orientation W_c (docs/12 §8.7; author 2026-10-29): the other thickness variants of group 10 with
# group 10's calibrated W_c (49.5 MPa on 10.1, docs/12 §8.7), campaign settings of fem_base_h020.i.
# Driver scripts/ni_study.py prepare run; then ni_specimen.py <case>/specimen.i curve fracture.
[Study]
  base         = ../final/fem_base_h020.i
  out_root     = ../../../results/ni_study/wc/variants
  concurrent   = 2
  include_base = false
[]
[Cases]
  [./10.3]
    group    = 10
    specimen = 10.3
    mesh     = ../../../results/ni_study/mesh/h020/10.3/mesh.e
    W_c      = 49.5
  []
  [./10.4]
    group    = 10
    specimen = 10.4
    mesh     = ../../../results/ni_study/mesh/h020/10.4/mesh.e
    W_c      = 49.5
  []
[]
