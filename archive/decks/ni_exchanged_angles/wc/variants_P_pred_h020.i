# Per-orientation W_c (docs/12 §8.7): the other thickness variants of group P with W_c = 82.68 MPa, on the
# PINNACLE fix of ADR 0034 addendum 2026-10-31 with the linear predictor.
# Driver scripts/ni_study.py prepare run; then ni_specimen.py <case>/specimen.i curve fracture.
[Study]
  base         = ../final/fem_base_h020.i
  out_root     = ../../../results/ni_study/wc/variants_pred
  concurrent   = 2
  include_base = false
[]
[Cases]
  [./P.3]
    group     = P
    specimen  = P.3
    mesh      = ../../../results/ni_study/mesh/h020/P.3/mesh.e
    W_c       = 82.68
    predictor = linear
  []
  [./P.4]
    group     = P
    specimen  = P.4
    mesh      = ../../../results/ni_study/mesh/h020/P.4/mesh.e
    W_c       = 82.68
    predictor = linear
  []
[]
