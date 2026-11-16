# 2026-11-03: rerun from scratch on the PINNACLE repeated-window damage-start fix (ADR 0034 addendum 2026-11-03); the earlier runs are in <out_root>/../superseded_20261103/.
# Per-orientation W_c (docs/12 §8.7): the other thickness variants of group L with W_c = 75 MPa, on the
# PINNACLE fix of ADR 0034 addendum 2026-10-31 with the linear predictor.
# Driver scripts/ni_study.py prepare run; then ni_specimen.py <case>/specimen.i curve fracture.
[Study]
  base         = ../final/fem_base_h020.i
  out_root     = ../../../results/ni_study/wc/variants_pred
  concurrent   = 2
  include_base = false
[]
[Cases]
  [./L.3]
    group     = L
    specimen  = L.3
    mesh      = ../../../results/ni_study/mesh/h020/L.3/mesh.e
    W_c       = 75
    predictor = linear
  []
  [./L.4]
    group     = L
    specimen  = L.4
    mesh      = ../../../results/ni_study/mesh/h020/L.4/mesh.e
    W_c       = 75
    predictor = linear
  []
[]
