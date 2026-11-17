# Last stage, route 'fem': the 15 specimens (5 orientation groups x 3 thickness variants) with the
# W_c = 0 set (G_c 30.95, chi 0, l0 0.6; docs/12 §8.3, author 2026-10-28:
# W_c = 75 calibrated on L transfers late on the other orientations). Driver scripts/ni_study.py
# (prepare | run | collect); then ni_specimen.py ... curve fracture per case for the failure
# strain and the fracture location. Compared with the paper's experiments only (Fig. 4(b), Fig. 5).
[Study]
  base       = fem_wc0_base_h020.i
  out_root   = ../../../results/ni_study/final/fem_wc0
  concurrent = 6
  include_base = false
[]
[Cases]
  [./9.1]
    group    = 9
    specimen = 9.1
    mesh     = ../../../results/ni_study/mesh/h020/9.1/mesh.e
  []
  [./9.3]
    group    = 9
    specimen = 9.3
    mesh     = ../../../results/ni_study/mesh/h020/9.3/mesh.e
  []
  [./9.4]
    group    = 9
    specimen = 9.4
    mesh     = ../../../results/ni_study/mesh/h020/9.4/mesh.e
  []
  [./10.1]
    group    = 10
    specimen = 10.1
    mesh     = ../../../results/ni_study/mesh/h020/10.1/mesh.e
  []
  [./10.3]
    group    = 10
    specimen = 10.3
    mesh     = ../../../results/ni_study/mesh/h020/10.3/mesh.e
  []
  [./10.4]
    group    = 10
    specimen = 10.4
    mesh     = ../../../results/ni_study/mesh/h020/10.4/mesh.e
  []
  [./H.1]
    group    = H
    specimen = H.1
    mesh     = ../../../results/ni_study/mesh/h020/H.1/mesh.e
  []
  [./H.3]
    group    = H
    specimen = H.3
    mesh     = ../../../results/ni_study/mesh/h020/H.3/mesh.e
  []
  [./H.4]
    group    = H
    specimen = H.4
    mesh     = ../../../results/ni_study/mesh/h020/H.4/mesh.e
  []
  [./L.1]
    group    = L
    specimen = L.1
    mesh     = ../../../results/ni_study/mesh/h020/L.1/mesh.e
  []
  [./L.3]
    group    = L
    specimen = L.3
    mesh     = ../../../results/ni_study/mesh/h020/L.3/mesh.e
  []
  [./L.4]
    group    = L
    specimen = L.4
    mesh     = ../../../results/ni_study/mesh/h020/L.4/mesh.e
  []
  [./P.1]
    group    = P
    specimen = P.1
    mesh     = ../../../results/ni_study/mesh/h020/P.1/mesh.e
  []
  [./P.3]
    group    = P
    specimen = P.3
    mesh     = ../../../results/ni_study/mesh/h020/P.3/mesh.e
  []
  [./P.4]
    group    = P
    specimen = P.4
    mesh     = ../../../results/ni_study/mesh/h020/P.4/mesh.e
  []
[]
