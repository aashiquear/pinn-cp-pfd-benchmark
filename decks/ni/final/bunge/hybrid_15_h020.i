# Rerun after the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): Bunge angles, set refitted on L.1.
# Fifteen specimens, route 'hybrid', the single parameter set (hybrid_base_h020.i). Driver
# scripts/ni_study.py; compared with the paper's experiments and with the FEM route.
[Study]
  base       = hybrid_base_h020.i
  out_root   = ../../../../results/ni_study/bunge/final/hybrid
  concurrent = 7
  include_base = false
[]
[Cases]
  [./9.1]
    group    = 9
    specimen = 9.1
    mesh     = ../../../../results/ni_study/mesh/h020/9.1/mesh.e
  []
  [./9.3]
    group    = 9
    specimen = 9.3
    mesh     = ../../../../results/ni_study/mesh/h020/9.3/mesh.e
  []
  [./9.4]
    group    = 9
    specimen = 9.4
    mesh     = ../../../../results/ni_study/mesh/h020/9.4/mesh.e
  []
  [./10.1]
    group    = 10
    specimen = 10.1
    mesh     = ../../../../results/ni_study/mesh/h020/10.1/mesh.e
  []
  [./10.3]
    group    = 10
    specimen = 10.3
    mesh     = ../../../../results/ni_study/mesh/h020/10.3/mesh.e
  []
  [./10.4]
    group    = 10
    specimen = 10.4
    mesh     = ../../../../results/ni_study/mesh/h020/10.4/mesh.e
  []
  [./H.1]
    group    = H
    specimen = H.1
    mesh     = ../../../../results/ni_study/mesh/h020/H.1/mesh.e
  []
  [./H.3]
    group    = H
    specimen = H.3
    mesh     = ../../../../results/ni_study/mesh/h020/H.3/mesh.e
  []
  [./H.4]
    group    = H
    specimen = H.4
    mesh     = ../../../../results/ni_study/mesh/h020/H.4/mesh.e
  []
  [./L.1]
    group    = L
    specimen = L.1
    mesh     = ../../../../results/ni_study/mesh/h020/L.1/mesh.e
  []
  [./L.3]
    group    = L
    specimen = L.3
    mesh     = ../../../../results/ni_study/mesh/h020/L.3/mesh.e
  []
  [./L.4]
    group    = L
    specimen = L.4
    mesh     = ../../../../results/ni_study/mesh/h020/L.4/mesh.e
  []
  [./P.1]
    group    = P
    specimen = P.1
    mesh     = ../../../../results/ni_study/mesh/h020/P.1/mesh.e
  []
  [./P.3]
    group    = P
    specimen = P.3
    mesh     = ../../../../results/ni_study/mesh/h020/P.3/mesh.e
  []
  [./P.4]
    group    = P
    specimen = P.4
    mesh     = ../../../../results/ni_study/mesh/h020/P.4/mesh.e
  []
[]
