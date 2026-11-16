# Check after the Euler-convention fix (PINNACLE ADR 0036): specimens 9.1 and H.1, the two groups
# that failed 19 to 26 % early with the exchanged angles, with the single set refitted on L.1
# (fem_base_h020.i). Author's order of 2026-11-15: these two first, then decide.
[Study]
  base       = fem_base_h020.i
  out_root   = ../../../../results/ni_study/bunge/final/fem
  concurrent = 2
  include_base = false
[]
[Cases]
  [./9.1]
    group    = 9
    specimen = 9.1
    mesh     = ../../../../results/ni_study/mesh/h020/9.1/mesh.e
  []
  [./H.1]
    group    = H
    specimen = H.1
    mesh     = ../../../../results/ni_study/mesh/h020/H.1/mesh.e
  []
[]
