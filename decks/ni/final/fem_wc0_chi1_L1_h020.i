# One-case check before fem_wc0_chi1_15_h020.i (specimen L.1 must reproduce 0.460) with the
# W_c = 0 baseline (G_c 30.95, chi 1, l0 0.6; author 2026-11-02: the first parameter set, reported
# before the W_c calibration on L and the per-orientation W_c). Driver scripts/ni_study.py
# (prepare | run | collect); then ni_specimen.py ... curve fracture per case for the failure
# strain and the fracture location. Compared with the paper's experiments only (Fig. 4(b), Fig. 5).
[Study]
  base       = fem_wc0_chi1_base_h020.i
  out_root   = ../../../results/ni_study/final/fem_wc0_chi1
  concurrent = 6
  include_base = false
[]
[Cases]
  [./L.1]
    group    = L
    specimen = L.1
    mesh     = ../../../results/ni_study/mesh/h020/L.1/mesh.e
  []
[]
