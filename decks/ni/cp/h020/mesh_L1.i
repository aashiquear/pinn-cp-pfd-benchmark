# Ni Part III: mesh sensitivity on specimen L.1 with the set calibrated on h020
# (decks/ni/cp/h020/calibrated_L1_h020.i). The h020 case is the calibrated run itself
# (results/ni_study/cal_h020/gn2/base, copied as <out_root>/base); the cases change only
# the mesh: h040 (0.4 mm gage cell, the reference the author keeps to show the coarse
# mesh's error), h010 (0.1 mm) and h020 with two elements through the thickness.
# Grip strain 30% in 300 steps (1e-3 per step, 5e-2/s on the bar). PINNACLE 9f31bb0.
# Threads 2 / 8 / 4; memory cap 16 GB per run, 24 GB for h010 (4x the points of h020).
[Study]
  base       = calibrated_L1_h020.i
  out_root   = ../../../../results/ni_study/mesh_study/L1_h020set
  markers    = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 3
[]
[Cases]
  [./h040]
    mesh    = ../../../../results/ni_study/mesh/h040/L.1/mesh.e
    threads = 2
  []
  [./h010]
    mesh          = ../../../../results/ni_study/mesh/h010/L.1/mesh.e
    threads       = 8
    memory_cap_gb = 24
  []
  [./h020z2]
    mesh    = ../../../../results/ni_study/mesh/h020z2/L.1/mesh.e
    threads = 4
  []
[]
