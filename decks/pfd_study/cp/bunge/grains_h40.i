# Recheck of 2026-11-16: base = base_h40.i in this folder (recalibrated CP set, Bunge angles,
# total-Lagrangian balance, implicit update, Gc 1122.66 at l0 5). Text below as in the first study.
# CP+PFD grain-instance check at h40 (docs/09 §5.1 step 5): the calibrated
# case (Gc 1106, l0 5) on the two other Part I grain structures (instances
# sB and sC: new Voronoi seeding grown the same way, new orientation draws,
# docs/07 / report §7), everything else unchanged. The base is instance A.
[Study]
  base       = base_h40.i
  out_root   = ../../../../results/pfd_study/cp/bunge/grains_h40
  concurrent = 2
[]
[Cases]
  [./sB]
    points = ../../../../results/cp_study/bunge/grain_instance/sB/orientations_points.csv
  []
  [./sC]
    points = ../../../../results/cp_study/bunge/grain_instance/sC/orientations_points.csv
  []
[]
