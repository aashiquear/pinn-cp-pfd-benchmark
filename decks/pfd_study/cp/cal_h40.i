# Gc(l0) calibration of the CP+PFD polycrystal plate at h40 (docs/09 §5.1 step 2): for
# each length scale, bisection in log Gc until fracture initiates (a node at
# d >= 0.8) at the last Vaughan strain, 17.2% nominal, within one load step.
# l0 levels: 2h, 4h and 8h of the h40 background (h = 2.5 um); the brittle
# rule l0 > 2h is the smallest. Runs: scripts/pfd_calibrate.py, on the host,
# three bisections side by side at 2 threads each; the l0 levels and the Gc bracket follow the J2 result.
[Calibrate]
  base          = base_h40.i
  out_root      = ../../../results/pfd_study/cp/cal_h40
  l0            = '5 10 20'
  Gc_lo         = 10
  Gc_hi         = 20000
  target_strain = 0.172
  criterion     = node
  tolerance     = 0.0008
  max_iter      = 9
  concurrent    = 3
  threads       = 2
[]
