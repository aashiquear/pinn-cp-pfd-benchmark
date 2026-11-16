# G_c calibrated to L.1's experimental failure strain (0.463, Fig. 4(b)) at l0 = 0.6 mm, on the
# total-Lagrangian balance. First round brackets the guess; then regula falsi in ln G_c.
# Driver: scripts/ni_pfd_calibrate.py.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../../results/ni_study/bunge/pfd/cal_L1_h020
  l0        = '0.6'
  Gc_guess  = '50'
  first     = '0.85 1.15'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 7
  threads   = 4
  chi       = 1
  W_c       = 0
[]
