# Rerun of wc75_gc30_chi0.55_L1_h020.i under a 40 GB resident cap (2026-10-23): the first run was oom-killed in
# the backward-Euler fallback's sub-stepping (exit -9), with d 0.011 (check) / at 1.5% strain (chi 0.8),
# i.e. a mechanics death, not a result. Same deck otherwise; the dead run is kept.
[Calibrate]
  base      = base_L1_h020_cap40.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/check_wc75_gc30_chi0.55_cap40
  l0        = '0.6'
  Gc_guess  = '30'
  first     = '1'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 1
  threads   = 2
  W_c       = 75
  chi       = 0.55
[]
