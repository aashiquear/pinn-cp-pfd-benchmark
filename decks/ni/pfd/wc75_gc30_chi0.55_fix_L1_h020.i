# Third attempt, after the PINNACLE fix of the fallback memory runaway (ADR 0033 addendum 2026-10-24),
# under the original 16 GB cap (the first two attempts, 16 and 40 GB, were oom-killed at step 72).
# Same deck otherwise; the dead runs are kept.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/check_wc75_gc30_chi0.55_fix
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
