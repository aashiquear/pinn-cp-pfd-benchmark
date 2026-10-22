# NOT RUN: l0 0.6 and 0.8 hit the target in round 3 (G_c 30.95 / 40.75), so no resume was needed.
# Resume of decks/ni/pfd/cal_L1_h020.i for l0 = 0.6 and 0.8 (2026-10-20): the first driver was
# stopped after it misjudged l0.4/gc15.56 (failed at 0.339, but its summary had no failure
# strain because separation and the mechanics failure came in the same step; fixed in
# ni_specimen.py). Same search, same out_root: the finished runs are reused, the search
# continues from them (started once their round-3 runs gc30.95 / gc40.75 ended and were refreshed).
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/cal_L1_h020
  l0        = '0.6 0.8'
  Gc_guess  = '16 22'
  first     = '0.5 2'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 6
  threads   = 2
[]
