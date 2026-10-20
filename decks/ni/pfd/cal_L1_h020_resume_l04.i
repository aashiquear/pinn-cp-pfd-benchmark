# Resume of decks/ni/pfd/cal_L1_h020.i for l0 = 0.4 only (2026-10-20): the first driver was
# stopped after it misjudged l0.4/gc15.56 (failed at 0.339, but its summary had no failure
# strain because separation and the mechanics failure came in the same step; fixed in
# ni_specimen.py). Same search, same out_root: the finished runs are reused, the search
# continues from them. l0 = 0.6 / 0.8 resume with the full deck once their round-3 runs end.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/cal_L1_h020
  l0        = '0.4'
  Gc_guess  = '11'
  first     = '0.5 2'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 6
  threads   = 2
[]
