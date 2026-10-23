# Quick check, not for reporting (author 2026-10-23): W_c 75, G_c 30 (the W_c = 0 order of G_c),
# chi 0.55, l0 0.6. With G_c / l0 = 50 MPa the elastic share cannot crack the fillets (AT2 strength
# ~ 0.55^0.5 x 750 MPa) and the plastic work past W_c reaches at most ~2 MPa by 47%, so the
# question is what the mixed driving force does in the gage at a large G_c. One run.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/check_wc75_gc30_chi0.55
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
