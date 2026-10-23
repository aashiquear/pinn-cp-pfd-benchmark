# chi sweep at the author's preferred pair (docs/12 §6.4, author 2026-10-23): W_c 75 MPa, G_c 0.1
# MPa mm, l0 0.6, with the elastic prefactor chi = 0.6 (0 = the calibrated purely ductile run,
# 1 = the first-round form that cracked the fillets at small G_c). Question: does a growing elastic
# contribution move the crack from the mid-gage (0.38) toward the fillets? One run, no search.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/chi_sweep_wc75/chi0.6
  l0        = '0.6'
  Gc_guess  = '0.1'
  first     = '1'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 1
  threads   = 2
  W_c       = 75
  chi       = 0.6
[]
