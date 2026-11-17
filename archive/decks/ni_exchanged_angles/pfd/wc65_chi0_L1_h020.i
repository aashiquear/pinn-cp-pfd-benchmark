# W_c stage, second round (docs/12 §6.4): purely ductile driving force, chi = 0 (PINNACLE
# elastic_prefactor; H = max(H, zeta <W_p - W_c>)), so the elastic energy no longer cracks the fillets
# at small G_c and W_c alone sets the onset. l0 0.6, W_c 65 MPa, G_c searched to 0.463 from 4
# (x0.5 / x2 first). NOT LAUNCHED until the author's go-ahead.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65
  l0        = '0.6'
  Gc_guess  = '4'
  first     = '0.5 2'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 6
  threads   = 2
  W_c       = 65
  chi       = 0
[]
