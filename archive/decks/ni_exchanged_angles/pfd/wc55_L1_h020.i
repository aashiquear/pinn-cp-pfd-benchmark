# W_c stage of the Ni CP+PFD calibration on L.1 h020 (docs/12 §6.3; author 2026-10-21): l0 fixed
# at 0.6 mm (l0/h = 3, the plateau), the plastic-work threshold W_c held at 55 MPa and G_c searched
# to the failure strain 0.463. Three such decks (W_c 55 / 65 / 75) run side by side to show the
# W_c trend; the gage-peak plastic work (check run, W_c = 0) is 61 MPa at 38%, 73.5 at 44%,
# 76 at 46.3%, so W_c = 55 starts damage near the strain where the peak reaches it, and a
# W_c above ~76 could not fail at the target. G_c guess from the W_c = 0 scaling
# (failure when W_p - W_c ~ 1.5 G_c / l0): 8 MPa mm, bracketed x0.5 / x2 first.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/wc55
  l0        = '0.6'
  Gc_guess  = '8'
  first     = '0.5 2'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 6
  threads   = 2
  W_c       = 55
[]
