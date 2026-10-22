# W_c stage of the Ni CP+PFD calibration on L.1 h020 (docs/12 §6.3; author 2026-10-20): l0 fixed at
# the plateau (0.6 mm, l0/h = 3), the plastic-work threshold W_c held at the Fig. 3 softening onset
# (the gage peak plastic work at the 45% marker, the last one the calibrated CP curve matches:
# about 75 MPa, check run fields), and G_c searched to the failure strain 0.463. With damage
# starting at 45% the section must fail within ~1.3% strain, so G_c is far below the W_c = 0
# value: AT2 homogeneous estimate d = 0.85 at 2 H l0 / G_c = 5.7 with H ~ 2 MPa gives G_c ~ 0.4.
# NOT LAUNCHED: the author decides the W_c value(s) first (docs/12 §6.3).
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/wc_L1_h020/wc75
  l0        = '0.6'
  Gc_guess  = '0.6'
  first     = '0.5 2'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 6
  threads   = 2
  W_c       = 75
[]
