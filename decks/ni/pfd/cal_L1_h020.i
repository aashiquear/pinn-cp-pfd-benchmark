# Ni CP+PFD calibration on specimen L.1, h020 (docs/12 §6; author 2026-10-19): G_c calibrated to
# L.1's experimental failure strain of Fig. 4(b) (0.463, gauge strain) at three length scales,
# l0 = 0.4 / 0.6 / 0.8 mm (2h / 3h / 4h of the 0.2 mm gage cell), each search in its own thread;
# l0 is then chosen by L.1's fracture shape (Fig. 5: inclined ~45 deg at ~0.20 of the gage).
# First estimates from the check run at G_c 125 (docs/12 §6.1: d = 0.85 at the peak point at
# 46.3% needs about 11 MPa mm at l0 0.4, about proportional to l0); the first round brackets each
# estimate by x0.5 and x2 (six runs side by side), then regula falsi in ln G_c.
# Driver: scripts/ni_pfd_calibrate.py.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/cal_L1_h020
  l0        = '0.4 0.6 0.8'
  Gc_guess  = '11 16 22'
  first     = '0.5 2'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 6
  threads   = 2
[]
