# J2+PFD Gc calibration at W_c = 200 MPa, l0 = 5 um (W_c scan of
# 2026-10-14: the author asked for a threshold that delays the softening to
# the end of the hardening range and a lower Gc; docs/10 §4.10).
[Calibrate]
  base          = base_h40_wc200.i
  out_root      = ../../../results/pfd_study/j2/cal_h40_wc200
  l0            = '5 10 20'   # l0 5 done in the scan (its trials are reused), 10 and 20 for the Gc(l0) figure
  Gc_lo         = 10
  Gc_hi         = 20000
  target_strain = 0.172
  criterion     = node
  tolerance     = 0.0008
  max_iter      = 9
  concurrent    = 3
  threads       = 2
[]
