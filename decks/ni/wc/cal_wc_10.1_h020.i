# Per-orientation W_c calibration (docs/12 §8.7; author 2026-10-29): one specimen per orientation
# (10.1); the calibrated W_c then runs the group's other thickness variants. G_c, l0 and chi are the
# carried set's; W_c is searched to the experiment's failure strain of 10.1 (Fig. 4(b): 0.5404).
# First guess 55 MPa from the plastic work of the stopped W_c = 75 campaign's curve at the target
# (scaled by L.1's 75 / 83); first round x0.9 and x1.1 side by side.
[Calibrate]
  base       = base_10.1_h020.i
  out_root   = ../../../results/ni_study/wc/10.1
  search     = W_c
  l0         = '0.6'
  W_c_guess  = '55'
  Gc         = 0.1
  chi        = 0
  first      = '0.9 1.1'
  outward    = 1.25
  target     = 0.5404
  tolerance  = 0.004
  max_runs   = 6
  threads    = 2
[]
