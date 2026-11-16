# Per-orientation W_c calibration (docs/12 §8.7; author 2026-10-29): one specimen per orientation
# (P.1); the calibrated W_c then runs the group's other thickness variants. G_c, l0 and chi are the
# carried set's; W_c is searched to the experiment's failure strain of P.1 (Fig. 4(b): 0.4657).
# First guess 77 MPa from the plastic work of the stopped W_c = 75 campaign's curve at the target
# (scaled by L.1's 75 / 83); first round x0.9 and x1.1 side by side.
[Calibrate]
  base       = base_P.1_h020.i
  out_root   = ../../../results/ni_study/wc/P.1
  search     = W_c
  l0         = '0.6'
  W_c_guess  = '77'
  Gc         = 0.1
  chi        = 0
  first      = '0.9 1.1'
  outward    = 1.25
  target     = 0.4657
  tolerance  = 0.004
  max_runs   = 6
  threads    = 2
[]
