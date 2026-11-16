# 2026-11-03: rerun from scratch on the PINNACLE repeated-window damage-start fix (ADR 0034 addendum 2026-11-03); the earlier runs are in <out_root>/../superseded_20261103/.
# Rerun of cal_wc_H.1_h020.i on the PINNACLE line-search probe fix with the linear predictor (docs/12 §8.7), from step 0: a checkpoint restart loses the warm-started state update.
# Per-orientation W_c calibration (docs/12 §8.7; author 2026-10-29): one specimen per orientation
# (H.1); the calibrated W_c then runs the group's other thickness variants. G_c, l0 and chi are the
# carried set's; W_c is searched to the experiment's failure strain of H.1 (Fig. 4(b): 0.5489).
# First guess 97 MPa from the plastic work of the stopped W_c = 75 campaign's curve at the target
# (scaled by L.1's 75 / 83); first round x0.9 and x1.1 side by side.
[Calibrate]
  base       = base_H.1_pred_h020.i
  out_root   = ../../../results/ni_study/wc/H.1_pred
  search     = W_c
  l0         = '0.6'
  W_c_guess  = '97'
  Gc         = 0.1
  chi        = 0
  first      = '0.9 1.1'
  outward    = 1.25
  target     = 0.5489
  tolerance  = 0.004
  max_runs   = 6
  threads    = 2
[]
