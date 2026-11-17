# One-run check of the PINNACLE repeated-window damage-start fix (ADR 0034 addendum 2026-11-03) on the
# calibrated P.1 point: W_c 82.68 failed at gauge 0.4667 on 7645793 (cal_wc_P.1_pred_h020.i) with
# 431 cutbacks, so the old damage ramp may have acted. If this run reproduces it, the P set stands.
[Calibrate]
  base       = base_P.1_pred_h020.i
  out_root   = ../../../results/ni_study/wc/P.1_check
  search     = W_c
  l0         = '0.6'
  W_c_guess  = '82.68'
  Gc         = 0.1
  chi        = 0
  first      = '1'
  outward    = 1.25
  target     = 0.4657
  tolerance  = 0.004
  max_runs   = 1
  threads    = 2
[]
