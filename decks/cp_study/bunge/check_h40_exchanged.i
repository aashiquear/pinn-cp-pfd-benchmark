# Al2219 recheck after the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): the calibrated
# Part I set (cal_center2.i) on the h40 plate, one case before any campaign.
# Exchanged angles (phi2, Phi, phi1) of seed 20261009: the crystals the Part I runs simulated;
# legacy balance, explicit state update. Reproduces results/cp_study/cal/gn2/base with the fixed code.
[Growth]
  run_dir = ../../../results/graingrowth_pilot/n32_bdf1_h40
  deck    = ../../graingrowth/pilot/n32_bdf1_h40.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  file = grains_seed20261009_exchanged.csv
[]

[Mechanics]
  out_dir      = ../../../results/cp_study/bunge/check/exchanged
  load_axis    = x
  strain       = 0.20
  t_end        = 200.0
  n_steps      = 256
  nx           = 40
  ny           = 20
  thickness    = 1.0
  newton_tol   = 1.0e-3
  newton_rel_tol = 0
  max_cutbacks = 10
  linear_method = direct
  snapshots_every = 4
  C11     = 108.0e3
  C12     = 61.3e3
  C44     = 28.5e3
  mu      = 25.7e3
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 5.29147
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 47.9872
  k1      = 0.198833
  y_c     = 0.00343758
[]
