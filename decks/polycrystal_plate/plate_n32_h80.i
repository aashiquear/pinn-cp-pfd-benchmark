# Pipeline deck: 32-grain phase-field polycrystal (pilot horizon h80) on the
# MOOSE-reference plate, crystal plasticity only. Consumed by
# scripts/polycrystal_plate.py; every stage reads its inputs from here.
#
# Mechanics values are the reference MOOSE deck's (moose/CPFE.i,
# CrystalPlasticityDislocationUpdate) in um / MPa / s:
#   ao = 1e-5 -> gamma0; xm = 0.12 -> n = 1/0.12; tau_c_0 = 0.112; k_0 -> k1;
#   y_c = 0.00286 (NOTE: Rezwan 2025 Table 2 gives 0.00268; the deck is taken
#   as the reference and the difference is recorded); init_rho_ssd = 220.
#   Loading: 0.15 um/s on the 50 um height for 10 s = 3.0% at 3e-3 /s, 40
#   steps of 0.25 s, as the MOOSE run that produced moose/CPFE_out.e.

[Growth]
  run_dir = ../../results/graingrowth_pilot/n32_bdf1_h80
  deck    = ../graingrowth/pilot/n32_bdf1_h80.i
[]

[Labeling]
  threshold = 0.5
  min_cells = 4
[]

[Orientations]
  seed = 20261009
[]

[Mechanics]
  out_dir      = ../../results/polycrystal_plate/n32_h80
  strain       = 0.03
  t_end        = 10.0
  n_steps      = 40
  nx           = 50
  ny           = 25
  thickness    = 1.0
  newton_tol   = 1.0e-6
  max_cutbacks = 10
  C11     = 108.0e3
  C12     = 61.3e3
  C44     = 28.5e3
  mu      = 25.7e3
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0  = 0.112
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 0.000286
  rho0    = 220.0
  k1      = 0.3432
  y_c     = 0.00286
[]
