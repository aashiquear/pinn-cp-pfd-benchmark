# BASE specimen deck, route 'hybrid' (PINN+FEM), restart of Part III (2026-11-15): fem_base_h020.i
# with the load-step network trained on the assembled weak-form residual (PINNACLE ADR 0029,
# mode fem_residual, on the total-Lagrangian assembler of ADR 0035) as the predictor and Newton
# as the corrector of every step. Single crystal: no grain features.
[Specimen]
  out_dir          = .
  mesh             = ../../../../mesh/h020/10.4/mesh.e
  orientations     = ../../../../../../data/extracted/ni_table1_orientations.csv
  group            = 10
  thickness        = ../../../../../../data/extracted/ni_table2_thickness.csv
  specimen         = 10.4
  gage_width       = 1.85
  length           = 25.4
  strain = 0.40
  t_end = 800          # gauge strain rate about 1e-3/s, as the CP calibration (600 s for 30 % bar strain)
  n_steps          = 400          # same 1e-3 grip strain per step
  snapshots_every  = 2
  backend          = torch
  use_gpu          = false
  threads = 4
  newton_tol       = 5.0e-2   # 2026-10-26: Newton floor of the implicit stress's FD tangent on multi-slip orientations (S01 1e-2 N); 5e-4 of the reaction; reactions unchanged at 5e-3 to 1e-7
  newton_rel_tol   = 1.0e-6
  max_cutbacks     = 10
  linear_method    = direct
  sparse_threshold = 500
  kinematics       = total_lagrangian
  state_update     = implicit
  state_iterations = 2
  C11     = 253000
  C12     = 152000
  C44     = 124000
  mu      = 26000
  gamma0  = 1.0e-5
  nexp    = 8.333333333333334
  tau_c0 = 0.26
  q_lat   = 1.4
  alpha0  = 0.4
  burgers = 2.5e-7
  rho0 = 2.0257e6
  k1 = 0.0243375
  y_c = 4.6689e-07
  predictor = fem_solver_grain
  predictor_features = ''
  predictor_device = cpu
  predictor_gpu = false
  corrector = newton
  memory_cap_gb    = 16
  Gc               = 49.37     # calibrated on L.1 (decks/ni/pfd/tl/cal_L1_h020.i: failure strain 0.4635 against 0.463)
  l0               = 0.6
  zeta             = 1          # plastic-work driving force prefactor
  W_c = 0
  eta              = 0
  chi = 1
  damage_model     = at2
  stop_damage      = 0
  checkpoint_every = 8
[]
