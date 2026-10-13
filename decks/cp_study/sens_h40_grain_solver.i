# Parameter sensitivity of the verified PINN-as-solver route with grain-aware
# inputs (mode = fem_residual, node_features, corrector = auto, tol 0.05) at
# h40: the same one-at-a-time set as sens_h40.i run through
# hybrid/h40_fem_solver_grain.i. The base case is
# results/cp_study/hybrid/h40_fem_solver_grain. Measures: curve and fields
# against the FE case of sens/h40, steps solved without Newton, training cost,
# against the parameter change.
[Study]
  base       = hybrid/h40_fem_solver_grain.i
  out_root   = ../../results/cp_study/sens_neural/h40_grain_solver
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 2
[]
[Cases]
  [./rho0_x0p5]
    rho0 = 23.9936
  []
  [./rho0_x2]
    rho0 = 95.9744
  []
  [./k1_x0p5]
    k1 = 0.0994165
  []
  [./k1_x2]
    k1 = 0.397666
  []
  [./y_c_x0p5]
    y_c = 0.00171879
  []
  [./y_c_x2]
    y_c = 0.00687516
  []
  [./tau_c0_x0p5]
    tau_c0 = 2.64574
  []
  [./tau_c0_x2]
    tau_c0 = 10.5829
  []
  [./rho0_x0p8]
    rho0 = 38.3898
  []
  [./rho0_x1p25]
    rho0 = 59.984
  []
  [./k1_x0p8]
    k1 = 0.159066
  []
  [./k1_x1p25]
    k1 = 0.248541
  []
  [./y_c_x0p8]
    y_c = 0.00275006
  []
  [./y_c_x1p25]
    y_c = 0.00429697
  []
  [./tau_c0_x0p8]
    tau_c0 = 4.23318
  []
  [./tau_c0_x1p25]
    tau_c0 = 6.61434
  []
[]
