# Parameter sensitivity of the PINN+FEM route (linear-extrapolation predictor, Newton corrector)
# at h40, recheck of 2026-11-16: the one-at-a-time set of sens_h40.i in this folder
# run through hybrid/h40_linear.i. The base case is results/cp_study/bunge/hybrid/h40_linear.
[Study]
  base       = hybrid/h40_linear.i
  out_root   = ../../../results/cp_study/bunge/sens_neural/h40_linear
  markers    = ../../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 8
  threads    = 2
[]
[Cases]
  [./rho0_x0p5]
    rho0 = 26.9653
  []
  [./rho0_x0p8]
    rho0 = 43.1445
  []
  [./rho0_x1p25]
    rho0 = 67.4133
  []
  [./rho0_x2]
    rho0 = 107.861
  []
  [./k1_x0p5]
    k1 = 0.0883855
  []
  [./k1_x0p8]
    k1 = 0.141417
  []
  [./k1_x1p25]
    k1 = 0.220964
  []
  [./k1_x2]
    k1 = 0.353542
  []
  [./y_c_x0p5]
    y_c = 0.00135712
  []
  [./y_c_x0p8]
    y_c = 0.00217138
  []
  [./y_c_x1p25]
    y_c = 0.00339279
  []
  [./y_c_x2]
    y_c = 0.00542846
  []
  [./tau_c0_x0p5]
    tau_c0 = 2.69478
  []
  [./tau_c0_x0p8]
    tau_c0 = 4.31165
  []
  [./tau_c0_x1p25]
    tau_c0 = 6.73695
  []
  [./tau_c0_x2]
    tau_c0 = 10.7791
  []
[]
