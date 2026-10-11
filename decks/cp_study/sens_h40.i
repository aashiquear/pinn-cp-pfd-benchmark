# Parameter sensitivity of the calibrated set at the chosen mesh (h40; author's
# decision 2026-10-11: within 5% of h160 at a fraction of the cost). One
# parameter at a time, four multiples spanning a factor 4 (x0.5, x0.8, x1.25,
# x2.0), the other three at their calibrated values; the base case is the
# calibrated run itself. Collected against the experiment markers.
[Study]
  base       = cal_center2.i
  out_root   = ../../results/cp_study/sens/h40
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 10
[]
[Cases]
  [./rho0_x0p5]
    rho0 = 20.576
  []
  [./rho0_x0p8]
    rho0 = 32.9217
  []
  [./rho0_x1p25]
    rho0 = 51.4401
  []
  [./rho0_x2]
    rho0 = 82.3042
  []
  [./k1_x0p5]
    k1 = 0.083855
  []
  [./k1_x0p8]
    k1 = 0.134168
  []
  [./k1_x1p25]
    k1 = 0.209638
  []
  [./k1_x2]
    k1 = 0.33542
  []
  [./y_c_x0p5]
    y_c = 0.00159872
  []
  [./y_c_x0p8]
    y_c = 0.00255795
  []
  [./y_c_x1p25]
    y_c = 0.0039968
  []
  [./y_c_x2]
    y_c = 0.00639488
  []
  [./tau_c0_x0p5]
    tau_c0 = 2.67462
  []
  [./tau_c0_x0p8]
    tau_c0 = 4.2794
  []
  [./tau_c0_x1p25]
    tau_c0 = 6.68656
  []
  [./tau_c0_x2]
    tau_c0 = 10.6985
  []
[]
