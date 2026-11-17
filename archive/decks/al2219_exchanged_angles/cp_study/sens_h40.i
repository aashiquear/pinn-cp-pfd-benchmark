# Parameter sensitivity of the calibrated set at the chosen mesh (h40; author's
# decision 2026-10-12: within 4.2% of h160 at 14 min per run). One
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
    rho0 = 23.9936
  []
  [./rho0_x0p8]
    rho0 = 38.3898
  []
  [./rho0_x1p25]
    rho0 = 59.984
  []
  [./rho0_x2]
    rho0 = 95.9744
  []
  [./k1_x0p5]
    k1 = 0.0994165
  []
  [./k1_x0p8]
    k1 = 0.159066
  []
  [./k1_x1p25]
    k1 = 0.248541
  []
  [./k1_x2]
    k1 = 0.397666
  []
  [./y_c_x0p5]
    y_c = 0.00171879
  []
  [./y_c_x0p8]
    y_c = 0.00275006
  []
  [./y_c_x1p25]
    y_c = 0.00429697
  []
  [./y_c_x2]
    y_c = 0.00687516
  []
  [./tau_c0_x0p5]
    tau_c0 = 2.64574
  []
  [./tau_c0_x0p8]
    tau_c0 = 4.23318
  []
  [./tau_c0_x1p25]
    tau_c0 = 6.61434
  []
  [./tau_c0_x2]
    tau_c0 = 10.5829
  []
[]
