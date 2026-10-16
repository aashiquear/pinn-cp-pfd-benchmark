# Ni CP calibration, Gauss-Newton iteration 3: centre from iteration 2's step and one
# multiplicative perturbation per free parameter (scripts/ni_gauss_newton.py).
[Study]
  base       = cal_center3.i
  out_root   = ../../../results/ni_study/cal/gn3
  markers    = ../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 5
[]
[Cases]
  [./rho0_p]
    rho0 = 916627
  []
  [./k1_p]
    k1 = 0.0129672
  []
  [./y_c_p]
    y_c = 1.77499e-06
  []
  [./tau_c0_p]
    tau_c0 = 0.373733
  []
[]
