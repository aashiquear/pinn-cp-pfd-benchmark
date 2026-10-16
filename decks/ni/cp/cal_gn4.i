# Ni CP calibration, Gauss-Newton iteration 4: centre from iteration 3's step and one
# multiplicative perturbation per free parameter (scripts/ni_gauss_newton.py).
[Study]
  base       = cal_center4.i
  out_root   = ../../../results/ni_study/cal/gn4
  markers    = ../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 5
[]
[Cases]
  [./rho0_p]
    rho0 = 822925
  []
  [./k1_p]
    k1 = 0.0133059
  []
  [./y_c_p]
    y_c = 2.71759e-06
  []
  [./tau_c0_p]
    tau_c0 = 0.368724
  []
[]
