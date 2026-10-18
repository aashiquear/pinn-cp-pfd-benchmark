# Ni CP calibration, Gauss-Newton iteration 2: centre from iteration 1's step and one
# multiplicative perturbation per free parameter (scripts/ni_gauss_newton.py).
[Study]
  base       = cal_center2.i
  out_root   = ../../../../results/ni_study/cal_h020/gn2
  markers    = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 5
[]
[Cases]
  [./rho0_p]
    rho0 = 870691
  []
  [./k1_p]
    k1 = 0.0137612
  []
  [./y_c_p]
    y_c = 1.85051e-06
  []
  [./tau_c0_p]
    tau_c0 = 0.365955
  []
[]
