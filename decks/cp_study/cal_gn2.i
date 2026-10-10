# FE Gauss-Newton iteration 2 at h40: center from iteration 1's step and one
# multiplicative perturbation per free parameter.
[Study]
  base       = cal_center2.i
  out_root   = ../../results/cp_study/cal/gn2
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 5
[]
[Cases]
  [./rho0_p]
    rho0 = 49.3825
  []
  [./k1_p]
    k1 = 0.201252
  []
  [./y_c_p]
    y_c = 0.00383692
  []
  [./tau_c0_p]
    tau_c0 = 8.02388
  []
[]
