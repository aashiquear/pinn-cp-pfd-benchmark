# FE Gauss-Newton iteration 1 at h40: the surrogate optimum (center) and one
# multiplicative perturbation per free parameter (x1.2, tau_c0 x1.5), from which
# the finite-difference Jacobian of the marker stresses is built.
[Study]
  base       = cal_center1.i
  out_root   = ../../results/cp_study/cal/gn1
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 5
[]
[Cases]
  [./rho0_p]   rho0   = 44.63   []
  [./k1_p]     k1     = 0.15074 []
  [./y_c_p]    y_c    = 0.0032  []
  [./tau_c0_p] tau_c0 = 8.1375  []
[]
