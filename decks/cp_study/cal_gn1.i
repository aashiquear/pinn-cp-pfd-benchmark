# FE Gauss-Newton iteration 1 at h40: the center of cal_center1.i and one
# multiplicative perturbation per free parameter (x1.2; tau_c0 x1.5), five
# runs. Collected against the experiment markers; scripts/cp_gauss_newton.py
# turns the five curves into the next center (cal_center2.i / cal_gn2.i).
[Study]
  base       = cal_center1.i
  out_root   = ../../results/cp_study/cal/gn1
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 5
[]
[Cases]
  [./rho0_p]
    rho0 = 49.38
  []
  [./k1_p]
    k1 = 0.20124
  []
  [./y_c_p]
    y_c = 0.00384
  []
  [./tau_c0_p]
    tau_c0 = 8.025
  []
[]
