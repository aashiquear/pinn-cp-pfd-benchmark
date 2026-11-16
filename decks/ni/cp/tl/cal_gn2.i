# Ni CP calibration on the total-Lagrangian balance, Gauss-Newton iteration 2.
[Study]
  base       = cal_center2.i
  out_root   = ../../../../results/ni_study/tl/cal/gn2
  markers    = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 4
[]
[Cases]
  [./rho0_p]
    rho0 = 2.43084e+06
  []
  [./k1_p]
    k1 = 0.029205
  []
  [./y_c_p]
    y_c = 5.60268e-07
  []
[]
