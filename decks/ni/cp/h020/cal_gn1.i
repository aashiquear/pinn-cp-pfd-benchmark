# Ni CP calibration on the h020 mesh, Gauss-Newton iteration 1: the h040-calibrated set as
# centre and one multiplicative perturbation per free parameter (scripts/ni_gauss_newton.py).
# The centre run is the mesh study's h020 case (same deck values), copied in as `base`.
[Study]
  base       = cal_center1.i
  out_root   = ../../../../results/ni_study/cal_h020/gn1
  markers    = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 4
[]
[Cases]
  [./rho0_p]
    rho0 = 822925
  []
  [./k1_p]
    k1 = 0.013306
  []
  [./y_c_p]
    y_c = 2.71759e-06
  []
  [./tau_c0_p]
    tau_c0 = 0.368724
  []
[]
