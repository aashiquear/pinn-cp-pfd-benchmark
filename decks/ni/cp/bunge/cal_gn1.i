# Ni CP calibration on L.1, Bunge convention (ADR 0036), Gauss-Newton iteration 1.
# The paper's set and one multiplicative perturbation (x1.2) per free parameter.
[Study]
  base       = cal_center1.i
  out_root   = ../../../../results/ni_study/bunge/cal/gn1
  markers    = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 4
[]
[Cases]
  [./rho0_p]
    rho0 = 2.82e6
  []
  [./k1_p]
    k1 = 0.0354
  []
  [./y_c_p]
    y_c = 0.948e-6
  []
[]
