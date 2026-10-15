# Ni CP calibration, Gauss-Newton iteration 1 (docs/11 §5 stage 2): the paper's
# Table 3 set as the centre (base_L1_h040.i) and one multiplicative perturbation
# per free parameter (x1.2; tau_c0 x1.5), for the finite-difference Jacobian of the
# engineering stress at the Fig. 3 experiment markers. scripts/ni_gauss_newton.py
# takes the Levenberg-Marquardt step and writes cal_center2.i / cal_gn2.i.
[Study]
  base       = base_L1_h040.i
  out_root   = ../../../results/ni_study/cal/gn1
  markers    = ../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 5
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
  [./tau_c0_p]
    tau_c0 = 0.39
  []
[]
