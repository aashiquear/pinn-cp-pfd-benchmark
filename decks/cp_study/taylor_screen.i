# Taylor material-point SURROGATE for the CP calibration (not the answer: an
# isostrain aggregate is an upper bound on the FE polycrystal; it locates the
# region and gives sensitivities cheaply, the FE at h40 decides).
#
# Bounds come from the experiment's own features, Rezwan (2025) Fig. 3
# markers (Vaughan 2024): yield ~345 MPa at 0.56% sets rho0; the ~600 MPa
# plateau sets k1/(2 y_c); the hardening between 2% and 7% separates k1 from
# y_c; tau_c0 is a uniform offset. Fixed (decision 2026-10-09): C11 C12 C44
# b mu q alpha0, and gamma0, n unless convergence forces a change.
[Surrogate]
  markers      = ../../data/digitized/al2219_fig3_experiment_markers.csv
  out_dir      = ../../results/cp_study/taylor_screen
  n_grains     = 30
  seed         = 20260910
  strain_rate  = 1.0e-3
  n_increments = 200
  workers      = 12
  stress       = engineering
[]
[Fixed]
  C11 = 108.0e3
  C12 = 61.3e3
  C44 = 28.5e3
  mu  = 25.7e3
  gamma0 = 1.0e-5
  rate_exponent_n = 8.333333333333334
  q_lat = 1.4
  alpha0 = 0.4
  burgers = 0.000286
[]
[Free]
  # name = 'lower upper start'
  rho0   = '20 80 40'
  k1     = '0.02 0.12 0.05'
  annihilation_diameter = '0.001 0.006 0.00268'
  tau_c0 = '0 50 5'
[]
