# Reference run (no PINNACLE): binary + orientation-field model, conditions of Fig. 4(c) of
# (granasy 2005 spherulites): S = 0.8, four-fold interface energy s0 = 0.1, isotropic mobility,
# single-minimum orientational energy (x = 0, n = 1/2), orientational mobility chi_l = 54
# (360 = the paper's value; 54 = reduced by the factor 0.15). 500 x 500, dx = 6.25e-3.
# Consumed by scripts/cryst_orientation_reference.py. ASSUMED: h_orientation (see the script
# and plan section 8.1), the seed, the step count, the noise as per-step increments.
[Grid]
  n = 500
  dx = 6.25e-3
  dt = 4.75e-6
  steps = 2000
  save_every = 100
  theta_substeps = 20
  precision = float64
[]
[Case]
  supersaturation = 0.8
  s0 = 0.1
  delta0 = 0.0
  k = 4
  x = 0.0
  m = 1
  n = 0.5
  chi_liquid = 54
  h_orientation = 1900
  noise_phi = 0.0025
  noise_theta = 0.0375
  seed_radius_dx = 12
  seed_theta = 0.0
  rng_seed = 20261117
  out_dir = results/crystallization/orientation_reference/fig4c_s0.8_chi54
[]
