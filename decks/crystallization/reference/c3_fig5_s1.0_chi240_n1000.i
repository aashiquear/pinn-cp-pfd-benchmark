# Reference run (no PINNACLE): binary + orientation-field model, after Fig. 5 of
# (granasy 2005 spherulites) at a quarter of its size: S = 1.0, s0 = 0.1, isotropic mobility,
# x = 0, chi_l = 240 (ASSUMED: the paper gives neither S nor chi_l for this figure).
# 1000 x 1000, dx = 6.25e-3.
# Consumed by scripts/cryst_orientation_reference.py. ASSUMED: h_orientation (see the script
# and plan section 8.1), the seed, the step count, the noise as per-step increments.
[Grid]
  n = 1000
  dx = 6.25e-3
  dt = 4.75e-6
  steps = 3200
  save_every = 200
  theta_substeps = 20
  precision = float64
[]
[Case]
  supersaturation = 1.0
  s0 = 0.1
  delta0 = 0.0
  k = 4
  x = 0.0
  m = 1
  n = 0.5
  chi_liquid = 240
  h_orientation = 1900
  noise_phi = 0.0025
  noise_theta = 0.0375
  seed_radius_dx = 12
  seed_theta = 0.0
  rng_seed = 20261117
  out_dir = results/crystallization/orientation_reference/fig5_s1.0_chi240_n1000
[]
