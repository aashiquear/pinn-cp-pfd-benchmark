# Reference run (no PINNACLE): binary + orientation-field model, conditions of Fig. 6 of
# (granasy 2005 spherulites): S = 0.75, two-fold kinetic coefficient with 99.5 % anisotropy
# (delta0 = 0.995, k = 2), isotropic interface free energy (s0 = 0), branching with a preferred
# misorientation of 30 degrees (m = 3, x = 0.15), seed at 30 degrees (theta = 1/6 of pi).
# 500 x 500, dx = 6.25e-3. Consumed by scripts/cryst_orientation_reference.py.
# ASSUMED: chi_l = 360 (not stated for this figure), n = 1/2, h_orientation, the seed, the
# step count, the noise as per-step increments.
[Grid]
  n = 500
  dx = 6.25e-3
  dt = 4.75e-6
  steps = 2600
  save_every = 130
  theta_substeps = 20
  precision = float64
[]
[Case]
  supersaturation = 0.75
  s0 = 0.0
  delta0 = 0.995
  k = 2
  x = 0.15
  m = 3
  n = 0.5
  chi_liquid = 360
  h_orientation = 1900
  noise_phi = 0.0025
  noise_theta = 0.0375
  seed_radius_dx = 12
  seed_theta = 0.1666667
  rng_seed = 20261117
  out_dir = results/crystallization/orientation_reference/fig6_s0.75
[]
