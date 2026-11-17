# Reference run (no PINNACLE): binary model of (granasy 2005 spherulites), single crystal,
# the conditions of Fig. 4(e): S = 1.0, interface-energy anisotropy s0 = 0.1 (k = 4), isotropic
# mobility, 500 x 500 grid, dx = 6.25e-3, dt = 4.75e-6 (6.6 um box). Consumed by
# scripts/cryst_binary_reference.py. Assumed (not in the paper): the seed (disc of 12 dx at the
# solidus composition), the number of steps, and the phase-field noise as a per-step
# increment of standard deviation 0.0025 p(phi) (liquid only).
[Grid]
  n = 500
  dx = 6.25e-3
  dt = 4.75e-6
  steps = 1800
  save_every = 100
[]
[Case]
  supersaturation = 1.0
  s0 = 0.1
  k = 4
  crystal_angle_deg = 0
  seed_radius_dx = 12
  noise_phi = 0.0025
  rng_seed = 20261117
  out_dir = results/crystallization/binary_reference/fig4e_aniso_s1.0
[]
