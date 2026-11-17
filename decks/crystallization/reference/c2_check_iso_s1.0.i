# Reference check (no PINNACLE): binary model, single crystal, isotropic, S = 1.0, small grid.
# Consumed by scripts/cryst_binary_reference.py. Purpose: seed size and stability before the
# 500 x 500 cases of (granasy 2005 spherulites) Fig. 4.
[Grid]
  n = 200
  dx = 6.25e-3
  dt = 4.75e-6
  steps = 3000
  save_every = 500
[]
[Case]
  supersaturation = 1.0
  s0 = 0.0
  k = 4
  seed_radius_dx = 12
  noise_phi = 0.0
  out_dir = results/crystallization/binary_reference/check_iso_s1.0
[]
