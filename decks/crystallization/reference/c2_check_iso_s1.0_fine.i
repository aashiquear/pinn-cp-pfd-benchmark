# Reference check (no PINNACLE): binary model, single crystal, isotropic, S = 1.0, small grid.
# Consumed by scripts/cryst_binary_reference.py. Purpose: seed size and stability before the
# 500 x 500 cases of (granasy 2005 spherulites) Fig. 4.
[Grid]
  n = 400
  dx = 3.125e-3
  dt = 1.1875e-6
  steps = 2000
  save_every = 400
[]
[Case]
  supersaturation = 1.0
  s0 = 0.0
  k = 4
  seed_radius_dx = 24
  noise_phi = 0.0
  out_dir = results/crystallization/binary_reference/check_iso_s1.0_fine
[]
