# Reference run (no PINNACLE): nucleation and growth of polycrystalline spherulites, after
# Fig. 8 of (granasy 2005 spherulites): conditions of its Fig. 4(f) (S = 1.0, s0 = 0.1,
# chi_l = 54), no seed: crystals nucleate from the thermal noise. Paper: 5000 x 5000 after a
# pre-heat of 10,000 steps at 1595 K; here 1000 x 1000 from a homogeneous liquid (no pre-heat).
# Phase-field noise: standard deviation 0.007 per step (`calibrated`: a nucleation rate that gives about 50 crystals in this box;
# at the transformation time of the paper, 0.0062, only six nucleate. The paper prints 0.0025 without saying what it multiplies: with 0.0025 as a
# per-step standard deviation nothing nucleates in 7000 steps, with 0.007 the box is full by step 1100).
# Consumed by scripts/cryst_orientation_reference.py.
[Grid]
  n = 1000
  dx = 6.25e-3
  dt = 4.75e-6
  steps = 2600
  save_every = 10
  frame_every = 650
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
  chi_liquid = 54
  h_orientation = 1900
  noise_phi = 0.007
  noise_theta = 0.0375
  seed_radius_dx = 0
  rng_seed = 20261112
  out_dir = results/crystallization/orientation_reference/c5_n1000_noise0.007_seed2
[]
