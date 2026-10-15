# Parameter sensitivity of the calibrated CP+PFD plate at h40 (docs/09 §5.1
# step 4), one parameter at a time: Gc x0.5 and x2, l0 x0.5 and x2 (both at
# the calibrated Gc), the fixed viscosity eta switched on at two levels to
# show it does not matter at the quasi-static rate, the plastic-work
# threshold W_c at two levels, and the staggered alternation count 5 against
# 3. The base is the calibrated run. Values marked @ are filled from
# calibrated.json by scripts/pfd_fill_base.py.
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/cp/sens_h40
  concurrent = 3
[]
[Sweep]
  Gc  = '276.53 553.059 2212.24 4424.48'   # x0.25, x0.5, x2, x4 of 1106.12 (author 2026-10-14: four per parameter)
  l0  = '2.5 3.75 7.5 10'
  eta = '1.0 3.0 10.0 30.0'
  W_c = '5.0 10.0 20.0 40.0'
  max_iterations = '5'
[]
