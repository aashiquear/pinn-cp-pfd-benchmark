# Parameter sensitivity of the calibrated J2+PFD plate at h40 (docs/09 §5.1
# step 4), one parameter at a time: Gc x0.5 and x2, l0 x0.5 and x2 (both at
# the calibrated Gc), the fixed viscosity eta switched on at two levels to
# show it does not matter at the quasi-static rate, the plastic-work
# threshold W_c at two levels, and the staggered alternation count 5 against
# 3. The base is the calibrated run. Values marked @ are filled from
# calibrated.json by scripts/pfd_fill_base.py.
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/j2/sens_h40_wc20
  concurrent = 7
[]
[Sweep]
  Gc  = '1001.53 4006.13'
  l0  = '2.5 10'
  eta = '1.0 10.0'
  W_c = '0.0 5.0'
  max_iterations = '5'
[]
