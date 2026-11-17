# Length-scale sensitivity at fixed Gc on the h40 plate (docs/09 §5.1 step 3,
# the Fig. 8 analogue of chakraborty 2016 multiscale): the calibrated base
# (l0 = 5 um, its calibrated Gc, filled in from results/pfd_study/cp/cal_h40/
# calibrated.json by scripts/pfd_fill_base.py) with l0 swept at 2h, 4h, 8h
# and 16h of the h40 background (h = 2.5 um), everything else unchanged.
# The base case is the calibrated run itself. Seven cases fit the host side
# by side (2 threads each).
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/cp/sweep_l_h40
  concurrent = 3
[]
[Sweep]
  l0 = '2.5 10 20 40'
[]
