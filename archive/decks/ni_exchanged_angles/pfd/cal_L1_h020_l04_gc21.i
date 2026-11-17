# One further run for l0 = 0.4 (2026-10-20): the search of cal_L1_h020_resume_l04.i spent its
# six runs without a hit. Bracket: gc20.17 lost its load at 0.440 (stress 215 -> 141 MPa in one
# step; its mechanics failed the next step, before d > 0.85 spanned the section), gc22 had not
# separated by 0.477 (max d 0.87; OOM-killed). Interpolated in ln G_c between them: 21.1.
# Same out_root; first = '1' runs Gc_guess alone.
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/cal_L1_h020
  l0        = '0.4'
  Gc_guess  = '21.1'
  first     = '1'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 1
  threads   = 2
[]
