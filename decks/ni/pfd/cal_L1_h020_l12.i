# l0/h = 6 point of the l0 plateau check on h020 (author 2026-10-20, option 1): l0 moved on the
# fixed mesh with G_c / l0 held at the calibrated l0 = 0.8 pair (40.75 / 0.8 = 50.9 MPa), so the
# material strength stays put and what changes is how well the mesh resolves the band. Plateau if
# curve and failure strain match l0/h = 4 (l0 0.8). One run; Gc_guess alone (first = '1').
[Calibrate]
  base      = base_L1_h020.i
  out_root  = ../../../results/ni_study/pfd/cal_L1_h020
  l0        = '1.2'
  Gc_guess  = '61.1'
  first     = '1'
  target    = 0.463
  tolerance = 0.004
  max_runs  = 1
  threads   = 2
[]
