# Fracture shape of the calibrated Ni CP+PFD runs on L.1 (W_c = 0; failure within 0.004 of 0.463)
# against Fig. 5 (docs/12 §6.2), read by scripts/ni_pfd_shape_figures.py; color = exx.
# l0 = 0.4 joins when its G_c = 21.1 run (decks/ni/pfd/cal_L1_h020_l04_gc21.i) ends.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_cal.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1.png
  panel_label = 'L.1: experiment (DIC, top) and the CP-PFD of Rezwan 2026 (bottom), Fig. 5'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  color       = exx
  dpi         = 200
[]
[Runs]
  [./l06]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$\ell_0$ = 0.6 mm'
  []
  [./l08]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc40.75
    label = '$\ell_0$ = 0.8 mm'
  []
[]
