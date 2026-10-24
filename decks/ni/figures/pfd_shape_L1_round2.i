# Fracture shape of the Ni CP+PFD calibration runs on L.1 against Fig. 5 (docs/12 §6), read by
# scripts/ni_pfd_shape_figures.py. Interim set: the round-2 runs of decks/ni/pfd/cal_L1_h020.i
# (failure below the target at every l0), for the figure layout; the calibrated runs replace them.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_round2.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  dpi         = 200
[]
[Runs]
  [./l04]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc11
    label = '$\ell_0$ = 0.4 mm'
  []
  [./l06]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc27.86
    label = '$\ell_0$ = 0.6 mm'
  []
  [./l08]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc38.18
    label = '$\ell_0$ = 0.8 mm'
  []
[]
