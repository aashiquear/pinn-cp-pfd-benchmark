# Fracture shape of the l0 plateau runs on L.1 (G_c / l0 near 51 MPa; l0/h = 2, 3, 4, 6) against
# Fig. 5 (docs/12 §6.2), read by scripts/ni_pfd_shape_figures.py; color = damage.
# Result (2026-10-21): the d > 0.85 zone widens with l0 (2.2 to 3.4 mm) and drifts along the gage (0.17 to 0.32).
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_plateau_damage.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1.png
  panel_label = 'L.1: experiment (DIC, top) and the CP-PFD of Rezwan 2026 (bottom), Fig. 5'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  color       = damage
  dpi         = 200
[]
[Runs]
  [./l04]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc21.1
    label = '$\ell_0$ = 0.4 mm, $\ell_0/h$ = 2'
  []
  [./l06]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$\ell_0$ = 0.6 mm, $\ell_0/h$ = 3'
  []
  [./l08]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc40.75
    label = '$\ell_0$ = 0.8 mm, $\ell_0/h$ = 4'
  []
  [./l12]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l1.2/gc61.1
    label = '$\ell_0$ = 1.2 mm, $\ell_0/h$ = 6'
  []
[]
