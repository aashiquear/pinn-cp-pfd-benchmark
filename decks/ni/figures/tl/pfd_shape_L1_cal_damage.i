# Restart of Part III (2026-11-15): fracture shape of the calibrated L.1 run (G_c 49.37 MPa mm,
# l0 0.6 mm, failure strain 0.4635 against 0.463) against Fig. 5; color = damage.
# Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../../figures/ni/pfd/shape_L1_cal_damage.png
  panel       = ../../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'
  color       = damage
  dpi         = 200
[]
[Runs]
  [./gc49]
    dir   = ../../../../results/ni_study/tl/pfd/cal_L1_h020/l0.6/gc49.37
    label = '$\ell_0$ = 0.6 mm'
  []
[]
