# After the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): Bunge angles, set refitted on L.1 (G_c 48.67).
# Restart of Part III (2026-11-15): fracture shape of the calibrated L.1 run (G_c 48.67 MPa mm,
# l0 0.6 mm, failure strain 0.4635 against 0.463) against Fig. 5; color = exx.
# Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../../figures/ni/pfd/shape_L1_cal.png
  panel       = ../../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'
  color       = exx
  dpi         = 200
[]
[Runs]
  [./gc48]
    dir   = ../../../../results/ni_study/bunge/pfd/cal_L1_h020/l0.6/gc48.67
    label = '$\ell_0$ = 0.6 mm'
  []
[]
