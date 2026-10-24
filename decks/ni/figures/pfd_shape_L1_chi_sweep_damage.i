# chi sweep at W_c 75 / G_c 0.1 (docs/12 §6.4): the phase field at failure for chi 0 (mid-gage) and
# chi 0.6 (fillet); chi 0.3 died before writing fields, chi 0.8 joins when done. Read by
# scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_chi_sweep_damage.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1_exp.png
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
  [./chi0]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$\chi$ = 0'
  []
  [./chi06]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi_sweep_wc75/chi0.6/l0.6/gc0.1
    label = '$\chi$ = 0.6'
  []
[]
