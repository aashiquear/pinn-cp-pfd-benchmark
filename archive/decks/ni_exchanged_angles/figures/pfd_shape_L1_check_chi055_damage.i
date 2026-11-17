# Quick check (author 2026-10-23): the phase field at the end of the path for W_c 75 / G_c 30 / chi 0.55
# (no failure), with the carried pair for comparison. Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_check_chi055_damage.png
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
  [./carried]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$W_c$ = 75 MPa, $\chi$ = 0, carried'
  []
  [./check]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/check_wc75_gc30_chi0.55_fix/l0.6/gc30
    label = '$W_c$ = 75 MPa, $\chi$ = 0.55'
  []
[]
