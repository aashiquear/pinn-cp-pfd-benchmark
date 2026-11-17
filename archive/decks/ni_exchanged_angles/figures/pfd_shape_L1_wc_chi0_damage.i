# Second round of the W_c stage, chi = 0 (docs/12 §6.4): fields at failure (separated) or at the last
# snapshot (so marked); color = damage. Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_wc_chi0_damage.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'          # the gage (undeformed x), where the rows are read
  color       = damage
  dpi         = 200
[]
[Runs]
  [./w55_gc4]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc55/l0.6/gc4
    label = '$W_c$ = 55 MPa, $\chi$ = 0'
  []
  [./w65_gc1384]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc1.384
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
  []
  [./w75_gc01]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
  []
  [./w75_gc003333]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.03333
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
  []
[]
