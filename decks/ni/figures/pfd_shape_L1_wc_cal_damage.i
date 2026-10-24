# Calibrated W_c > 0 set (chi = 0, l0 0.6; docs/12 §6.4): fields at failure with the d = 0.85 outline
# and the crack line (max-d ridge across the gage width); color = damage. Read by
# scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_wc_cal_damage.png
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
  [./w55]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc55/l0.6/gc3.255
    label = '$W_c$ = 55 MPa, $\chi$ = 0'
  []
  [./w65]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc1.384
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
  []
  [./w75]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
  []
[]
