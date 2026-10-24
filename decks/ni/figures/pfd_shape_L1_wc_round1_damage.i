# First round of the W_c stage (docs/12 §6.4): the field at the last snapshot of the runs that reached
# the gage (no separation: mechanics failure at d -> 1) and of one that cracked the fillet (a killed run has no snapshots); color = damage.
# Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_wc_round1_damage.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  color       = damage
  dpi         = 200
[]
[Runs]
  [./wc55_gc16]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc55/l0.6/gc16
    label = '$W_c$ = 55 MPa'
  []
  [./wc65_gc8]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc65/l0.6/gc8
    label = '$W_c$ = 65 MPa'
  []
  [./wc65_gc2]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc65/l0.6/gc2
    label = '$W_c$ = 65 MPa'
  []
[]
