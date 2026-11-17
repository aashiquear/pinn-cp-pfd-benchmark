# Second round of the W_c stage, chi = 0 (docs/12 §6.4): the finished runs at l0 = 0.6 against the
# experiment, the CP-only curve and the W_c = 0 carried pair; one panel per W_c, max d below.
# Read by scripts/ni_pfd_curve_figures.py. Runs join as they finish.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_wc_chi0.png
  markers        = ../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.463
  damage         = true
  panel_title    = '$W_c$ = 55 MPa | $W_c$ = 65 MPa | $W_c$ = 75 MPa'
  legend         = below
  strain_max     = 60
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./w55_ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0, $\chi$ = 1 (carried pair)'
    color = '#9a9a9a'
    panel = 0
  []
  [./w55_gc1333]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc55/l0.6/gc1.333
    label = '$W_c$ = 55 MPa, $\chi$ = 0'
    color = '#d1495b'
    panel = 0
  []
  [./w55_gc4]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc55/l0.6/gc4
    label = '$W_c$ = 55 MPa, $\chi$ = 0'
    color = '#1f4e9c'
    panel = 0
  []
  [./w55_gc16]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc55/l0.6/gc16
    label = '$W_c$ = 55 MPa, $\chi$ = 0'
    color = '#00798c'
    panel = 0
  []
  [./w65_ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0, $\chi$ = 1 (carried pair)'
    color = '#9a9a9a'
    panel = 1
  []
  [./w65_gc06667]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc0.6667
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
    color = '#d1495b'
    panel = 1
  []
  [./w65_gc1384]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc1.384
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
    color = '#1f4e9c'
    panel = 1
  []
  [./w65_gc2]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc2
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
    color = '#00798c'
    panel = 1
  []
  [./w65_gc8]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc8
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
    color = '#edae49'
    panel = 1
  []
  [./w75_ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0, $\chi$ = 1 (carried pair)'
    color = '#9a9a9a'
    panel = 2
  []
  [./w75_gc12]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc1.2
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
    color = '#d1495b'
    panel = 2
  []
  [./w75_gc03]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.3
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
    color = '#1f4e9c'
    panel = 2
  []
  [./w75_gc01]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
    color = '#00798c'
    panel = 2
  []
  [./w75_gc003333]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.03333
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
    color = '#edae49'
    panel = 2
  []
[]
