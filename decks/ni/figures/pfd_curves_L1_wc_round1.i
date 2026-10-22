# First round of the W_c stage (docs/12 §6.4): the six runs at l0 = 0.6 (W_c 55 / 65 / 75 MPa, two G_c
# each) against the experiment, the CP-only curve and the W_c = 0 carried pair; one panel per W_c,
# max d below. Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_wc_round1.png
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
  [./wc55_ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0 (carried pair)'
    color = '#9a9a9a'
    panel = 0
  []
  [./wc55_gc4]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc55/l0.6/gc4
    label = '$W_c$ = 55 MPa'
    color = '#d1495b'
    panel = 0
  []
  [./wc55_gc16]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc55/l0.6/gc16
    label = '$W_c$ = 55 MPa'
    color = '#1f4e9c'
    panel = 0
  []
  [./wc65_ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0 (carried pair)'
    color = '#9a9a9a'
    panel = 1
  []
  [./wc65_gc2]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc65/l0.6/gc2
    label = '$W_c$ = 65 MPa'
    color = '#d1495b'
    panel = 1
  []
  [./wc65_gc8]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc65/l0.6/gc8
    label = '$W_c$ = 65 MPa'
    color = '#1f4e9c'
    panel = 1
  []
  [./wc75_ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0 (carried pair)'
    color = '#9a9a9a'
    panel = 2
  []
  [./wc75_gc03]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc75/l0.6/gc0.3
    label = '$W_c$ = 75 MPa'
    color = '#d1495b'
    panel = 2
  []
  [./wc75_gc12]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/wc75/l0.6/gc1.2
    label = '$W_c$ = 75 MPa'
    color = '#1f4e9c'
    panel = 2
  []
[]
