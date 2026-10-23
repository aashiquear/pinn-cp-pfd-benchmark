# Calibrated W_c > 0 set (chi = 0, l0 0.6; docs/12 §6.4): the three calibrated pairs against the
# experiment, the CP-only curve and the W_c = 0 carried pair, max d below. Read by
# scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_wc_cal.png
  markers        = ../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.463
  damage         = true
  strain_max     = 60
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./ref]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0, $\chi$ = 1 (carried pair)'
    color = '#9a9a9a'
  []
  [./w55]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc55/l0.6/gc3.255
    label = '$W_c$ = 55 MPa, $\chi$ = 0'
    color = '#d1495b'
  []
  [./w65]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc65/l0.6/gc1.384
    label = '$W_c$ = 65 MPa, $\chi$ = 0'
    color = '#edae49'
  []
  [./w75]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$W_c$ = 75 MPa, $\chi$ = 0'
    color = '#00798c'
  []
[]
