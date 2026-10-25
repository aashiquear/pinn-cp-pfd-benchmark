# Quick check, not for reporting (author 2026-10-23): W_c 75 / G_c 30 / chi 0.55 at l0 0.6, against the
# carried pair (chi 0 / W_c 75 / G_c 0.1) and the W_c = 0 pair. Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_check_chi055.png
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
  [./wc0]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$W_c$ = 0, $\chi$ = 1'
    color = '#9a9a9a'
  []
  [./carried]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$W_c$ = 75 MPa, $\chi$ = 0 (carried)'
    color = '#00798c'
  []
  [./check]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/check_wc75_gc30_chi0.55_fix/l0.6/gc30
    label = '$W_c$ = 75 MPa, $\chi$ = 0.55'
    color = '#d1495b'
  []
[]
