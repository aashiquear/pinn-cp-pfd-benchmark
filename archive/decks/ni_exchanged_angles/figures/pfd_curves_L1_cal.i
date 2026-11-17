# Stress-strain curves of the calibrated Ni CP+PFD runs on L.1, W_c = 0 (docs/12 §6.2): the runs of
# the shape figure decks/ni/figures/pfd_shape_L1_cal.i, against the Fig. 3 experiment and the
# calibrated CP-only curve. Read by scripts/ni_pfd_curve_figures.py. l0 = 0.4 (G_c 21.1) and the
# l0/h = 6 check (l0 1.2, G_c 61.1) join when their runs end.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_cal.png
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
  [./l06]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$\ell_0$ = 0.6 mm'
    color = '#edae49'
  []
  [./l08]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc40.75
    label = '$\ell_0$ = 0.8 mm'
    color = '#00798c'
  []
[]
