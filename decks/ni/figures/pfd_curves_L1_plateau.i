# l0 plateau on h020 at G_c / l0 held near 51 MPa (docs/12 §6.2, author's option 1, 2026-10-20): the
# calibrated pairs l0/h = 2, 3, 4 and the l0/h = 6 check. Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_plateau.png
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
  [./l04]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc21.1
    label = '$\ell_0$ = 0.4 mm ($\ell_0/h$ = 2)'
    color = '#d1495b'
  []
  [./l06]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    label = '$\ell_0$ = 0.6 mm ($\ell_0/h$ = 3)'
    color = '#edae49'
  []
  [./l08]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc40.75
    label = '$\ell_0$ = 0.8 mm ($\ell_0/h$ = 4)'
    color = '#00798c'
  []
  [./l12]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l1.2/gc61.1
    label = '$\ell_0$ = 1.2 mm ($\ell_0/h$ = 6)'
    color = '#1f4e9c'
  []
[]
