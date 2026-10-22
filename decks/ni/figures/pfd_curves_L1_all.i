# Every run of the W_c = 0 G_c(l0) calibration on L.1 h020 (docs/12 §6.2), one panel per l0, with
# the max-d row: how the curve and the damage history move with G_c at each length scale. Read by
# scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_all.png
  markers        = ../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.463
  damage         = true
  panel_title    = '$\ell_0$ = 0.4 mm ($\ell_0/h$ = 2) | $\ell_0$ = 0.6 mm ($\ell_0/h$ = 3) | $\ell_0$ = 0.8 mm ($\ell_0/h$ = 4)'
  strain_max     = 60
  stress_max     = 300
  stress_tick    = 50
  legend         = below
  dpi            = 200
[]
[Runs]
  [./l04_gc55]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc5.5
    panel = 0
  []
  [./l04_gc11]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc11
    panel = 0
  []
  [./l04_gc1556]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc15.56
    panel = 0
  []
  [./l04_gc185]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc18.5
    panel = 0
  []
  [./l04_gc2017]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc20.17
    panel = 0
  []
  [./l04_gc211]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc21.1
    panel = 0
  []
  [./l04_gc22]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc22
    panel = 0
  []
  [./l06_gc8]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc8
    panel = 1
  []
  [./l06_gc2786]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc27.86
    panel = 1
  []
  [./l06_gc3095]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95
    panel = 1
  []
  [./l06_gc32]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc32
    panel = 1
  []
  [./l08_gc11]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc11
    panel = 2
  []
  [./l08_gc3818]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc38.18
    panel = 2
  []
  [./l08_gc4075]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc40.75
    panel = 2
  []
  [./l08_gc44]
    dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.8/gc44
    panel = 2
  []
[]
