# chi sweep at W_c 75 / G_c 0.1 / l0 0.6 (docs/12 §6.4, author 2026-10-23): does the elastic
# contribution move the crack? Read by scripts/ni_pfd_curve_figures.py; chi 0.8 joins when done.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_chi_sweep.png
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
  [./chi0]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = '$\chi$ = 0'
    color = '#00798c'
  []
  [./chi03]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi_sweep_wc75/chi0.3/l0.6/gc0.1
    label = '$\chi$ = 0.3'
    color = '#d1495b'
  []
  [./chi06]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi_sweep_wc75/chi0.6/l0.6/gc0.1
    label = '$\chi$ = 0.6'
    color = '#edae49'
  []
[]
