# FEM route, specimen 9.4 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out            = ../fem_curve_9.4.png
  markers        = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_experiment_markers.csv
  experiment     = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.626
  hide_experiment = true
  damage         = true
  strain_max     = 80
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./run]
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/final/fem/9.4
    label = 'specimen 9.4, FEM'
    color = '#1f4e9c'
  []
[]
