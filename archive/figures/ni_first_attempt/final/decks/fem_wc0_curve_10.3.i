# FEM, W_c = 0 route, specimen 10.3 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out            = ../fem_wc0_curve_10.3.png
  markers        = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_experiment_markers.csv
  experiment     = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.5284
  hide_experiment = true
  damage         = true
  strain_max     = 80
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./run]
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/final/fem_wc0/10.3
    label = 'specimen 10.3, FEM, W_c = 0'
    color = '#1f4e9c'
  []
[]
