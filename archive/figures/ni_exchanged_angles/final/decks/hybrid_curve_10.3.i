# PINN+FEM route, specimen 10.3 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out            = ../hybrid_curve_10.3.png
  markers        = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_experiment_markers.csv
  experiment     = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/tl/cal/gn2/base/results.csv
  failure_target = 0.5284
  hide_experiment = true
  damage         = true
  d_markers      = false
  strain_max     = 70
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./run]
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/tl/final/hybrid/10.3
    label = 'specimen 10.3, PINN+FEM'
    color = '#1f4e9c'
  []
[]
