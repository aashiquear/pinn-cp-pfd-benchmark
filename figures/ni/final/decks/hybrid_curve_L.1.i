# PINN+FEM route, specimen L.1 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out            = ../hybrid_curve_L.1.png
  markers        = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_experiment_markers.csv
  experiment     = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/tl/cal/gn2/base/results.csv
  failure_target = 0.4633
  hide_experiment = false
  damage         = true
  d_markers      = false
  strain_max     = 70
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./run]
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/tl/final/hybrid/l.1
    label = 'specimen L.1, PINN+FEM'
    color = '#1f4e9c'
  []
[]
