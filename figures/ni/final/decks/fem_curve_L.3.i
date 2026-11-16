# FEM route, specimen L.3 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out            = ../fem_curve_L.3.png
  markers        = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_experiment_markers.csv
  experiment     = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/bunge/cal/gn2/center_a/results.csv
  failure_target = 0.4814
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
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/bunge/final/fem/l.3
    label = 'specimen L.3, FEM'
    color = '#1f4e9c'
  []
[]
