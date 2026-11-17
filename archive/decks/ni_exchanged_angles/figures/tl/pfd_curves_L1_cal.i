# Restart of Part III (2026-11-15), total-Lagrangian balance: the G_c search on L.1 h020 with the
# paper's fracture model (AT2, chi 1, zeta 1, no threshold, l0 0.6 mm) against the Fig. 3
# experiment and the CP-only curve of the carried set. Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../../figures/ni/pfd/curves_L1_cal.png
  markers        = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../../results/ni_study/tl/cal/gn2/base/results.csv
  failure_target = 0.463
  damage         = true
  strain_max     = 60
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./gc30]
    dir   = ../../../../results/ni_study/tl/pfd/cal_L1_h020/l0.6/gc30
    color = '#edae49'
  []
  [./gc49]
    dir   = ../../../../results/ni_study/tl/pfd/cal_L1_h020/l0.6/gc49.37
    label = 'calibrated'
    color = '#d1495b'
  []
  [./gc60]
    dir   = ../../../../results/ni_study/tl/pfd/cal_L1_h020/l0.6/gc60
    color = '#00798c'
  []
  [./gc120]
    dir   = ../../../../results/ni_study/tl/pfd/cal_L1_h020/l0.6/gc120
    color = '#1f4e9c'
  []
[]
