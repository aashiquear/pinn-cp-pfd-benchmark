# After the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): the G_c search on L.1 h020 with the
# paper's fracture model (AT2, chi 1, zeta 1, no threshold, l0 0.6 mm) against the Fig. 3
# experiment and the CP-only curve of the carried set. Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../../figures/ni/pfd/curves_L1_cal.png
  markers        = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../../results/ni_study/bunge/cal/gn2/base/results.csv
  failure_target = 0.463
  damage         = true
  strain_max     = 60
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./gc42]
    dir   = ../../../../results/ni_study/bunge/pfd/cal_L1_h020/l0.6/gc42.5
    color = '#edae49'
  []
  [./gc48]
    dir   = ../../../../results/ni_study/bunge/pfd/cal_L1_h020/l0.6/gc48.67
    label = 'calibrated'
    color = '#d1495b'
  []
  [./gc57]
    dir   = ../../../../results/ni_study/bunge/pfd/cal_L1_h020/l0.6/gc57.5
    color = '#00798c'
  []
[]
