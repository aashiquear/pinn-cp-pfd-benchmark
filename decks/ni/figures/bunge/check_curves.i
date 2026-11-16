# Check after the Euler-convention fix (PINNACLE ADR 0036): L.1 (calibration), 9.1 and H.1 with the
# single set refitted on L.1 (G_c 48.67 MPa mm), against the same specimens run with the exchanged
# angles (restart set, G_c 49.37). Experiment: Fig. 3 trace of L, up to its failure strain.
# Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../../figures/ni/bunge/check_curves.png
  markers        = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../../results/ni_study/bunge/cal/gn2/center_a/results.csv
  failure_target = 0.463
  damage         = true
  strain_max     = 70
  stress_max     = 300
  stress_tick    = 50
  legend         = below
  dpi            = 200
[]
[Runs]
  [./l1]
    dir   = ../../../../results/ni_study/bunge/pfd/cal_L1_h020/l0.6/gc48.67
    label = 'L.1, Bunge angles (experiment fails at 0.463)'
    color = '#d1495b'
    panel = 0
  []
  [./n91]
    dir   = ../../../../results/ni_study/bunge/final/check/fem_9.1_H.1/9.1
    label = '9.1, Bunge angles (experiment 0.546)'
    color = '#1f4e9c'
    panel = 0
  []
  [./h1]
    dir   = ../../../../results/ni_study/bunge/final/check/fem_9.1_H.1/h.1
    label = 'H.1, Bunge angles (experiment 0.549)'
    color = '#8d5a97'
    panel = 0
  []
  [./ol1]
    dir   = ../../../../results/ni_study/tl/final/fem/l.1
    label = 'L.1, exchanged angles'
    color = '#d1495b'
    panel = 1
  []
  [./o91]
    dir   = ../../../../results/ni_study/tl/final/fem/9.1
    label = '9.1, exchanged angles'
    color = '#1f4e9c'
    panel = 1
  []
  [./oh1]
    dir   = ../../../../results/ni_study/tl/final/fem/h.1
    label = 'H.1, exchanged angles'
    color = '#8d5a97'
    panel = 1
  []
[]
