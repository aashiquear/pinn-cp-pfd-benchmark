# After the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): crystal plasticity on L.1 h020,
# no damage, gauge strain rate about 1e-3/s: the paper's Table 3 set and the carried set
# (decks/ni/cp/bunge/cal_center2a.i) against the Fig. 3 experiment of specimen L.
# Read by scripts/ni_pfd_curve_figures.py ([Curves]: runs without damage, drawn to their end).
[Figures]
  out            = ../../../../figures/ni/cal/cal_curve.png
  markers        = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../../data/digitized/ni_fig3_specimenL_calibration.csv
  failure_target = 0.463
  damage         = false
  strain_max     = 70
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Curves]
  [./paper]
    path  = ../../../../results/ni_study/bunge/cal/gn1/base/results.csv
    label = 'Table 3 of the paper, RMSE 14.3 MPa'
    color = '#00798c'
  []
  [./carried]
    path  = ../../../../results/ni_study/bunge/cal/gn2/center_a/results.csv
    label = 'carried set, RMSE 4.1 MPa'
    color = '#d1495b'
  []
[]
