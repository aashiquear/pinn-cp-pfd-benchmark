# After the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): Bunge angles, set refitted on L.1 (G_c 48.67).
# Restart of Part III (2026-11-15): L.1 on the three solver routes with the single parameter set
# (G_c 48.67 MPa mm, l0 0.6 mm): FEM, PINN+FEM (network start, Newton corrector) and PINN
# (network field accepted when its residual is within 5 % of the reaction).
# Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../../figures/ni/final/routes_curves_L1.png
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
  [./fem]
    dir   = ../../../../results/ni_study/bunge/final/fem/l.1
    label = 'FEM'
    color = '#1f4e9c'
  []
  [./hybrid]
    dir   = ../../../../results/ni_study/bunge/final/hybrid/l.1
    label = 'PINN+FEM'
    color = '#d1495b'
  []
  [./pinn]
    dir   = ../../../../results/ni_study/bunge/final/pinn/L.1
    label = 'PINN'
    color = '#00798c'
  []
[]
