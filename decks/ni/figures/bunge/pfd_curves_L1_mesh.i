# After the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): Bunge angles, set refitted on L.1 (G_c 48.67).
# Restart of Part III (2026-11-15): mesh check of the single parameter set on L.1. h020 (792
# elements, the calibration run) against h010 (3,168 elements), same parameters.
# Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../../figures/ni/pfd/curves_L1_mesh.png
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
  [./h020]
    dir   = ../../../../results/ni_study/bunge/pfd/cal_L1_h020/l0.6/gc48.67
    label = 'h020, 792 elements'
    color = '#d1495b'
  []
  [./h010]
    dir   = ../../../../results/ni_study/bunge/pfd/mesh_L1_h010
    label = 'h010, 3,168 elements'
    color = '#00798c'
  []
[]
