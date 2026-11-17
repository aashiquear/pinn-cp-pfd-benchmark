# Recheck of 2026-11-16 (recalibrated set, Bunge angles, total-Lagrangian balance, implicit update).
# CP+PFD plate surrogates against the solver routes (docs/10 §7.4): held-out
# accuracy, data and training cost, time per query, break-even against the
# FE, linear and verified-solver routes of decks/pfd_study/cp/bunge/neural_h40.i.
# Read by scripts/surrogate_compare.py after the surrogate decks have run.
[Compare]
  out_dir  = ../../../../../results/pfd_study/cp/bunge/surrogate
  figures  = ../../../../../figures/pfd
  prefix   = cp_surrogate_h40
  design   = ../../../../../results/pfd_study/cp/bunge/surrogate/h40_data
  routes   = ../../../../../results/pfd_study/cp/bunge/neural_h40/neural_h40_summary.csv
  mesh_vtk = ../../../../../results/pfd_study/cp/bunge/surrogate/h40_field/fields.vtk
[]
[Surrogates]
  [./curve]
    dir   = ../../../../../results/pfd_study/cp/bunge/surrogate/h40_curve
    label = 'curve surrogate (stress and max damage vs strain), CPU'
  []
  [./field]
    dir   = ../../../../../results/pfd_study/cp/bunge/surrogate/h40_field
    label = 'damage-field surrogate at the end of the path, CPU'
  []
[]
