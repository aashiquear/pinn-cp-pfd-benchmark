# Figures of the neural_h40 study (scripts/pfd_study_figures.py).
[Figures]
  study         = neural_h40.i
  figures       = ../../../figures/pfd
  prefix        = cp_neural_h40
  param         = case
  xlabel        = ''
  logx          = true
  markers       = ../../../data/digitized/al2219_fig3_experiment_markers.csv
  target_strain = 0.172
[]
[Labels]
  base = 'FE (cold Newton)'
  linear = 'PINN+FEM, linear extrapolation'
  solver_grain = 'PINN as solver, grain-aware, verified'
[]
