# Parameter sensitivity of the two neural routes carried forward (author's
# decision 2026-10-12) against the FE sensitivity study on h40. Read by
# scripts/neural_sensitivity.py.
[Compare]
  fe_root = ../../results/cp_study/sens/h40
  fe_base = ../../results/cp_study/cal/gn2/base
  out_dir = ../../results/cp_study/sens_neural
  figures = ../../figures/report
  prefix  = sens_neural_h40
  markers = ../../data/digitized/al2219_fig3_experiment_markers.csv
[]
[Routes]
  [./linear]
    root  = ../../results/cp_study/sens_neural/h40_linear
    base  = ../../results/cp_study/hybrid/h40_linear
    label = 'PINN+FEM, linear-extrapolation predictor'
  []
  [./grain]
    root  = ../../results/cp_study/sens_neural/h40_grain_solver
    base  = ../../results/cp_study/hybrid/h40_fem_solver_grain
    label = 'PINN as solver, grain-aware inputs (verified)'
  []
[]
