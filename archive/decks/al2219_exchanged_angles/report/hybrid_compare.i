# Neural routes at h40: the cold FE run against the linear-extrapolation and
# FE-residual-trained network predictors and the verified PINN-as-solver run (decks/cp_study/hybrid/h40_*.i). Read by
# scripts/hybrid_compare.py; the first run is the reference.
[Compare]
  out_dir = ../../results/cp_study/hybrid
  figures = ../../figures/report
  prefix  = hybrid_h40
  corrector_tol = 0.05               # drawn on the residual-ratio figure
[]
[Runs]
  [./cold]
    dir   = ../../results/cp_study/cal/gn2/base
    label = 'cold FE (previous step as initial guess)'
  []
  [./linear]
    dir   = ../../results/cp_study/hybrid/h40_linear
    label = 'linear extrapolation predictor'
  []
  [./fem_residual]
    dir   = ../../results/cp_study/hybrid/h40_fem_residual
    label = 'FE-residual-trained network predictor (torch GPU)'
  []
  [./fem_solver]
    dir   = ../../results/cp_study/hybrid/h40_fem_solver
    label = 'PINN as solver, verified (corrector = auto)'
  []
  [./fem_residual_grain]
    dir   = ../../results/cp_study/hybrid/h40_fem_residual_grain
    label = 'FE-residual-trained predictor, grain-aware inputs'
  []
  [./fem_solver_grain]
    dir   = ../../results/cp_study/hybrid/h40_fem_solver_grain
    label = 'PINN as solver, verified, grain-aware inputs'
  []
[]
