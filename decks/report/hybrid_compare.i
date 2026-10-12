# Hybrid first check at h40: the cold FE run against the linear and network
# load-step predictors (decks/cp_study/hybrid/h40_*.i). Read by
# scripts/hybrid_compare.py; the first run is the reference.
[Compare]
  out_dir = ../../results/cp_study/hybrid
  figures = ../../figures/report
  prefix  = hybrid_h40
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
  [./physics]
    dir   = ../../results/cp_study/hybrid/h40_physics
    label = 'physics-trained network predictor (torch GPU)'
  []
[]
