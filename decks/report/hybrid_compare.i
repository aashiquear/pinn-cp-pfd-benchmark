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
    dir   = ../../results/cp_study/hybrid/h40_cold
    label = 'cold FE (previous step as initial guess)'
  []
  [./linear]
    dir   = ../../results/cp_study/hybrid/h40_linear
    label = 'linear extrapolation predictor'
  []
  [./network]
    dir   = ../../results/cp_study/hybrid/h40_network
    label = 'network predictor (torch GPU)'
  []
[]
