# Restart of Part III (2026-11-15): the specimen surrogate of decks/ni/final/tl/surrogate/surrogate_h020.i.
# Held-out synthetic specimens against their FE runs (curves and scalars), and the fifteen real
# specimens (never training data) against the FE route and the experiment.
# Read by scripts/ni_surrogate_figures.py. Experiment columns only.
[Figures]
  out_dir  = ../../../figures/ni/final
  prefix   = surrogate
  model    = ../../../results/ni_study/tl/final/surrogate/model_h020
  targets  = ../../../results/ni_study/tl/final/surrogate/targets_h020
  fem      = ../../../results/ni_study/tl/final/fem
  parity   = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  fracture = ../../../data/digitized/ni_fig5_fracture_descriptors.csv
  dpi      = 200
[]
