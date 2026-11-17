# Plate surrogate against the solver routes (report §11): held-out accuracy,
# data and training cost, time per query, break-even. Read by
# scripts/surrogate_compare.py after the surrogate decks have run on the host.
[Compare]
  out_dir = ../../results/cp_study/surrogate
  figures = ../../figures/report
  prefix  = surrogate_h40
  design  = ../../results/cp_study/surrogate/h40_data
  markers = ../../data/digitized/al2219_fig3_experiment_markers.csv
  routes  = ../../results/cp_study/hybrid/hybrid_h40_summary.csv
[]
[Surrogates]
  [./curve]
    dir   = ../../results/cp_study/surrogate/h40_curve
    label = 'curve surrogate, trained and queried on the GPU'
  []
  [./curve_cpu]
    dir   = ../../results/cp_study/surrogate/h40_curve_cpu
    label = 'curve surrogate, CPU'
  []
  [./field]
    dir   = ../../results/cp_study/surrogate/h40_field
    label = 'displacement-field surrogate at 20%, GPU'
  []
[]
