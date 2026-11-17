# Restart of Part III (2026-11-15): the fifteen specimens on the FEM route with the single
# parameter set (decks/ni/final/tl/fem_base_h020.i): strain-to-failure parity and fracture
# location against the experiment. Read by scripts/ni_final_figures.py. Experiment columns only.
[Figures]
  out_dir  = ../../../figures/ni/final
  prefix   = fem
  parity   = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  fracture = ../../../data/digitized/ni_fig5_fracture_descriptors.csv
  dpi      = 200
[]
[Routes]
  [./fem]
    root  = ../../../results/ni_study/tl/final/fem
    label = 'FEM'
  []
[]
