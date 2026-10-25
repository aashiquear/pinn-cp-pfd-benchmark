# Last stage, FEM route: strain-to-failure parity and fracture location of the 15 specimens against
# the experiment (docs/12 §8). Read by scripts/ni_final_figures.py. Experiment columns only.
[Figures]
  out_dir  = ../../../figures/ni/final
  prefix   = fem
  parity   = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  fracture = ../../../data/digitized/ni_fig5_fracture_descriptors.csv
  dpi      = 200
[]
[Routes]
  [./fem]
    root  = ../../../results/ni_study/final/fem
    label = 'FEM'
  []
[]
