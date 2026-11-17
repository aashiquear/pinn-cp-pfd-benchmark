# Last stage, FEM route with W_c = 0, chi = 0 and the G_c of the W_c = 0 calibration (G_c 30.95,
# l0 0.6; docs/12 §8.3): strain-to-failure parity and fracture location of the 15 specimens against
# the experiment. Read by scripts/ni_final_figures.py. Experiment columns only. No [Extra] block:
# the W_c = 0 calibration run of L.1 used chi = 1 and is not this set; L.1 is one of the 15.
[Figures]
  out_dir  = ../../../figures/ni/final
  prefix   = fem_wc0
  parity   = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  fracture = ../../../data/digitized/ni_fig5_fracture_descriptors.csv
  dpi      = 200
[]
[Routes]
  [./fem_wc0]
    root  = ../../../results/ni_study/final/fem_wc0
    label = 'FEM, W_c = 0, chi = 0'
  []
[]
