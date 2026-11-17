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
[Extra]
  # the L.1 run of the W_c = 75 calibration (docs/12 §6.4; 30% path, newton_tol 1e-3): the
  # specimen the carried set was calibrated on (author 2026-10-28)
  [./L1]
    route    = fem
    specimen = L.1
    dir      = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
  []
[]
