# PINN+FEM route, specimen 10.1 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out         = ../hybrid_shape_10.1.png
  panel       = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig5_panels/10.1_exp.png
  panel_label = '10.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.16 1.28'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_always = True
  ridge_x     = '8.7 16.7'
  color       = exx
  dpi         = 200
[]
[Runs]
  [./run]
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/tl/final/hybrid/10.1
    label = 'specimen 10.1, PINN+FEM'
  []
[]
