# FEM route, specimen P.3 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out         = ../fem_shape_P.3.png
  panel       = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig5_panels/P.3_exp.png
  panel_label = 'P.3: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.18 0.76'
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
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/bunge/final/fem/p.3
    label = 'specimen P.3, FEM'
  []
[]
