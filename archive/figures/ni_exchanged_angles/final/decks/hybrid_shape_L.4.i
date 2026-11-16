# PINN+FEM route, specimen L.4 (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out         = ../hybrid_shape_L.4.png
  panel       = /mnt/extended/projects/pinn-cp-pfd-benchmark/data/digitized/ni_fig5_panels/L.4_exp.png
  panel_label = 'L.4: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.2 0.91'
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
    dir   = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/ni_study/tl/final/hybrid/l.4
    label = 'specimen L.4, PINN+FEM'
  []
[]
