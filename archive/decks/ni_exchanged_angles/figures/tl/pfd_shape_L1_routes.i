# Restart of Part III (2026-11-15): fracture field of L.1 on the three solver routes with the
# single parameter set, against Fig. 5; color = exx, d = 0.85 outline and ridge line.
# Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../../figures/ni/final/routes_shape_L1.png
  panel       = ../../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'
  color       = exx
  dpi         = 200
[]
[Runs]
  [./fem]
    dir   = ../../../../results/ni_study/tl/final/fem/l.1
    label = 'FEM'
  []
  [./hybrid]
    dir   = ../../../../results/ni_study/tl/final/hybrid/l.1
    label = 'PINN+FEM'
  []
  [./pinn]
    dir   = ../../../../results/ni_study/tl/final/pinn/L.1
    label = 'PINN'
  []
[]
