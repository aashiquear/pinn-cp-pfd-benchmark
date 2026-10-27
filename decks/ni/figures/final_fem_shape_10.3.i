# Last stage, FEM route, specimen 10.3: E_xx at failure with the 0.85 outline and the crack line,
# against the experiment's DIC panel of 10.3. Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/final/fem_shape_10.3.png
  panel       = ../../../data/digitized/ni_fig5_panels/10.3_exp.png
  panel_label = '10.3: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.16 1.27'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'
  color       = exx
  dpi         = 200
[]
[Runs]
  [./s103]
    dir   = ../../../results/ni_study/final/fem/10.3
    label = 'specimen 10.3, FEM'
  []
[]
