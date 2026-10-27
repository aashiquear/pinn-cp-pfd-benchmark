# Last stage, FEM route, specimen 10.3 (first finished; docs/12 §8): stress-strain curve and max d
# with the experiment's failure strain of 10.3 (Fig. 4(b): 0.528) as the dotted line. The Fig. 3
# markers are specimen L's and are left out. Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/final/fem_curve_10.3.png
  markers        = ../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.528
  hide_experiment = true
  damage         = true
  strain_max     = 80
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./s103]
    dir   = ../../../results/ni_study/final/fem/10.3
    label = 'specimen 10.3, FEM'
    color = '#2a9d8f'
  []
[]
