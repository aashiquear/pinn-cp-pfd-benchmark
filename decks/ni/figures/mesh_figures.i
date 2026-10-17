# Figures of the Ni mesh study on specimen L.1 (docs/12 §3) with the calibrated set.
# Read by scripts/ni_mesh_figures.py. Same style as the calibration figures
# (decks/ni/figures/cal_figures.i). Curves end at the first load cutback; a run
# the host lost is shown up to its last accepted step (`partial` in the label).
[Figures]
  out_dir     = ../../../figures/ni/mesh
  markers     = ../../../data/digitized/ni_fig3_experiment_markers.csv
  strain_max  = 60
  stress_max  = 300
  stress_tick = 50
  dpi         = 200
  strains     = '10 20 30 40'      # gauge strains (%) of the convergence figure
[]
# h is the gage cell (mm); layers the elements through the thickness.
[Mesh]
  [./h040]
    results = ../../../results/ni_study/mesh_study/L1/base/results.csv
    label   = 'h040 (calibration mesh)'
    h       = 0.4
    layers  = 1
  []
  [./h020]
    results = ../../../results/ni_study/mesh_study/L1/h020/results.csv
    label   = 'h020'
    h       = 0.2
    layers  = 1
  []
  [./h010]
    results = ../../../results/ni_study/mesh_study/L1/h010/results.csv
    label   = 'h010 (partial: host lost at step 112)'
    h       = 0.1
    layers  = 1
  []
  [./h020z2]
    results = ../../../results/ni_study/mesh_study/L1/h020z2/results.csv
    label   = 'h020, two layers'
    h       = 0.2
    layers  = 2
  []
[]
