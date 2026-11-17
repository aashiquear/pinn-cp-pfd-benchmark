# Figures of the Ni mesh study on specimen L.1 with the h020-calibrated set (docs/12 §5;
# author 2026-10-17: h040 shown as the reference, the finer meshes converge on h020).
# Read by scripts/ni_mesh_figures.py. Same style as the calibration figures
# (decks/ni/figures/cal_figures.i). Every run postdates PINNACLE ddc76c3, so curves run
# over every converged step up to the first unconverged one (ni_specimen.py curve); a
# run that spent its cutback budget early ends there (labeled with its end strain).
[Figures]
  out_dir     = ../../../figures/ni/mesh_h020
  markers     = ../../../data/digitized/ni_fig3_experiment_markers.csv
  strain_max  = 60
  stress_max  = 300
  stress_tick = 50
  dpi         = 200
  strains     = '10 20 30 40'      # gauge strains (%) of the convergence figure
  cut_at_first_cutback = false
[]
# h is the gage cell (mm); layers the elements through the thickness.
[Mesh]
  [./h040]
    results = ../../../results/ni_study/mesh_study/L1_h020set/h040/results.csv
    label   = 'h040 (reference)'
    h       = 0.4
    layers  = 1
  []
  [./h020]
    results = ../../../results/ni_study/cal_h020/gn2/base/results.csv
    label   = 'h020 (production mesh, calibrated)'
    h       = 0.2
    layers  = 1
  []
  [./h010]
    results = ../../../results/ni_study/mesh_study/L1_h020set/h010/results.csv
    label   = 'h010'
    h       = 0.1
    layers  = 1
  []
  [./h020z2]
    results = ../../../results/ni_study/mesh_study/L1_h020set/h020z2/results.csv
    label   = 'h020, two layers'
    h       = 0.2
    layers  = 2
  []
[]
