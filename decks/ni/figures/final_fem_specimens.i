# Last stage, FEM route: per-specimen curve and field figures (docs/12 §8), one pair per finished
# specimen, each on the experiment panel's own color scale. Read by scripts/ni_final_specimen_figures.py.
[Figures]
  root      = ../../../results/ni_study/final/fem
  out_dir   = ../../../figures/ni/final
  prefix    = fem
  label     = 'FEM'
  parity    = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  scales    = ../../../data/digitized/ni_fig5_panel_scales.csv
  panels    = ../../../data/digitized/ni_fig5_panels
  cp_only   = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  specimens = all
  dpi       = 200
[]
