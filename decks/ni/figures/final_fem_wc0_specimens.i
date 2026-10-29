# Last stage, FEM route with W_c = 0, chi = 0, G_c 30.95 (docs/12 §8.3): per-specimen curve and
# field figures, one pair per finished specimen, each on the experiment panel's own color scale.
# The crack outline and ridge are drawn only once the failure definition is met (author 2026-10-28).
# Read by scripts/ni_final_specimen_figures.py.
[Figures]
  root      = ../../../results/ni_study/final/fem_wc0
  out_dir   = ../../../figures/ni/final
  prefix    = fem_wc0
  label     = 'FEM, W_c = 0'
  parity    = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  scales    = ../../../data/digitized/ni_fig5_panel_scales.csv
  panels    = ../../../data/digitized/ni_fig5_panels
  cp_only   = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  specimens = all
  dpi       = 200
[]
