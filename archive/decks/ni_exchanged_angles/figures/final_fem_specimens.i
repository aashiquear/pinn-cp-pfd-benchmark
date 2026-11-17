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
[Extra]
  # the L.1 run of the W_c = 75 calibration (docs/12 §6.4; 30% path, newton_tol 1e-3): the
  # specimen the carried set was calibrated on (author 2026-10-28)
  [./L1]
    specimen = L.1
    dir      = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
  []
[]
