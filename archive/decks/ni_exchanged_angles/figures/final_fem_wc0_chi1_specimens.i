# W_c = 0 baseline (chi 1, G_c 30.95, l0 0.6; author 2026-11-02): per-specimen curve and field
# figures, on the experiment panel's own color scale. The ridge is drawn on every run, failed or
# not; the curve carries the through-section d > 0.55 / 0.65 / 0.75 / 0.85 markers
# (ni_specimen.py thresholds). Read by scripts/ni_final_specimen_figures.py.
[Figures]
  root      = ../../../results/ni_study/final/fem_wc0_chi1
  out_dir   = ../../../figures/ni/final
  prefix    = fem_wc0_chi1
  label     = 'FEM, W_c = 0'
  parity    = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  scales    = ../../../data/digitized/ni_fig5_panel_scales.csv
  panels    = ../../../data/digitized/ni_fig5_panels
  cp_only   = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  specimens = all
  ridge_always = true
  d_markers = true
  strain_max = 100
  dpi       = 200
[]
