# After the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): Bunge angles, set refitted on L.1 (G_c 48.67).
# Restart of Part III (2026-11-15): per-specimen curve and field figures of the PINN+FEM route with the
# single parameter set, on the experiment panel's own color scale. Read by
# scripts/ni_final_specimen_figures.py.
[Figures]
  root      = ../../../results/ni_study/bunge/final/hybrid
  out_dir   = ../../../figures/ni/final
  prefix    = hybrid
  label     = 'PINN+FEM'
  parity    = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv
  scales    = ../../../data/digitized/ni_fig5_panel_scales.csv
  panels    = ../../../data/digitized/ni_fig5_panels
  cp_only   = ../../../results/ni_study/bunge/cal/gn2/base/results.csv
  specimens = all
  ridge_always = true
  strain_max = 70
  dpi       = 200
[]
