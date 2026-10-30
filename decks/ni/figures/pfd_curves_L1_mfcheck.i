# Multi-fidelity check (docs/12 §8.5, docs/13 P7.4): L.1 with the carried CP+PFD set on h040 against
# the h020 calibration run; experiment to the adopted failure strain, CP-only h020 curve; max d below.
# Read by scripts/ni_pfd_curve_figures.py.
[Figures]
  out            = ../../../figures/ni/pfd/curves_L1_mfcheck.png
  markers        = ../../../data/digitized/ni_fig3_experiment_markers.csv
  experiment     = ../../../data/digitized/ni_fig3_specimenL_calibration.csv
  cp_only        = ../../../results/ni_study/cal_h020/gn2/base/results.csv
  failure_target = 0.463
  damage         = true
  panel_title    = 'L.1, carried set ($\chi$ = 0, $W_c$ = 75 MPa, $G_c$ = 0.1 MPa mm, $l_0$ = 0.6 mm)'
  legend         = below
  strain_max     = 60
  stress_max     = 300
  stress_tick    = 50
  dpi            = 200
[]
[Runs]
  [./h020]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = 'h020 (calibration run)'
    color = '#1f4e9c'
    panel = 0
  []
  [./h040]
    dir   = ../../../results/ni_study/pfd/mfcheck_L1_h040
    label = 'h040 (same set)'
    color = '#d1495b'
    panel = 0
  []
[]
