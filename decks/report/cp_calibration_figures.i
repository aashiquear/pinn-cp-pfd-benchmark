# Report figures for the CP FE calibration (docs/reports/cp_fe_calibration_report.md).
# Read by scripts/report_figures.py; every curve, marker file and axis limit the
# figures use is declared here, nothing is hard-wired in the script.
#
# Style (author's decision 2026-10-10): every stress-strain curve starts at the
# initial condition (0, 0), which is prepended; parameters are written as
# mathematical symbols; no grid lines; the full axes box is drawn; the box ends
# one tick block above the largest computed stress (700 MPa here) and at 20.5%
# strain. Engineering stress on both sides.

[Figures]
  out_dir     = ../../figures/report
  markers     = ../../data/digitized/al2219_fig3_experiment_markers.csv
  strain_max  = 20.5
  stress_max  = 700
  stress_tick = 100
  dpi         = 200
[]

[Calibrated]
  # Final calibrated set at h40 (results/cp_study/cal/gn2/base, deck decks/cp_study/cal_center2.i)
  results     = ../../results/cp_study/cal/gn2/base/results.csv
  grains      = ../../results/cp_study/cal/gn2/base/grains.npz
  grains_csv  = ../../results/cp_study/cal/gn2/base/grains.csv
  gauss_newton = ../../results/cp_study/cal/gn2/gauss_newton.json
  rho0        = 41.15
  k1          = 0.1677
  y_c         = 0.0032
  tau_c0      = 5.35
[]

# One-at-a-time perturbations around the calibrated set (same study, cal_gn2.i).
# symbol is a matplotlib mathtext string; factor is the multiplier on the center.
[Perturbations]
  [./rho0]
    results = ../../results/cp_study/cal/gn2/rho0_p/results.csv
    symbol  = $\rho_0$
    factor  = 1.2
    value   = 49.38
    unit    = $\mu$m$^{-2}$
  []
  [./k1]
    results = ../../results/cp_study/cal/gn2/k1_p/results.csv
    symbol  = $k_1$
    factor  = 1.2
    value   = 0.2013
    unit    = '-'
  []
  [./y_c]
    results = ../../results/cp_study/cal/gn2/y_c_p/results.csv
    symbol  = $y_c$
    factor  = 1.2
    value   = 0.00384
    unit    = $\mu$m
  []
  [./tau_c0]
    results = ../../results/cp_study/cal/gn2/tau_c0_p/results.csv
    symbol  = $\tau_{c0}$
    factor  = 1.5
    value   = 8.02
    unit    = MPa
  []
[]

# Mesh sensitivity with the calibrated set (decks/cp_study/mesh_h{20,80,160}.i;
# h40 is the calibrated run itself). label is the legend text; h is the
# element size in micrometers.
[Mesh]
  [./h20]
    results = ../../results/cp_study/mesh/h20/results.csv
    label   = h20
    h       = 5.0
  []
  [./h40]
    results = ../../results/cp_study/cal/gn2/base/results.csv
    label   = h40
    h       = 2.5
  []
  [./h80]
    results = ../../results/cp_study/mesh/h80/results.csv
    label   = h80
    h       = 1.25
  []
  [./h160]
    results = ../../results/cp_study/mesh/h160/results.csv
    label   = h160
    h       = 0.625
  []
[]
