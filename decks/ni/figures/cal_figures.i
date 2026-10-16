# Figures of the Ni CP calibration on specimen L.1 (docs/12 §2) for the Part III
# report: the full set (curve, error, one-at-a-time panels and overlay,
# sensitivity) for the calibrated iteration named by `full_set`, and the
# convergence figure over all iterations (author, 2026-10-16: the earlier
# iterations appear only in the convergence figure). Read by
# scripts/ni_calibration_figures.py; every run, marker file and axis limit is
# declared here.
#
# Style (author's rules for the Al2219 report, 2026-10-10): every curve starts at
# (0, 0); parameters as mathematical symbols; no grid; closed axes box from
# (0, 0), one tick block above the largest plotted stress; engineering stress on
# both sides (Ni: reaction / (gage width x mean station thickness), against the
# gauge strain of the gage- / gage+ node sets, as the paper's DIC gauge).
# Curves end at the last step before the first load cutback: beyond 54-60%
# gauge strain the global Newton fails and the remaining steps are not a
# solution (docs/12 §2.4).

[Figures]
  out_dir     = ../../../figures/ni/cal
  markers     = ../../../data/digitized/ni_fig3_experiment_markers.csv
  strain_max  = 60
  stress_max  = 300
  stress_tick = 50
  dpi         = 200
  full_set    = gn4
[]

# Iterations: root holds base/ and <param>_p/ case directories and the
# gauss_newton.json the step script wrote (centre, per-marker errors, Jacobian).
[Iterations]
  [./gn1]
    root  = ../../../results/ni_study/cal/gn1
    label = 'iteration 1 (paper Table 3 set)'
  []
  [./gn2]
    root  = ../../../results/ni_study/cal/gn2
    label = 'iteration 2'
  []
  [./gn3]
    root  = ../../../results/ni_study/cal/gn3
    label = 'iteration 3'
  []
  [./gn4]
    root  = ../../../results/ni_study/cal/gn4
    label = 'iteration 4 (calibrated set)'
  []
[]

# Parameter symbols and units (mm, MPa, s); factor is the perturbation the
# study used around each centre (decks/ni/cp/cal_gn*.i).
[Parameters]
  [./rho0]
    symbol = $\rho_0$
    unit   = mm$^{-2}$
    factor = 1.2
  []
  [./k1]
    symbol = $k_m$
    unit   = '-'
    factor = 1.2
  []
  [./y_c]
    symbol = $y_c$
    unit   = mm
    factor = 1.2
  []
  [./tau_c0]
    symbol = $\tau_{c0}$
    unit   = MPa
    factor = 1.5
  []
[]
