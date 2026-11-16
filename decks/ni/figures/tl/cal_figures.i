# Restart of Part III (2026-11-15): crystal plasticity on L.1 h020 with the total-Lagrangian
# balance, gauge strain rate about 1e-3/s. Iteration 1 is centered on the paper's Table 3,
# iteration 2 on the carried set. Read by scripts/ni_calibration_figures.py.
[Figures]
  out_dir     = ../../../../figures/ni/cal
  markers     = ../../../../data/digitized/ni_fig3_experiment_markers.csv
  strain_max  = 60
  stress_max  = 300
  stress_tick = 50
  dpi         = 200
  full_set    = gn2
  cut_at_first_cutback = false
[]
[Iterations]
  [./gn1]
    root  = ../../../../results/ni_study/tl/cal/gn1
    label = 'iteration 1 (Table 3 of the paper)'
  []
  [./gn2]
    root  = ../../../../results/ni_study/tl/cal/gn2
    label = 'iteration 2 (carried set)'
  []
[]
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
