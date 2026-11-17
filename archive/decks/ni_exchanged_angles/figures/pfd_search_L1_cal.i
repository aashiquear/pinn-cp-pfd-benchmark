# The W_c = 0 G_c searches of the Ni CP+PFD calibration on L.1 h020 (docs/12 §6.2) as one figure:
# failure strain against G_c per l0, and against G_c / l0 (the collapse). Read by
# scripts/ni_pfd_search_figures.py; every run under results/ni_study/pfd/cal_L1_h020 (void_* skipped).
[Figures]
  out            = ../../../figures/ni/pfd/search_L1_cal.png
  root           = ../../../results/ni_study/pfd/cal_L1_h020
  failure_target = 0.463
  tolerance      = 0.004
  h              = 0.2
  dpi            = 200
[]
