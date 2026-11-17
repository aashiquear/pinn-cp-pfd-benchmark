# FE samples for the plate surrogate (ADR 0030 data source), h40, the
# calibrated structure and orientations: 24 Latin-hypercube rows over the four
# calibrated parameters, each spanning a factor 4 in log scale (x0.5 .. x2 of
# the calibrated value, the range of the sensitivity study), run through the
# study pipeline four at a time. Read by scripts/surrogate_design.py.
[Design]
  base       = ../cal_center2.i
  out_root   = ../../../results/cp_study/surrogate/h40_data
  samples    = 24
  sampler    = lhs
  seed       = 0
  concurrent = 4
  columns    = eng_stress
  abscissa   = nominal_strain
  fields     = 'displacement equivalent_plastic_strain'
[]
[Parameters]
  [./rho0]
    path  = mechanics/rho0
    range = '23.9936 95.9744'
    scale = log
  []
  [./k1]
    path  = mechanics/k1
    range = '0.0994165 0.397666'
    scale = log
  []
  [./y_c]
    path  = mechanics/y_c
    range = '0.00171879 0.00687516'
    scale = log
  []
  [./tau_c0]
    path  = mechanics/tau_c0
    range = '2.64574 10.5829'
    scale = log
  []
[]
