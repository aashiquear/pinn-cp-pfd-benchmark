# Reference table: steady planar crystallization front of the single-field model in Ni
# (case A2). Consumed by scripts/cryst_reference.py (no PINNACLE run): travelling-wave
# solution by shooting on the speed; closed form v = -6 M Delta_f delta for p = phi^2 (3 - 2 phi).
[Material]
  # provenance: reported_table (Table 2 of granasy 2019 nucleation)
  molar_mass = 58.6934e-3   # kg/mol
  H_f        = 17.29e3      # J/mol
  T_f        = 1728.3       # K
  density    = 8357         # kg/m^3
  gamma_eq   = 0.302        # J/m^2
  d_10_90    = 1.0e-9       # m
[]
[Study]
  kind = front
  reduced_undercoolings = '0.10 0.20 0.30'
  interpolations = 'tanh standard'
  out_dir = results/crystallization/reference/ni_front
[]
