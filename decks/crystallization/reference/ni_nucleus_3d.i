# Reference table: critical nucleus of the single-field model in undercooled Ni
# (granasy 2019 nucleation, Fig. 3 and Fig. 5; material data of its Table 2).
# Consumed by scripts/cryst_reference.py (no PINNACLE run): shooting solution of Eq. 15.
[Material]
  # provenance: reported_table (Table 2 of the review)
  molar_mass = 58.6934e-3   # kg/mol
  H_f        = 17.29e3      # J/mol, heat of fusion
  T_f        = 1728.3       # K
  density    = 8357         # kg/m^3
  gamma_eq   = 0.302        # J/m^2, solid-liquid interface free energy at T_f
  d_10_90    = 1.0e-9       # m, 10-90 % interface thickness
[]
[Study]
  dim = 3
  # Fig. 3: equidistant reduced undercoolings between 0.15 and 0.30 (seven solid lines);
  # 0.10 and 0.125 added for the low end of Fig. 3(c)
  reduced_undercoolings = '0.10 0.125 0.15 0.175 0.20 0.225 0.25 0.275 0.30'
  # well:interpolation pairs (Table 1 of the review); Fig. 5 legend:
  # standard PF, g_fcc & p_PF, GL g_bcc & p_bcc, GL g_fcc & p_fcc
  functions = 'quartic:standard gl_fcc:standard quartic:gl_bcc gl_fcc:gl_fcc'
  # proportional: eps^2 and w proportional to T (gamma = gamma_eq T/T_f, constant width);
  # constant: eps^2 and w fixed at their melting-point values. Which one Fig. 3 used is not
  # stated for the PF line; both are tabulated and compared with the digitized curve.
  temperature_scaling = 'proportional constant'
  out_dir = results/crystallization/reference/ni_nucleus_3d
[]
