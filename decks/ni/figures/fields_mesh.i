# Equivalent plastic strain fields of the Ni mesh study on L.1 (docs/12 §3), read by
# scripts/ni_fields.py. Plan view (x-y) of the deformed specimen around the gage, one
# row per run, one column per time; the colour scale is shared within a column. Only
# finished runs carry a field series (h010 and h020z2 were lost with the host).
[Figures]
  out      = ../../../figures/ni/mesh/mesh_eqps.png
  field    = equivalent_plastic_strain
  label    = 'equivalent plastic strain'
  times    = '2.0 4.0 5.6'          # s; bar strain = 0.05 x t
  x_range  = '5.0 20.4'             # mm, the gage, fillets and grip ends
  dpi      = 200
[]
[Runs]
  [./h040]
    pvd     = ../../../results/ni_study/cal/gn4/base/solver/fields.pvd
    results = ../../../results/ni_study/cal/gn4/base/results.csv
    label   = 'h040 (calibration mesh)'
  []
  [./h020]
    pvd     = ../../../results/ni_study/mesh_study/L1/h020/solver/fields.pvd
    results = ../../../results/ni_study/mesh_study/L1/h020/results.csv
    label   = 'h020'
  []
[]
