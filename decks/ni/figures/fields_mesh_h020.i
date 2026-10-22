# Equivalent plastic strain fields of the Ni mesh study on L.1 with the h020-calibrated set
# (docs/12 §5), read by scripts/ni_fields.py. Plan view (x-y) of the deformed specimen
# around the gage, one row per run, one column per time; the color scale is shared
# within a column. PINNACLE writes the snapshots (every 20 steps) at the end of a run and
# records them only at converged steps: h010's are valid to 4.8 s (step 240; its curve
# ends at step 242); h020z2 was stopped by the memory cap and has none.
[Figures]
  out      = ../../../figures/ni/mesh_h020/mesh_eqps.png
  field    = equivalent_plastic_strain
  label    = 'equivalent plastic strain'
  times    = '1.6 3.2 4.8'          # s; bar strain = 0.05 x t (8 / 16 / 24%)
  x_range  = '5.0 20.4'             # mm, the gage, fillets and grip ends
  dpi      = 200
[]
[Runs]
  [./h040]
    pvd     = ../../../results/ni_study/mesh_study/L1_h020set/h040/solver/fields.pvd
    results = ../../../results/ni_study/mesh_study/L1_h020set/h040/results.csv
    label   = 'h040 (reference)'
  []
  [./h020]
    pvd     = ../../../results/ni_study/cal_h020/gn2/base/solver/fields.pvd
    results = ../../../results/ni_study/cal_h020/gn2/base/results.csv
    label   = 'h020 (production mesh)'
  []
  [./h010]
    pvd     = ../../../results/ni_study/mesh_study/L1_h020set/h010/solver/fields.pvd
    results = ../../../results/ni_study/mesh_study/L1_h020set/h010/results.csv
    label   = 'h010'
  []
[]
