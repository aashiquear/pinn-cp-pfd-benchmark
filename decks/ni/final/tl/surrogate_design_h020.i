# Surrogate of the Ni specimens, restart of Part III (2026-11-15): training set of 30
# Latin-hypercube synthetic single crystals (Euler angles uniform over the rotation group; gage
# thickness as mean and linear taper over the five stations), run on the FEM route with the single
# parameter set and the total-Lagrangian balance (fem_base_h020.i) on h020 meshes. Same seed and
# ranges as the withdrawn first attempt (decks/ni/final/surrogate_design_h020.i), so the design
# rows are the same specimens. The 15 real specimens are held out.
# Read by scripts/ni_surrogate_design.py.
[Design]
  out_root   = ../../../../results/ni_study/tl/final/surrogate/data
  samples    = 30
  seed       = 0
  t_mean     = '0.325 0.400'
  t_taper    = '-0.04 0.06'
  mesh_deck  = ../../mesh/h020.i
  study_base = fem_base_h020.i
  concurrent = 7
[]
