# Rerun after the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16): Bunge angles, set refitted on L.1.
# Surrogate of the Ni specimens, restart of Part III (2026-11-15): training set of 30
# Latin-hypercube synthetic single crystals (Euler angles uniform over the rotation group; gage
# thickness as mean and linear taper over the five stations), run on the FEM route with the single
# parameter set and the total-Lagrangian balance (fem_base_h020.i) on h020 meshes. Seed as before;
# the thickness ranges are widened to contain the fifteen real specimens (measured means 0.322 to
# 0.395 mm, tapers -0.073 to +0.067 mm; the first design's -0.04 to +0.06 left four outside).
# The 15 real specimens are held out.
# Read by scripts/ni_surrogate_design.py.
[Design]
  out_root   = ../../../../results/ni_study/bunge/final/surrogate/data
  samples    = 30
  seed       = 0
  t_mean     = '0.320 0.400'
  t_taper    = '-0.08 0.07'
  mesh_deck  = ../../mesh/h020.i
  study_base = fem_base_h020.i
  concurrent = 7
[]
