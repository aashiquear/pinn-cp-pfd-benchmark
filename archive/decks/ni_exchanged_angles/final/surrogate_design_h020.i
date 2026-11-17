# Surrogate of the Ni specimens (docs/12 §8; author 2026-10-25): training set of 30 Latin-hypercube
# synthetic single crystals (Euler angles uniform over the rotation group; gage thickness as mean
# and linear taper over the five stations), run with the carried CP+PFD set on h020-style meshes;
# the 15 real specimens are held out. Read by scripts/ni_surrogate_design.py.
[Design]
  out_root   = ../../../results/ni_study/final/surrogate/data
  samples    = 30
  seed       = 0
  t_mean     = '0.325 0.400'
  t_taper    = '-0.04 0.06'
  mesh_deck  = ../mesh/h020.i
  study_base = fem_base_h020.i
  concurrent = 6
[]
