# Mesh sensitivity at fixed l0 = 5 um and the calibrated Gc (docs/09 §5.1
# step 3): the h40 base against the h80 and h160 hole meshes (l0/h = 2, 4, 8),
# same load path, same staggered settings. The brittle rule l0 > 2h is the
# coarsest point. Three cases, 2 threads each.
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/j2/sweep_h_l5
  concurrent = 3
[]
[Cases]
  [./h80]
    mesh = ../../../results/pfd_study/mesh/h80_hole5/mesh.e
  []
  [./h160]
    mesh = ../../../results/pfd_study/mesh/h160_hole5/mesh.e
  []
[]
