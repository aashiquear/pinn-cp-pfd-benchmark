# CP+PFD mesh sensitivity at fixed l0 = 5 um and the calibrated Gc = 1106
# MPa um (docs/09 §5.1 step 3, CP): the h40 base against the h80 and h160
# hole meshes (l0/h = 2, 4, 8), same load path and staggered settings.
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/cp/sweep_h_l5
  concurrent = 2
[]
[Cases]
  [./h80]
    mesh = ../../../results/pfd_study/mesh/h80_hole5/mesh.e
  []
  [./h160]
    mesh = ../../../results/pfd_study/mesh/h160_hole5/mesh.e
  []
[]
