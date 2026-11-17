# Result files shipped with the repository (author's decision 2026-10-15,
# reverting the "no field files in git" policy for the files the report
# presents): the field series behind the report figures and animations, the
# meshes and the trained networks. Read by scripts/ship_results.py.
[Ship]
  manifest = ../../results/SHIPPED.md
[]
[Series]
  # Part I (report §5, §8): plastic-strain fields at 2, 8, 15 and 20%
  [./part1_cal]
    dir = ../../results/cp_study/cal/gn2/base
    frames = '2 8 15 20'
  []
  [./part1_mesh]
    dir = ../../results/cp_study/mesh/*
    frames = '2 8 15 20'
  []
  [./part1_neural]
    dir = ../../results/cp_study/hybrid/h40_*
    frames = '2 8 15 20'
  []
  # Part II (report §15, §16): the two calibrated plates animated, the shape study every second snapshot
  [./part2_cal_cp]
    dir = ../../results/pfd_study/cp/neural_h40/base
    frames = all
  []
  [./part2_cal_j2]
    dir = ../../results/pfd_study/j2/sens_h40_wc200/base
    frames = all
  []
  [./part2_shape]
    dir = ../../results/pfd_study/cp/shape_h40/*
    frames = every:2
  []
[]
[Files]
  [./meshes]
    glob = ../../results/pfd_study/mesh/*/mesh.e
  []
  [./mesh_vtk]
    glob = ../../results/pfd_study/mesh/*/mesh.vtk
  []
  [./surrogate_fields]
    glob = ../../results/pfd_study/cp/surrogate/h40_field/fields.vtk
  []
  [./surrogate_nets]
    glob = ../../results/*/surrogate/h40_*/surrogate_network.npz
  []
  [./surrogate_data]
    glob = ../../results/pfd_study/cp/surrogate/h40_data/samples*/*/target.npz
  []
  [./part1_surr_data]
    glob = ../../results/cp_study/surrogate/h40_data/samples*/*/target.npz
  []
[]
