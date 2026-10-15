# Shipped result files

Field series, meshes and networks behind the report's figures and tables,
force-added to git by `scripts/ship_results.py decks/report/ship_results.i`
(the .gitignore keeps every other field file out). Series list the
snapshots the figures read; `all` series feed an animation.

| set | run | files | size (MB) |
|---|---|---|---|
| part1_cal: 2 8 15 20 | `results/cp_study/cal/gn2/base` | fields.pvd, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 1.9 |
| part1_mesh: 2 8 15 20 | `results/cp_study/mesh/h160` | fields.pvd, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 29.7 |
| part1_mesh: 2 8 15 20 | `results/cp_study/mesh/h20` | fields.pvd, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 0.5 |
| part1_mesh: 2 8 15 20 | `results/cp_study/mesh/h80` | fields.pvd, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 7.4 |
| part1_neural: 2 8 15 20 | `results/cp_study/hybrid/h40_fem_residual` | fields.pvd, fields.vtk, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 2.7 |
| part1_neural: 2 8 15 20 | `results/cp_study/hybrid/h40_fem_residual_grain` | fields.pvd, fields.vtk, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 2.6 |
| part1_neural: 2 8 15 20 | `results/cp_study/hybrid/h40_fem_solver` | fields.pvd, fields.vtk, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 2.6 |
| part1_neural: 2 8 15 20 | `results/cp_study/hybrid/h40_fem_solver_grain` | fields.pvd, fields.vtk, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 2.6 |
| part1_neural: 2 8 15 20 | `results/cp_study/hybrid/h40_linear` | fields.pvd, fields_step0005.vtk, fields_step0025.vtk, fields_step0047.vtk, fields_step0063.vtk | 2.2 |
| part2_cal_cp: all | `results/pfd_study/cp/neural_h40/base` | 58 files (fields.pvd .. fields_step0055.vtk) | 31.7 |
| part2_cal_j2: all | `results/pfd_study/j2/sens_h40_wc200/base` | 60 files (fields.pvd .. fields_step0057.vtk) | 30.6 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/a_x1.25_y5` | 32 files (fields.pvd .. fields_step0058.vtk) | 17.4 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/a_x2.5_y5` | 32 files (fields.pvd .. fields_step0057.vtk) | 17.4 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/b_x2.5_y10` | 28 files (fields.pvd .. fields_step0049.vtk) | 15.2 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/b_x3.54_y7.07` | 29 files (fields.pvd .. fields_step0052.vtk) | 16.3 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/base` | 31 files (fields.pvd .. fields_step0055.vtk) | 16.8 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/c_a45` | 33 files (fields.pvd .. fields_step0059.vtk) | 18.4 |
| part2_shape: every:2 | `results/pfd_study/cp/shape_h40/c_a90` | 33 files (fields.pvd .. fields_step0060.vtk) | 17.9 |
| meshes | `.` | 10 files (mesh.e .. mesh.e) | 2.7 |
| mesh_vtk | `.` | 10 files (mesh.vtk .. mesh.vtk) | 1.7 |
| surrogate_fields | `.` | fields.vtk | 0.3 |
| surrogate_nets | `.` | surrogate_network.npz, surrogate_network.npz, surrogate_network.npz | 0.2 |
| surrogate_data | `.` | 48 files (target.npz .. target.npz) | 0.5 |

Total: 239 MB in 457 files.
