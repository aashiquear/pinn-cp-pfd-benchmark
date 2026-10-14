# PFD benchmark mesh, level h160 (0.625 um background, two refinements below the CP study's
# h40 mesh), 100 x 50 x 1 um plate, central circular hole of diameter 5 um
# (author's choice 2026-10-13: the Irwin plastic zone of the hole stays clear
# of the plate edges). One HEX8 through the thickness; conforming O-grid ring
# around the hole (hexaform `plate`, see docs/workflow_pfd_plate.md). Grain
# ids in the figure come from the CP study's 30-grain point file, assigned
# at element centroids exactly as PINNACLE's point_file source does.
[Mesh]
  out_dir    = ../../../results/pfd_study/mesh/h160_hole5
  width      = 100
  height     = 50
  thickness  = 1
  h          = 0.625
  hole       = circle
  diameter   = 5
  ring_cells = 2
  n_r        = auto
  grading    = auto
  points     = ../../../results/cp_study/cal/gn2/base/orientations_points.csv
  figure     = ../../../figures/pfd/mesh_h160_hole5.png
[]
