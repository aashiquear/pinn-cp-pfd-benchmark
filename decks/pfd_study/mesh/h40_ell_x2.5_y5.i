# Stress-concentrator shape study (docs/09 §6, author's go-ahead 2026-10-14),
# level h40: series A/C: 2a = 5 um transverse (y), 2b = 2.5 um along the load; K_t 5. The hexaform plate generator's ellipse hole: diameter =
# 'x-extent y-extent' at angle 0 (the first axis rotates from +x by angle).
[Mesh]
  out_dir    = ../../../results/pfd_study/mesh/h40_ell_x2.5_y5
  width      = 100
  height     = 50
  thickness  = 1
  h          = 2.5
  hole       = ellipse
  diameter   = '2.5 5'
  angle      = 0
  ring_cells = 2
  n_r        = auto
  grading    = auto
  points     = ../../../results/cp_study/cal/gn2/base/orientations_points.csv
  figure     = ../../../figures/pfd/mesh_h40_ell_x2.5_y5.png
[]
