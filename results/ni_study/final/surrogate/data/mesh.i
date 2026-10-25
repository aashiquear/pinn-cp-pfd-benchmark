# Meshes of the surrogate design's synthetic specimens (written by scripts/ni_surrogate_design.py).
[Mesh]
  out_root        = ./mesh
  planform        = ../../../../../data/extracted/ni_fig2_planform.csv
  thickness       = ./thickness.csv
  specimens       = all
  station_spacing = 2.0
  grip_thickness  = 0.35
  h               = 0.2
  h_grip          = 0.4
  n_z             = 1
  figure          = ./mesh.png
[]
