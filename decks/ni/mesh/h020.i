# Ni single-crystal specimens (Rezwan 2026), dog-bone meshes at gage cell 0.2 mm:
# planform from data/extracted/ni_fig2_planform.csv, thickness stations from
# Table 2, grips at the common nominal 0.35 mm (author's instruction of
# 2026-10-15: identical ends, measured variation in the 8 mm gage).
[Mesh]
  out_root        = ../../../results/ni_study/mesh/h020
  planform        = ../../../data/extracted/ni_fig2_planform.csv
  thickness       = ../../../data/extracted/ni_table2_thickness.csv
  specimens       = all
  station_spacing = 2.0
  grip_thickness  = 0.35
  h               = 0.2
  h_grip          = 0.4
  n_z             = 1
  figure          = ../../../figures/ni/mesh_h020.png
[]
