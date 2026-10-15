# Figures of the stress-concentrator shape study (docs/10 §9): curves and
# panels per case (categorical), the circle as the base.
[Figures]
  study         = shape_h40.i
  figures       = ../../../figures/pfd
  prefix        = cp_shape_h40
  param         = case
  markers       = ../../../data/digitized/al2219_fig3_experiment_markers.csv
  target_strain = 0.172
[]
[Labels]
  base           = 'circle 5 um (K_t 3)'
  A_x2.5_y5      = 'A/C: 5 x 2.5 um, a = 2.5 (K_t 5)'
  A_x1.25_y5     = 'A: 5 x 1.25 um, a = 2.5 (K_t 9)'
  B_x3.54_y7.07  = 'B: 7.07 x 3.54 um, same area (K_t 5)'
  B_x2.5_y10     = 'B: 10 x 2.5 um, same area (K_t 9)'
  C_a45          = 'C: 5 x 2.5 um at 45 deg'
  C_a90          = 'C: 5 x 2.5 um at 90 deg (K_t 2)'
[]
