# Field maps of the shape study (docs/10 §9): damage and equivalent plastic
# strain at fixed strains, the crack path with the d >= 0.8 elements solid,
# and the initiation/propagation GIF, circle first. Run after
# scripts/pfd_study.py decks/pfd_study/cp/shape_h40.i has collected.
[Fields]
  figures    = ../../../figures/pfd
  prefix     = cp_shape_h40
  strains    = '10 15 17.2 19'
  threshold  = 0.8
  field      = equivalent_plastic_strain
  vmax       = 0.5
  deform     = 1.0
  gif        = true
  gif_every  = 2
[]
[Runs]
  [./circle]
    dir   = ../../../results/pfd_study/cp/neural_h40/base
    label = 'circle 5 um'
  []
  [./A_x2.5_y5]
    dir   = ../../../results/pfd_study/cp/shape_h40/A_x2.5_y5
    label = 'A/C: 5 x 2.5 um'
  []
  [./A_x1.25_y5]
    dir   = ../../../results/pfd_study/cp/shape_h40/A_x1.25_y5
    label = 'A: 5 x 1.25 um'
  []
  [./B_x3.54_y7.07]
    dir   = ../../../results/pfd_study/cp/shape_h40/B_x3.54_y7.07
    label = 'B: 7.07 x 3.54 um'
  []
  [./B_x2.5_y10]
    dir   = ../../../results/pfd_study/cp/shape_h40/B_x2.5_y10
    label = 'B: 10 x 2.5 um'
  []
  [./C_a45]
    dir   = ../../../results/pfd_study/cp/shape_h40/C_a45
    label = 'C: 45 deg'
  []
  [./C_a90]
    dir   = ../../../results/pfd_study/cp/shape_h40/C_a90
    label = 'C: 90 deg'
  []
[]
