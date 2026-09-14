# Circular-grain shrinkage, the quantitative calibration of the kinetics.
# A circular grain of one order parameter embedded in another obeys
#   dA/dt = -2 pi M_GB sigma
# exactly (independent of interface width). With M_GB = sigma = 1 that is
# -6.2832 area units per unit time. Same material parameters, same interface
# width and same grid spacing as the n=16 polycrystal relaxation, so the
# measured rate applies directly to it.
[Global]
  backend          = numpy
  use_gpu          = false
  precision        = float64
  seed             = 0
  output_dir       = ./results/graingrowth/shrink/solver
[]
[Mesh]
  type             = box_mesh
  dim              = 2
  domain           = '0 0 50 50'
  divisions        = '79 79'
[]
[Variables]
  [./eta]
    family         = fem
    components     = 2
    order          = 2
  []
[]
[Materials]
  [./gb]
    kind           = gb_evolution
    M_GB           = 1.0
    w_GB           = 1.903763
    sigma          = 1.0
    gamma          = 1.5
  []
[]
[Kernels]
  [./deta_dt]
    type           = time_derivative
    variable       = eta
    coefficient    = 1.0
  []
  [./bulk]
    type           = ac_bulk
    variable       = eta
    material       = gb
    spatial_dim    = 2
  []
  [./interface]
    type           = ac_interface
    variable       = eta
    material       = gb
    spatial_dim    = 2
  []
[]
[ICs]
  [./circle]
    variable       = eta
    kind           = circle
    center         = '25.0 25.0'
    radius         = 15.0
    interface_width = 1.903763
    inside         = 1.0
    outside        = 0.0
  []
[]
[Utilities]
  [./count]
    kind           = feature_flood_count
    threshold      = 0.5
  []
[]
[Executioner]
  type             = pf_fem_transient
  scheme           = convex_splitting
  n_steps          = 20
  dt               = 3.62
  linear_solver    = direct
[]
[Output]
  formats          = 'csv json'
  fields           = 'eta'
[]
