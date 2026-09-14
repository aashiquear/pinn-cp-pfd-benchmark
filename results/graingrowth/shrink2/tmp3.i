# Single-grain shrinkage on the FEM path, reproducing the bundled PINN case
# examples/pinn/phase_field/grain_growth/grain_shrink_law_2d.i EXACTLY in its
# physics and parameters, so the expected answer is a published number rather
# than my own derivation:
#
#     dA/dt = -2 pi L kappa = -2 pi M_GB sigma = -2 pi (1.0)(0.05) = -0.31416
#
# One order parameter shrinking against the eta = 0 background, double-well
# bulk term plus gb_evolution interface term. If the FEM path reproduces
# -0.31416 the kinetics are delivered correctly and the slow polycrystal
# relaxation is physics; if it does not, the kinetics are the defect.
[Global]
  backend          = numpy
  use_gpu          = false
  precision        = float64
  seed             = 0
  output_dir       = ./results/graingrowth/shrink2/tmp3
[]
[Mesh]
  type             = box_mesh
  dim              = 2
  domain           = '0 0 1 1'
  divisions        = '96 96'
[]
[Variables]
  [./eta]
    family         = fem
    components     = 1
    order          = 2
  []
[]
[Materials]
  [./dw]
    kind           = double_well
    barrier_height = 0.1
    min_a          = 0.0
    min_b          = 1.0
  []
  [./gb]
    kind           = gb_evolution
    M_GB           = 1.0
    w_GB           = 0.05
    sigma          = 0.05
    gamma          = 1.5
  []
[]
[Kernels]
  [./td]
    type           = time_derivative
    variable       = eta
    coefficient    = 1.0
  []
  [./bulk]
    type           = ac_bulk
    variable       = eta
    material       = dw
    mobility       = 1.0
    gradient_coefficient = 1.0
    spatial_dim    = 2
  []
  [./iface]
    type           = ac_interface
    variable       = eta
    material       = gb
    mobility       = 1.0
    gradient_coefficient = 1.0
    spatial_dim    = 2
  []
[]
[Functions]
  [./circle_ic]
    kind       = parsed
    expression = 0.5 * (1 + tanh((0.3 - sqrt((x-0.5)**2 + (y-0.5)**2)) / 0.0474))
  []
[]
[ICs]
  [./seed]
    variable = eta
    kind     = from_function
    function = circle_ic
  []
[]
[Executioner]
  type             = pf_fem_transient
  scheme           = convex_splitting
  n_steps          = 120
  dt               = 0.005
  linear_solver    = direct
[]
[Output]
  formats          = 'csv json'
  fields           = 'eta'
[]
