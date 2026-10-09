# Grain-growth PILOT, 32-grain Voronoi seeding on the 100 x 50 um plate (same
# seeds, grid 203 x 102 and gb_evolution material as results/graingrowth/n32/deck.i).
# Fully implicit bdf1 at dt = 0.5/(L mu) = 0.13774 (1/(L mu) = 0.27549, t_int = w_GB^2 = 2.204),
# the step the verification set decks/graingrowth/verify/ showed dt-converged to 0.2%.
# Horizon 40 steps = 5.51 time units = 2.5 interface times.
# Output: final nodal eta field in fields.vtk; energy/gb_area history in csv.
[Global]
  backend          = numpy
  use_gpu          = false
  precision        = float64
  seed             = 0
  output_dir       = results/graingrowth_pilot/n32_bdf1_h40
[]

[Mesh]
  type             = box_mesh
  dim              = 2
  domain           = '0 0 100 50'
  divisions        = '203 102'
[]

[Variables]
  [./eta]
    family         = fem
    components     = 8
    # NOTE: ``order`` here is the Gauss QUADRATURE rule, not the polynomial
    # order of the element. Only order = 2 is implemented; order = 1 raises
    # from deep inside quadrature construction.
    order          = 2
  []
[]

[Materials]
  [./gb]
    kind           = gb_evolution
    M_GB           = 1.0
    w_GB           = 1.48456
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
  [./grains]
    variable       = eta
    kind           = polycrystal_coloring
    n_grains       = 32
    n_order_params = 8
    # The study's own Voronoi seeds, so the relaxed structure is the same
    # polycrystal the single-crystal orientation sets were built on.
    seeds          = '21.375858 23.045362 75.718029 45.159014 63.156194 29.624736 66.365345 17.009656 24.174361 15.065615 73.147485 23.932946 97.353649 34.135987 64.704549 28.935469 51.892754 37.376832 61.399728 26.761573 64.888254 37.133289 97.137358 40.047967 90.186312 29.428829 18.431662 26.014469 35.953548 35.925936 80.988959 14.328497 31.125088 23.619702 28.641826 25.062145 44.007717 22.199897 98.014032 2.291408 33.174592 3.587130 6.374748 32.996269 85.558100 44.833441 14.472919 21.309166 47.722829 22.734432 97.714474 35.063633 79.690805 22.250247 34.506507 47.157881 48.955963 23.341934 49.203973 24.121772 51.921684 1.723780 95.044895 24.627059'
    periodic       = false
    int_width      = 1.48456
    inside         = 1.0
    outside        = 0.0
  []
[]

[Utilities]
  [./count]
    kind           = feature_flood_count
    threshold      = 0.5
  []
  # NOTE: ``grain_tracker`` is deliberately NOT used here. Its ``n_order_params``
  # is the size of an OP POOL that bounds how many grains can coexist, so with
  # grains coloured onto fewer order parameters than there are grains it fails
  # with "no free order parameter left for a new grain". The tracker's persistent
  # identity and remap bookkeeping are for long coarsening runs; this is a short
  # relaxation at constant grain count, and ``feature_flood_count`` already
  # reports the connected-component count that is the stopping criterion.
[]

[Executioner]
  type             = pf_fem_transient
  # convex_splitting, not bdf1. Measured on this problem at 861 nodes, 2 steps:
  # bdf1 at dt = 0.05 costs 41.6 s/step, convex_splitting 18.9 s/step, and
  # convex_splitting is insensitive to dt (18.3 s/step at dt = 0.25) because it
  # linearises the convex part implicitly and gives unconditional energy decay
  # -- which is what the executioner's own documentation recommends when Newton
  # stalls on a stiff well. A first attempt used bdf1 at dt = 1.16, which is
  # 0.32 of the interface relaxation time w_GB^2/(M_GB sigma); the reference
  # MOOSE example's ratio is 0.078, and at 4x that a Newton scheme spends its
  # whole budget on cutbacks -- 55 minutes produced no completed step.
  #
  # Cost is dominated by per-step overhead, not the linear solve: 3.9x the
  # nodes cost 1.1x the time, so the fine grid the diffuse boundary needs is
  # essentially free.
  scheme           = bdf1
  max_newton_iter  = 20
  newton_tol       = 1.0e-8
  max_cutbacks     = 6
  n_steps          = 40
  dt               = 0.137745
  linear_solver    = direct
[]

[Postprocessors]
  [./gb_area]
    type            = grain_boundary_area
    fields          = 'eta'
    grains_per_side = 2
    delta_eta       = 1.0
  []
[]

[Output]
  formats          = 'vtk csv json'
  fields           = 'eta'
[]
