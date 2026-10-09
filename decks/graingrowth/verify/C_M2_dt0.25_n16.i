# Grain-growth verification: is the Allen-Cahn mobility connected, and is the
# previous dt = 8/(L mu) saturated under convex splitting?  Same 16-grain
# Voronoi seeding, grid and material as results/graingrowth/n16/deck.i; only
# M_GB, dt and n_steps differ.  Case C_M2_dt0.25_n16: M_GB = 2.0, dt = 0.113259, n_steps = 16.
# Physical horizon T = dt*n_steps; the bulk relaxation time is 1/(L mu) = 0.45304 at M_GB = 1.
[Global]
  backend          = numpy
  use_gpu          = false
  precision        = float64
  seed             = 0
  output_dir       = results/graingrowth_verify/C_M2_dt0.25_n16
[]

[Mesh]
  type             = box_mesh
  dim              = 2
  domain           = '0 0 100 50'
  divisions        = '158 79'
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
    M_GB           = 2.0
    w_GB           = 1.90376
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
    n_grains       = 16
    n_order_params = 8
    # The study's own Voronoi seeds, so the relaxed structure is the same
    # polycrystal the single-crystal orientation sets were built on.
    seeds          = '21.375858 23.045362 75.718029 45.159014 63.156194 29.624736 66.365345 17.009656 24.174361 15.065615 73.147485 23.932946 97.353649 34.135987 64.704549 28.935469 51.892754 37.376832 61.399728 26.761573 64.888254 37.133289 97.137358 40.047967 90.186312 29.428829 18.431662 26.014469 35.953548 35.925936 80.988959 14.328497'
    periodic       = false
    int_width      = 1.90376
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
  scheme           = convex_splitting
  n_steps          = 16
  dt               = 0.113259
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
  formats          = 'csv json'
  fields           = 'eta'
[]
