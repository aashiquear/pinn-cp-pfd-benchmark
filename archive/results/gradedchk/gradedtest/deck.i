# Pseudo-3D perforated plate, 100 x 50 x 10 um, one element through z.
#
# PLACEHOLDERS, substituted by scripts/run_plate.py:
#   /mnt/extended/projects/pinn-cp-pfd-benchmark/meshes/plate/plate_plain_d1.e [0.0, 0.0, 1.0] 1e+09 1 40 7 1 0.2 /mnt/extended/projects/pinn-cp-pfd-benchmark/results/gradedchk/gradedtest/solver
#   0.26490000000000002 25 0.34320000000000001 0.0026800000000000001 12   (calibrated CP set)
#
# STRESS STATE. Generalized plane stress: both z faces are traction-free, so
# sigma_zz = 0 and the single linear element through the thickness carries a
# constant through-thickness strain. u_z is held on zmin ONLY -- that removes
# the rigid-body translation without constraining the thickness change, and
# because zmax is free the constant sigma_zz must vanish, so the pin costs
# nothing physically.
#
# IN-PLANE CONSTRAINT. u_x = 0 on the whole xmin face and prescribed on xmax.
# Holding a whole face in x still leaves a y-translation, so exactly ONE node
# is pinned in y (node set `pin_y`, the mid-height node on xmin). A y roller
# on a whole face would suppress Poisson contraction there and inflate the
# transverse stress -- the error already diagnosed on the cube study.
#
# DAMAGE BOUNDARY. d = 0 on `grips` (the loaded ends), NOT on the full outer
# boundary and NOT on a set containing the hole. The crack initiates at the
# hole and runs transversely to the free edges, so pinning damage on the
# outer boundary would forbid it reaching the edge and cap the run short of
# final fracture; pinning it on the hole would forbid initiation outright.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/gradedchk/gradedtest/solver
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]

[Mesh]
  type   = file
  path   = /mnt/extended/projects/pinn-cp-pfd-benchmark/meshes/plate/plate_plain_d1.e
  format = exodus
  scale  = 1.0
  verify = true
[]

[Variables]
  [./u]
    family     = fem
    components = 3
  []
  [./d]
    family     = fem
    components = 1
  []
[]

[Functions]
  [./axial]
    kind       = parsed
    # t/200 -- NOT t/10.0. This ramp and the executioner's t_end are
    # the same quantity; hard-coding 10.0 here while t_end is derived from
    # the strain rate makes the ramp t_end/10 times too steep, so the first
    # load step demands the whole deformation at once.
    expression = 0.2*100.0*t/200
  []
[]

[Materials]
  [./cp]
    type                  = dislocation_density_cp
    C11                   = 108.0e3
    C12                   = 61.3e3
    C44                   = 28.5e3
    mu                    = 25.7e3
    gamma0                = 1.0e-5
    # n = 1/m_s. Published m_s = 0.12 gives 8.333; a higher n is less
    # rate sensitive, so it needs LESS overstress to reach a given slip
    # rate (1.26x tau_c at n=20 against 1.74x at n=8.33).
    rate_exponent_n       = 12
    tau_c0                = 0.26490000000000002
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 0.000286
    rho0                  = 25
    k1                    = 0.34320000000000001
    annihilation_diameter = 0.0026800000000000001
    [./crystal]
      structure         = fcc
      loading_direction = [0.0, 0.0, 1.0]
      free_rotation_deg = 0.0
    []
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = 1e+09
    length_scale = 1
    damage_model = at2
    anisotropy   = isotropic
  []
  [./drive]
    kind                   = pfd_driving_force
    plastic_work_prefactor = 1.0
  []
  [./deg]
    kind = pfd_degradation
    k    = 0.001
  []
[]

[Kernels]
  [./sd_cp]
    type          = stress_divergence_cp
    variable      = u
    material      = cp
    dim           = 3
    finite_strain = true
    gnd           = false
    degradation   = deg
  []
  [./pfd]
    type          = pfd_evolution
    variable      = d
    dim           = 3
    fracture      = frac
    driving_force = drive
    degradation   = deg
  []
[]

[BCs]
  # Face rollers with one fixed end: xmin holds u_x = 0, ymin holds u_y = 0,
  # zmin holds u_z = 0, xmax carries the prescribed u_x.
  #
  # Differs from BC_MINIMAL only in the transverse restraint -- a whole FACE
  # of y-rollers instead of a single pin line. It removes the same rigid-body
  # rotation but also forbids the ymin face from narrowing, so the specimen
  # cannot contract freely on that side. Whether that is benign is
  # measurable, not assumed: a thin plate in uniaxial tension must carry
  # sigma_yy ~ 0, so transverse stress near ymin is the test.
  #
  # Differs from BC_CLAMPED in NOT holding u_y or u_z on the LOADED face, so
  # no triaxial grip state is created where the load is applied.
  [./xmin_fixed]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./ymin_roller]
    variable   = u
    kind       = dirichlet
    node_set   = ymin
    components = [1]
    values     = [0.0]
  []
  [./zmin_roller]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
  [./xmax_tension]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
    function   = axial
  []
[]




[Executioner]
  type            = staggered
  n_steps         = 40
  # Graded OUTER load steps. This block owns the load factor that is sampled
  # into results.csv, so grading the mechanics sub-solve does nothing for
  # output resolution. Uniform steps put the first reported point past yield.
  graded_first_fraction = 0.0001
  graded_growth         = 1.1
  max_iterations  = 1
  tolerance       = 1.0e-4
  snapshots_every = 1
  [./mechanics]
    variable      = u
    type          = cp_quasi_static
    n_steps       = 7
    newton_tol    = 1.0e-6
    max_iter      = 50
    # The crack tip drives a strain-rate concentration into the undamaged
    # material beside the band, which exhausts the constitutive sub-step
    # budget. The executioner responds by halving the increment, but the
    # default of 6 halvings bottoms out at 1/64 of a load step and the step
    # is then abandoned. Propagation needs a finer floor, and a cutback only
    # costs time on the steps that actually need it.
    max_cutbacks  = 15
    linear_method = cg
    # Always sparse: see the parameter doc. The 12,000-DOF default would put
    # every density level below 5 on a dense tangent, 2.3-5.2x slower.
    sparse_threshold = 0
    # Derived by the harness as strain/strain_rate. The constitutive law
    # is rate dependent, so this must match the rate the CP parameters
    # were calibrated at.
    t_end         = 200
  []
  [./damage]
    variable          = d
    type              = pfd_steady
    n_load_steps      = 1
    boundary_node_set = grips
    boundary_damage   = 0.0
  []
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
  fields          = 'displacement stress strain damage'
[]
