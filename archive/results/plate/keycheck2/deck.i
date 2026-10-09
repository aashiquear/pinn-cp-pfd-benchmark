# Pseudo-3D perforated plate, 100 x 50 x 10 um, one element through z.
#
# PLACEHOLDERS, substituted by scripts/run_plate.py:
#   /mnt/extended/projects/pinn-cp-pfd-benchmark/meshes/plate/plate_hole_d1.e [1.0, 0.0, 0.0] 30 1.5 3 2 2 0.006 /mnt/extended/projects/pinn-cp-pfd-benchmark/results/plate/keycheck2/solver
#   60.31771394231891 2.393387495846754 0.08114503403820844 0.0020581167904802733   (calibrated CP set)
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
  output_dir = /mnt/extended/projects/pinn-cp-pfd-benchmark/results/plate/keycheck2/solver
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]

[Mesh]
  type   = file
  path   = /mnt/extended/projects/pinn-cp-pfd-benchmark/meshes/plate/plate_hole_d1.e
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
    expression = 0.006*100.0*t/10.0
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
    rate_exponent_n       = 8.333333333333334
    tau_c0                = 60.31771394231891
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 0.000286
    rho0                  = 2.393387495846754
    k1                    = 0.08114503403820844
    annihilation_diameter = 0.0020581167904802733
    [./crystal]
      structure         = fcc
      loading_direction = [1.0, 0.0, 0.0]
      free_rotation_deg = 0.0
    []
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = 30
    length_scale = 1.5
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
  [./xmin_fixed]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./pin_y]
    variable   = u
    kind       = dirichlet
    node_set   = pin_y
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
  n_steps         = 3
  max_iterations  = 2
  tolerance       = 1.0e-4
  snapshots_every = 1
  [./mechanics]
    variable      = u
    type          = cp_quasi_static
    n_steps       = 2
    newton_tol    = 1.0e-6
    max_iter      = 50
    # The crack tip drives a strain-rate concentration into the undamaged
    # material beside the band, which exhausts the constitutive sub-step
    # budget. The executioner responds by halving the increment, but the
    # default of 6 halvings bottoms out at 1/64 of a load step and the step
    # is then abandoned. Propagation needs a finer floor, and a cutback only
    # costs time on the steps that actually need it.
    max_cutbacks  = 6
    linear_method = cg
    # Always sparse: see the parameter doc. The 12,000-DOF default would put
    # every density level below 5 on a dense tangent, 2.3-5.2x slower.
    sparse_threshold = 0
    t_end         = 10.0
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
