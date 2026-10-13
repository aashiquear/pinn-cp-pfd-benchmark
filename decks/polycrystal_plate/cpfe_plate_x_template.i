# Pseudo-3D polycrystal plate, crystal plasticity only (no damage): 100 x 50
# x 1 um, one HEX8 element through the thickness, pulled along +x (the 100 um
# length) from the x = 100 face. Study geometry decided 2026-10-09.
#
# PLACEHOLDERS (@NAME@) are filled by scripts/polycrystal_plate.py from the
# pipeline deck; every value that reaches this file is recorded there.
#
# BOUNDARY CONDITIONS (the MOOSE deck's set, rotated to x loading): u_x = 0 on
# the x = 0 face, u_y = 0 on the y = 0 face, u_z = 0 on the z = 0 face, and
# u_x prescribed on the x = 100 face by the ramp `pull`, which reaches
# @STRAIN@ x 100 um at t_end. The y = 50, z = 1 and loaded faces are
# otherwise free, so the plate contracts laterally and through-thickness
# (generalized plane stress with one element through z).
#
# ORIENTATIONS. One Bunge triple per element, taken from the point file the
# pipeline writes after labelling the phase-field grain structure
# (`[OrientationField] source`, ADR 0026 in PINNACLE). The crystal block's
# own Euler angles are the per-QP field's fallback and are never used when
# the field covers the mesh.
#
# UNITS. micrometer, megapascal, second (the MOOSE deck's unit system).

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = @OUTDIR@
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]

[Mesh]
  type      = box_mesh
  domain    = '0 0 0 100 50 @THICKNESS@'
  divisions = '@NX@ @NY@ 1'
[]

[Variables]
  [./u]
    family     = fem
    components = 3
  []
[]

[Functions]
  [./pull]
    # Top-face displacement: strain x height x t / t_end.
    kind       = parsed
    expression = @STRAIN@*100.0*t/@TEND@
  []
[]

[Materials]
  [./cp]
    type                  = dislocation_density_cp
    C11                   = @C11@
    C12                   = @C12@
    C44                   = @C44@
    mu                    = @MU@
    gamma0                = @GAMMA0@
    rate_exponent_n       = @NEXP@
    tau_c0                = @TAUC0@
    q_lat                 = @QLAT@
    alpha0                = @ALPHA0@
    burgers               = @BURGERS@
    rho0                  = @RHO0@
    k1                    = @K1@
    annihilation_diameter = @YC@
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'
    []
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
  []
[]

[OrientationField]
  [./source]
    type          = point_file
    path          = @POINTFILE@
    interpolation = nearest
  []
[]

[BCs]
  [./x_left]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./y_bottom]
    variable   = u
    kind       = dirichlet
    node_set   = ymin
    components = [1]
    values     = [0.0]
  []
  [./z_back]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
  [./x_pull]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
    function   = pull
  []
[]

[Executioner]
  type             = cp_quasi_static
  n_steps          = @NSTEPS@
  t_end            = @TEND@
  newton_tol       = @NEWTONTOL@
  # Relative test off (0): a warm start from a predictor must be held to the
  # SAME absolute criterion as a cold start; the relative-to-first-residual
  # test penalises a good initial guess (docs/07 §10.0 of the study).
  newton_rel_tol   = @NEWTONRELTOL@
  max_iter         = 50
  max_cutbacks     = @MAXCUTBACKS@
  # direct (sparse LU) by default: conjugate gradients stalls once the
  # finite-strain CP tangent loses positive definiteness (7.8-9.8% strain on
  # this plate); gmres is the iterative alternative for larger meshes.
  linear_method    = @LINEAR@
  sparse_threshold = 0
  # Per-step field snapshots (0 = final state only); N writes
  # fields_step{NNNN}.vtk + fields.pvd every N accepted steps.
  snapshots_every  = @SNAPSHOTS@
  # newton (default) corrects every step; none accepts a [Utilities]
  # predictor's proposal as the solution (PINN as solver, ADR 0028).
  corrector        = @CORRECTOR@
  corrector_tol    = @CORRECTORTOL@
  predictor_reject_ratio = @REJECTRATIO@
[]

@UTILITIES@
[Output]
  formats         = 'vtk csv json'
  save_checkpoint = false
  fields          = 'displacement stress strain equivalent_plastic_strain slip_activity mean_rho'
[]
