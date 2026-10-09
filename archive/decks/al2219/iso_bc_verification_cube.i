# Verification of the isochoric (constant-volume) loading idiom — no void.
#
# PURPOSE: gap G4 in docs/00_capability_audit.md was recorded as "no isochoric
#   loading driver". This deck shows that no new PINNACLE component is needed:
#   the source paper's constant-volume loading is expressible by COMPOSING the
#   registered `parsed` function with an ordinary `dirichlet` BC, which is what
#   the study plan requires to be attempted before any framework work.
#
# WHAT THE PAPER SPECIFIES (Rezwan 2025 intragranular, Sec. 2.4 / Sec. 3.1.1,
#   independently extracted): "tension applied along the crystallographic
#   direction on one face; compression on the two reciprocal directions; three
#   faces loaded; load ratios adjusted dynamically so that det(F) = 1 globally;
#   remaining faces fixed."
#
# HOW THAT MAPS ONTO A CLOSED FORM: with the three unloaded faces held on their
#   own normals (removing rigid-body motion) and the three loaded faces given
#   prescribed normal displacements, the macroscopic deformation is a uniform
#   stretch diag(lx, ly, lz) with li = 1 + ui/L. Volume preservation is then
#
#       det(F) = lx * ly * lz = 1.
#
#   Splitting the compensation equally between the two lateral directions,
#   ly = lz = lx^(-1/2), satisfies it EXACTLY at every load level — and, since
#   u_lat/u_ax = (lx^(-1/2) - 1)/(lx - 1) is not constant, the displacement
#   ratio must be recomputed at every increment. That is precisely the
#   paper's "load ratios adjusted dynamically", so the phrase describes a
#   nonlinear closed-form relation rather than an equilibrium iteration.
#   No solve-time constraint, and hence no new component, is required.
#
#   The general split ly = lx^(-p), lz = lx^-(1-p) is also exactly isochoric
#   for any p, because the exponents sum to 1. p = 1/2 is the equal-compression
#   reading of "the two reciprocal directions"; a study of unequal p is a
#   sensitivity, not a different model.
#
# UNITS: micrometre - megapascal - second, matching the source paper.
#   Domain edge L = 50 um; nominal axial engineering strain 12% at t_end,
#   which is the strain at which the paper reports its void-growth ratio.
#
# NOT THE PRODUCTION DECK: no void, coarse 6^3 mesh, elastic-plastic CP only
#   with damage switched off, because the quantity under test here is the
#   BOUNDARY CONDITION, and a void or a damage field would only obscure
#   whether det(F) = 1 holds. The production decks add both.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ./out/iso_bc_verification_cube
[]

[Problem]
  type      = composed
  dimension = 3
  # Finite strain: the volume constraint det(F) = 1 is a finite-strain
  # statement, and the 12% axial strain is well outside the small-strain range
  # where det(F) ~ 1 + tr(eps) would be an acceptable stand-in.
  finite_strain = true
[]

[Mesh]
  # `domain` is min-corner then max-corner, not per-axis pairs.
  type      = box_mesh
  domain    = [0.0, 0.0, 0.0, 50.0, 50.0, 50.0]
  divisions = [6, 6, 6]
[]

[Variables]
  [./u]
    family     = fem
    components = 3
  []
[]

[Functions]
  # Axial: linear ramp to 12% engineering strain over t_end = 10 s, i.e. a
  # prescribed displacement of 6.0 um on a 50 um edge.
  [./axial]
    kind       = parsed
    expression = 6.0*t/10.0
  []
  # Lateral: the exact volume-preserving companion,
  #   u_lat(t) = L * ( (1 + eps_ax(t))^(-1/2) - 1 ),  eps_ax(t) = 0.12*t/10.
  # Written out rather than hidden in a component so the constraint is
  # readable in the deck and checkable by hand.
  [./lateral_isochoric]
    kind       = parsed
    expression = 50.0*((1.0+0.12*t/10.0)**(-0.5)-1.0)
  []
[]

[Materials]
  [./cp]
    kind              = dislocation_density_cp
    C11               = 108.0e3
    C12               = 61.3e3
    C44               = 28.5e3
    mu                = 25.7e3
    gamma0            = 1.0e-5
    rate_exponent_n   = 8.333333333333334
    tau_c0            = 0.112
    q_lat             = 1.4
    alpha0            = 0.4
    burgers           = 0.000286
    rho0              = 220.0
    k1                = 0.3432
    annihilation_diameter = 0.00268
    [./crystal]
      structure         = fcc
      loading_direction = [0.0, 0.0, 1.0]
      free_rotation_deg = 0.0
    []
  []
[]

[Kernels]
  [./sd_cp]
    type     = stress_divergence_cp
    variable = u
    material = cp
    dim      = 3
  []
[]

[BCs]
  # Three faces held on their own normal: removes rigid-body motion and makes
  # the origin corner the reference point of the uniform stretch. Each face is
  # constrained in ONE component only, so the faces are free to slide in plane
  # (a roller). Fully fixing them would over-constrain the transverse
  # contraction and the deformation would no longer be the uniform stretch the
  # volume constraint is derived for.
  [./xmin_roller]
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
  # Loaded faces. `values` is 1.0 and the whole displacement comes from the
  # function, because a BC's value is SCALED by its function
  # (BoundaryCondition.function_scale), so putting the magnitude in the
  # function keeps the constraint algebra in one place.
  [./zmax_tension]
    variable   = u
    kind       = dirichlet
    node_set   = zmax
    components = [2]
    values     = [1.0]
    function   = axial
  []
  [./xmax_lateral]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
    function   = lateral_isochoric
  []
  [./ymax_lateral]
    variable   = u
    kind       = dirichlet
    node_set   = ymax
    components = [1]
    values     = [1.0]
    function   = lateral_isochoric
  []
[]

[Executioner]
  type          = cp_quasi_static
  variable      = u
  n_steps       = 10
  newton_tol    = 1.0e-6
  max_iter      = 50
  linear_method = cg
  t_end         = 10.0
[]

[Output]
  formats = 'json'
  fields  = 'displacement stress strain'
[]
