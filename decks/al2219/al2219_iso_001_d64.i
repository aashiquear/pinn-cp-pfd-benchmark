# Al2219 single crystal, spherical void, coupled CP-PFD, ISOCHORIC loading, [001].
#
# Reference configuration for the convergence study. One deck per refinement
# level; only `[Mesh] path` and the mesh-dependent solver counts differ, so the
# levels are comparable by construction.
#
# PARAMETER PROVENANCE: every value below traces to a row in
# data/extracted/benchmark_dataset.csv, extracted independently from
# (Rezwan 2025 intragranular) rather than from the companion data sheet.
# Provenance tags there distinguish reported_table / reported_text / derived /
# calibrated. The two that are NOT simply reported:
#   * k1 = 0.3432 and annihilation_diameter = 0.00268 um are DERIVED from the
#     paper's (k0, y_c) Kocks-Mecking form; see docs/03_material_point_verification.md.
#   * Gc = 22.75 MPa-um is the paper's own CALIBRATED value, and gap C.16 records
#     that it is far below the range expected for a metal, so it is treated as an
#     effective regularized value and swept in the sensitivity study, not as a
#     transferable material property.
#
# UNITS: micrometre - megapascal - second, the source paper's own system. The
# mesh is written in um and read with scale = 1.0.
#
# GEOMETRY: 50 um cube, 5 um void. The void CENTRE is an assumption, not a
# reported value -- gap C.5, `x_void` is unavailable -- so it is centred at
# (25, 25, 25) and that choice is flagged wherever it affects a result.
#
# LOADING: isochoric, verified in iso_bc_verification_cube.i to hold
# det(F) = 1.000000000000. Axial ramp to 12% engineering strain on zmax with the
# two lateral faces following the exact volume-preserving companion; the three
# opposite faces are rollers. The paper reports fracture at 7.0% strain for this
# case, so 12% carries the run past fracture.
#
# DAMAGE BOUNDARY CONDITION -- READ THIS BEFORE CHANGING IT. `pfd_steady`
# defaults `boundary_node_set` to `all_boundary`, and on an imported mesh that
# set is TOPOLOGICAL: it includes the void surface. Using the default would pin
# d = 0 on the void surface, which is the one place the physics requires damage
# to nucleate, and the run would complete and look converged. The mesh generator
# therefore writes an explicit `outer_boundary` set (bounding box only) and a
# `void_surface` set, and this deck names the former. The void surface is left
# free, which is what makes it a stress riser.

[Global]
  # numpy, not torch: the measured backend envelope (docs/02_backend_envelope.md)
  # found the CP field solve 9.9x faster on NumPy CPU than on this box's GPU,
  # because the constitutive update is a per-quadrature-point loop rather than a
  # large dense kernel. The GPU wins for PINN training and is used there. This is
  # a measurement-driven deviation from the project's GPU-first default, recorded
  # as research question R2.
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ./runs/al2219_iso_001_d64
[]

[Problem]
  type          = composed
  dimension     = 3
  # det(F) = 1 is a finite-strain statement, and 12% axial strain is well
  # outside where 1 + tr(eps) would substitute.
  finite_strain = true
[]

[Mesh]
  type   = file
  path   = meshes/al2219/void_r5_L50_d64.e
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
    # 12% of the 50 um edge, ramped over t_end. NEVER quote a parsed
    # expression: the .i grammar tokenises quoted values and the expression
    # arrives as a token list.
    kind       = parsed
    expression = 6.0*t/10.0
  []
  [./lateral_isochoric]
    # u_lat(t) = L * ( (1 + eps_ax(t))^(-1/2) - 1 ), the exact volume-preserving
    # companion to the axial ramp.
    kind       = parsed
    expression = 50.0*((1.0+0.12*t/10.0)**(-0.5)-1.0)
  []
[]

[Materials]
  [./cp]
    type                  = dislocation_density_cp
    C11                   = 108.0e3   # MPa
    C12                   = 61.3e3    # MPa
    C44                   = 28.5e3    # MPa
    mu                    = 25.7e3    # MPa, Taylor-term shear modulus
    gamma0                = 1.0e-5    # reference slip rate, 1/s
    rate_exponent_n       = 8.333333333333334   # 1/m with m = 0.12
    tau_c0                = 0.112     # MPa, athermal slip resistance
    q_lat                 = 1.4       # latent-hardening ratio
    alpha0                = 0.4       # Taylor prefactor
    burgers               = 0.000286  # um
    rho0                  = 220.0     # 1/um^2, initial SSD density
    k1                    = 0.3432    # storage coefficient (derived)
    annihilation_diameter = 0.00268   # um, = 2*y_c (derived)
    [./crystal]
      structure         = fcc
      # [001] loading axis. The rotation ABOUT that axis is undetermined by a
      # Miller direction (gap C.6); it is fixed explicitly here rather than
      # left implicit. Measured effect: 0.06% on axial stress over 0-60 deg,
      # 56% on the shear response -- negligible for stress-strain, material for
      # damage location and fracture angle.
      loading_direction = [0.0, 0.0, 1.0]
      free_rotation_deg = 0.0
    []
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = 22.75     # MPa-um (calibrated; see gap C.16)
    length_scale = 1.0       # um, the paper's l
    damage_model = at2
    anisotropy   = isotropic
  []
  [./drive]
    # Ductile driving force: H = max(H, psi+_el + zeta * W_p). The plastic-work
    # term is what makes this the papers' model rather than a brittle one --
    # without it slip localization never feeds fracture (gap G1a).
    kind                   = pfd_driving_force
    plastic_work_prefactor = 1.0   # zeta, Eq. (31)
  []
  [./deg]
    kind = pfd_degradation
    k    = 0.001   # positive-definiteness parameter, Table 3
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
  # Rollers: ONE component each. Fully fixing these faces would over-constrain
  # the transverse contraction the volume constraint is derived for.
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
  # Loaded faces. `values` is 1.0 and the magnitude lives in the function, so
  # the constraint algebra stays in one place.
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
  type           = staggered
  n_steps        = 20
  max_iterations = 10
  tolerance      = 1.0e-4
  # Record the damage field at EVERY load step. Mandatory for this study: an
  # end-state scalar cannot show when damage initiated or whether the band is
  # localized at the void, and both are reported results.
  snapshots_every = 1
  [./mechanics]
    variable      = u
    type          = cp_quasi_static
    n_steps       = 1
    newton_tol    = 1.0e-6
    max_iter      = 50
    linear_method = direct
    t_end         = 10.0
  []
  [./damage]
    variable          = d
    type              = pfd_steady
    n_load_steps      = 1
    # NOT all_boundary -- see the header. This is the bounding box only.
    boundary_node_set = outer_boundary
    boundary_damage   = 0.0
  []
[]

[Output]
  formats         = 'json csv'
  save_checkpoint = false
  fields          = 'displacement stress strain damage'
[]
