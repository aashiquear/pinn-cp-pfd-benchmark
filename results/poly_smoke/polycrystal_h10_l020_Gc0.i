# Al2219 30-grain polycrystal calibration specimen (Rezwan 2025, Sec. 2.4).
#
# This is the specimen the Fig. 3 experiment measured: 50 x 50 x 100 um, 30
# grains, uniaxial tension. The test was stopped at 17% strain by catastrophic
# fracture, so 17% is the experimental fracture-onset target for Gc.
#
# PLACEHOLDERS (@NAME@) are filled by scripts/run_polycrystal.py. The
# orientation field is injected programmatically after parsing: it is one Euler
# triple per quadrature point, which is ~128k numbers at h = 2.5 um and does
# not belong in a text deck.
#
# BOUNDARY CONDITIONS. Uniaxial STRESS, not uniaxial strain: three roller faces
# remove rigid-body motion, the loaded face carries the prescribed axial
# displacement, and the two lateral faces are FREE so the specimen contracts as
# a tensile bar does. This is the difference from the void deck, which imposes
# the paper's isochoric condition on the lateral faces.
#
# DAMAGE BOUNDARY. Pinned at the grip end only. Pinning d = 0 over the whole
# outer boundary -- the void deck's choice, correct there because fracture must
# initiate at the void -- would forbid surface initiation here, which is where
# a tensile specimen usually fails.
#
# CRYSTAL BLOCK carries no loading_direction: every quadrature point takes its
# orientation from the injected field, and a global orientation would silently
# override it.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/poly_smoke/run_h10_l020_Gc0
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]

[Mesh]
  type      = box_mesh
  domain    = [0.0, 0.0, 0.0, 50.0, 50.0, 100.0]
  divisions = [5, 5, 10]
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
    expression = 10*t/10.0
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
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 0.000286
    # Calibrated against the nine Vaughan (2024) Fig. 3 markers, NOT the
    # published values: this model forms the Taylor term from the
    # latent-hardened forest sum, so the published set gives tau_c(0) =
    # 176.71 MPa against the paper's reported 43.72.
    tau_c0                = 60.317714
    rho0                  = 2.393387
    k1                    = 0.081145
    annihilation_diameter = 0.002058
    [./crystal]
      structure = fcc
    []
  []
  [./frac]
    kind         = pfd_fracture
    Gc           = 1e+09
    length_scale = 20
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
  [./zmax_tension]
    variable   = u
    kind       = dirichlet
    node_set   = zmax
    components = [2]
    values     = [1.0]
    function   = axial
  []
[]

[Executioner]
  type            = staggered
  n_steps         = 4
  max_iterations  = 1
  tolerance       = 1.0e-4
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
    boundary_node_set = zmax
    boundary_damage   = 0.0
  []
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
  fields          = 'displacement stress strain damage'
[]
