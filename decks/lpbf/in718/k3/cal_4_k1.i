# Part V, case K3: tensile test of a polycrystal cube of alloy 718, starting parameters.
#
# PROBLEM: a cube of 120 um with 216 Voronoi grains of 20 um and uniform random
#   orientations (the measured texture of the build is close to random: fractions
#   within 15 degrees of <100>, <110>, <111> of 0.10, 0.21, 0.10 against 0.10, 0.20,
#   0.14) is pulled along z, the build direction, to 20 % engineering strain at 1e-3 /s.
#   Target: the mean engineering curve of the eight NIST specimens
#   (benzing 2025 tensiledata; data/lpbf/in718/k3_tensile_mean_curve.csv).
#
# WHY THIS DECK EXISTS: it is the forward model of the crystal-plasticity calibration of
#   case K3 (scripts/lpbf_k3.py calibrate). This deck holds the STARTING set; the
#   calibration writes one deck per iteration next to it.
#
# MATERIAL: `dislocation_density_cp`, the law of Parts I to III of this benchmark
#   (finite strain, backward-Euler slip update, consistent tangent):
#       gamma_dot = gamma0 |tau / tau_c|^n sign(tau),
#       tau_c = tau_c0 + alpha0 mu b sqrt(sum_b q_ab rho_b),
#       d rho = (k1 sqrt(rho) - 2 y_c rho) |d gamma| / b.
#   STARTING VALUES are derived from the calibrated crystal-plasticity set of the
#   author's toolkit (AM-Orientation-Shift-TOOLKIT, calibrations/in718/tensile_cp.json:
#   a Voce-type law with g0 = 227.2 MPa, saturation 666.3 MPa, h0 = 1,407 MPa,
#   exponent 3.10, q = 1.4, gamma0 = 1e-3 /s, m = 0.05, fitted to the same NIST curve
#   with a Taylor model of 400 grains): the same initial resistance, the same
#   saturation resistance and the same initial hardening rate h0 (1 - g0/gs)^a = 386 MPa:
#       tau_c0 = g0 - alpha0 mu b sqrt(S rho0),   S = 1 + 11 q,
#       k1 = 2 sqrt(S) (386 MPa) / (alpha0 mu),   y_c = k1 / (2 sqrt(rho_sat)).
#   Elastic constants, rate constants and q are the toolkit's and are held. Assumed:
#   alpha0 = 0.3, mu = 81 GPa (mean of the Voigt and Reuss shear moduli of the elastic
#   constants), b = 0.254 nm, rho0 = 1e10 /m2.
#   FITTED by the calibration: tau_c0, k1, annihilation_diameter.
#
# LOADING: displacement of the top face (zmax) along z, linear in time to 24 um (20 %);
#   the bottom face is held in z, the faces xmin and ymin are symmetry planes, the faces
#   xmax, ymax are free: uniaxial stress on average.
#
# TIME: physical; 100 load steps of 2 s (0.2 % strain each).
#
# UNITS: SI (m, Pa, s).

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ../../../../results/lpbf/in718/k3/cal/it4_k1
[]

[Mesh]
  type      = box_mesh
  domain    = '0.0 0.0 0.0 120.0e-6 120.0e-6 120.0e-6'
  divisions = '6 6 6'
[]

[Variables]
  [./u]
    family     = fem
    components = 3
    order      = 2
  []
[]

[Materials]
  [./cp]
    kind                  = dislocation_density_cp
    C11                   = 250.0e9
    C12                   = 160.0e9
    C44                   = 120.0e9
    mu                    = 81.0e9
    gamma0                = 1.0e-3
    rate_exponent_n       = 20.0           # 1 / m, m = 0.05
    tau_c0                = 2.633131e+08
    q_lat                 = 1.4
    alpha0                = 0.3
    burgers               = 2.54e-10
    rho0                  = 1.0e10
    k1                    = 5.543287e-02
    annihilation_diameter = 1.497206e-09
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'       # replaced per point by [OrientationField]
    []
  []
[]

[OrientationField]
  [./source]
    type          = point_file
    path          = grains_seed1.csv
    interpolation = nearest
  []
[]

[Kernels]
  [./sd_cp]
    type          = stress_divergence_cp
    variable      = u
    material      = cp
    dim           = 3
    finite_strain = true
    kinematics    = total_lagrangian
    state_update  = implicit
  []
[]

[BCs]
  [./xmin]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./ymin]
    variable   = u
    kind       = dirichlet
    node_set   = ymin
    components = [1]
    values     = [0.0]
  []
  [./zmin]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
  [./pull]
    variable   = u
    kind       = dirichlet
    node_set   = zmax
    components = [2]
    values     = [24.0e-6]
  []
[]

[Executioner]
  type            = cp_quasi_static
  n_steps         = 100
  t_end           = 200.0
  newton_tol      = 1.0e-9
  newton_rel_tol  = 1.0e-6
  max_iter        = 30
  max_cutbacks    = 50
  linear_method   = direct
  line_search     = backtracking
[]

[Output]
  formats         = 'json'
  save_checkpoint = false
[]
