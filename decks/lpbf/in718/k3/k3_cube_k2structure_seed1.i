# Part V, case K3: tensile test of a cube of the SIMULATED build (grain structure of case K2).
#
# PROBLEM: the process-to-property step of the chain. A cube of 120 um is cut out of the
#   built material of the four-layer column of case K2 (k2_column_h10.i: x and y from 240
#   to 360 um, z from 220 to 340 um, above the plate surface at 200 um) and pulled along
#   z, the build direction, to 20 % engineering strain at 1e-3 /s. Its grains and their
#   orientations are those the process simulation produced
#   (`scripts/lpbf_k3.py cube`, file grains_k2_column_seed1.csv); the mesh is the
#   column's (10 um, 12 x 12 x 12 elements). The crystal-plasticity set is the calibrated
#   one of k3_cube_n6_seed1.i; nothing is fitted here.
#   Grain realization: column seed 1.
#   Target: the mean engineering curve of the eight NIST specimens
#   (benzing 2025 tensiledata; data/lpbf/in718/k3_tensile_mean_curve.csv).
#
# WHY THIS DECK EXISTS: it holds the CALIBRATED crystal-plasticity set of case K3:
#   tau_c0 = 263.37 MPa, k1 = 0.052690, annihilation_diameter = 1.4891 nm, the result of
#   five Gauss-Newton steps of `scripts/lpbf_k3.py calibrate k3_cube_n6_start.i` (2.66 MPa
#   rms against the measured curve between 1 % and 20 % strain; the starting set of
#   k3_cube_n6_start.i gave 66.6 MPa). The MATERIAL text below describes how the
#   STARTING values were derived.
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
  output_dir = ../../../../results/lpbf/in718/k3/k2_structure_seed1
[]

[Mesh]
  type      = box_mesh
  domain    = '0.0 0.0 0.0 120.0e-6 120.0e-6 120.0e-6'
  divisions = '12 12 12'
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
    tau_c0                = 2.633653e+08
    q_lat                 = 1.4
    alpha0                = 0.3
    burgers               = 2.54e-10
    rho0                  = 1.0e10
    k1                    = 5.268964e-02
    annihilation_diameter = 1.489077e-09
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'       # replaced per point by [OrientationField]
    []
  []
[]

[OrientationField]
  [./source]
    type          = point_file
    path          = grains_k2_column_seed1.csv
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
