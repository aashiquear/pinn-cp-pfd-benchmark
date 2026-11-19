# Part IV, case D1 (FE route): which phase forms from a random state at eps = 0.3, mean density -0.12.
#
# PHYSICS. Phase-field crystal model (granasy 2019 nucleation, Eq. 30 and 36, no noise;
# elder 2004 plastic): psi is the reduced particle density,
#
#     F = int { psi/2 [-eps + (1 + lap)^2] psi + psi^4/4 } dA,
#     dpsi/dt = lap(mu),   mu = dF/dpsi = -eps psi + psi^3 + (1 + lap)^2 psi.
#
# The homogeneous state is the liquid; a triangular lattice of density peaks is the crystal
# (in two dimensions), stripes the third phase. All quantities are dimensionless.
#
# HOW IT IS SOLVED. Sixth order in space: mixed form with the chemical potential mu and
# nu = lap(psi) on bilinear elements (PINNACLE ADR 0046), implicit steps with Newton's
# method. The box sides are mirror lines (natural conditions n.grad psi = n.grad nu =
# n.grad mu = 0).
#
# EXPECTED RESULT (one-mode phase diagram, elder 2004 plastic; values in the report):
# stripes around psi0 = 0, the triangular phase at intermediate |psi0|, the liquid (uniform
# field) beyond the liquidus.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/pfc/d1_phase_eps0.3_psi-0.12_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 58.0415796555 43.9822971503'
  divisions = '96 73'
[]

[Variables]
  [./psi]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu]
    family     = fem
    components = 1
    order      = 2
  []
  [./nu]
    family     = fem
    components = 1
    order      = 2
  []
[]

[ICs]
  [./psi0]
    variable = psi
    kind     = uniform_random
    minimum  = -0.17
    maximum  = -0.07
  []
[]

[Kernels]
  [./dpsi_dt]
    type     = coupled_time_derivative
    variable = psi
  []
  [./transport]
    type     = split_ch_wres
    variable = mu
    mobility = 1.0
  []
  [./potential]
    type               = sh_potential
    variable           = psi
    chemical_potential = mu
    laplacian          = nu
    epsilon            = 0.3
  []
  [./laplacian]
    type     = weak_laplacian
    variable = nu
    of       = psi
  []
[]

[Executioner]
  type            = pf_fem_transient
  scheme          = bdf1
  n_steps         = 150
  dt              = 2.0
  max_newton_iter = 25
  newton_tol      = 1.0e-8
  linear_solver   = gmres_lu     # GMRES preconditioned by a reused LU factorization
  snapshots_every = 10
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'psi'
[]
