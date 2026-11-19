# Part IV, case D1 (FE route): the triangular crystal at eps = 0.3, mean density -0.25, 8 elements per lattice constant.
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
# EXPECTED RESULT. The one-mode pattern relaxes to the crystal: amplitude and free-energy
# density (`energy` of solver_history.csv divided by the box area 182.343000) against the
# spectral reference scripts/cryst_pfc_reference.py; the error falls as h^2.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/pfc/d1_triangular_eps0.3_psi-0.25_h8_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 14.5103949139 12.5663706144'
  divisions = '16 14'
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

[Functions]
  [./start]
    kind       = parsed
    expression = -0.25+-0.6*(cos(0.866025403784439*(1*(x-0)+0*(y-0)))*cos(0.5*(-0*(x-0)+1*(y-0)))-0.5*cos((-0*(x-0)+1*(y-0))))
  []
[]

[ICs]
  [./psi0]
    variable = psi
    kind     = from_function
    function = start
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
  n_steps         = 60
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
