# Part IV, case D2 (FE route): growth of a triangular crystal from a seed, eps = 0.3, mean density -0.34.
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
# SET-UP. A quarter of the crystal: box 94.318 x 94.248 (13 lattice constants by
# 15 row periods), 10 elements per lattice constant, seed of radius 16 (one-mode
# pattern of amplitude -0.5052, that of the crystal at its coexistence density, under a tanh
# envelope of width 3) in the corner. The density is conserved: the
# crystal is denser than the liquid, so a depletion zone forms ahead of the front and the
# growth stops at the lever-rule fraction.
#
# EXPECTED RESULT. The crystalline fraction against time follows the spectral reference
# (scripts/cryst_pfc_reference.py grow 0.3 -0.34 26 15 16 1500 0.125 -0.505164); values in the report.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/pfc/d2_growth_eps0.3_psi-0.34_r16_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 94.3175669402 94.2477796077'
  divisions = '130 130'
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
    expression = -0.34+0.5*(1-tanh((sqrt((x-0)**2+(y-0)**2)-16)/3))*-0.5051642599*(cos(0.866025403784439*(1*(x-0)+0*(y-0)))*cos(0.5*(-0*(x-0)+1*(y-0)))-0.5*cos((-0*(x-0)+1*(y-0))))
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
  n_steps         = 750
  dt              = 2.0
  max_newton_iter = 25
  newton_tol      = 1.0e-8
  linear_solver   = gmres_lu     # GMRES preconditioned by a reused LU factorization
  snapshots_every = 15
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'psi'
[]
