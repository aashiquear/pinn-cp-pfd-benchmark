# Part IV, case D3 (FE route): five seeds of different lattice orientation grow and impinge, eps = 0.3, mean density -0.32.
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
# SET-UP. Box 130.594 x 131.947 (18 lattice constants), 10 elements per lattice constant;
# seeds of radius 12 at (30, 28) rotated 0 deg, (98, 34) rotated 12 deg, (62, 70) rotated 27 deg, (24, 104) rotated 18 deg, (104, 100) rotated -9 deg.
# The box sides are mirror lines: a grain that reaches a side meets its mirror image there.
#
# EXPECTED RESULT. The grains grow until they meet; because their lattices do not match,
# the boundaries between them are rows of dislocations (particles with five or seven
# neighbours) whose spacing falls with growing misorientation. Counts in the report.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/pfc/d3_polycrystal_eps0.3_psi-0.32_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 130.593554225 131.946891451'
  divisions = '180 182'
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
    expression = -0.32+0.5*(1-tanh((sqrt((x-30)**2+(y-28)**2)-12)/3))*-0.5051642599*(cos(0.866025403784439*(1*(x-30)+0*(y-28)))*cos(0.5*(-0*(x-30)+1*(y-28)))-0.5*cos((-0*(x-30)+1*(y-28))))+0.5*(1-tanh((sqrt((x-98)**2+(y-34)**2)-12)/3))*-0.5051642599*(cos(0.866025403784439*(0.978147600734*(x-98)+0.207911690818*(y-34)))*cos(0.5*(-0.207911690818*(x-98)+0.978147600734*(y-34)))-0.5*cos((-0.207911690818*(x-98)+0.978147600734*(y-34))))+0.5*(1-tanh((sqrt((x-62)**2+(y-70)**2)-12)/3))*-0.5051642599*(cos(0.866025403784439*(0.891006524188*(x-62)+0.45399049974*(y-70)))*cos(0.5*(-0.45399049974*(x-62)+0.891006524188*(y-70)))-0.5*cos((-0.45399049974*(x-62)+0.891006524188*(y-70))))+0.5*(1-tanh((sqrt((x-24)**2+(y-104)**2)-12)/3))*-0.5051642599*(cos(0.866025403784439*(0.951056516295*(x-24)+0.309016994375*(y-104)))*cos(0.5*(-0.309016994375*(x-24)+0.951056516295*(y-104)))-0.5*cos((-0.309016994375*(x-24)+0.951056516295*(y-104))))+0.5*(1-tanh((sqrt((x-104)**2+(y-100)**2)-12)/3))*-0.5051642599*(cos(0.866025403784439*(0.987688340595*(x-104)+-0.15643446504*(y-100)))*cos(0.5*(0.15643446504*(x-104)+0.987688340595*(y-100)))-0.5*cos((0.15643446504*(x-104)+0.987688340595*(y-100))))
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
  n_steps         = 1000
  dt              = 2.0
  max_newton_iter = 25
  newton_tol      = 1.0e-8
  linear_solver   = gmres_lu     # GMRES preconditioned by a reused LU factorization
  snapshots_every = 20
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'psi'
[]
