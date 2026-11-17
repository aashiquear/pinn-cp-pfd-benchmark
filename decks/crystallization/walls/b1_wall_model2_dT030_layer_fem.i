# Part IV, case B1 (FE route): heterogeneous crystal nucleus on a foreign wall, Model II,
# wall value phi_0 = 2/3 (contact angle 61.2 deg), undercooling 30 K: the wall layer WITHOUT a nucleus (reference level)
# (review Fig. 6, central row and dashed line).
#
# PHYSICS. Single-field model (granasy 2019 nucleation, Eq. 12, 17-19; granasy 2007
# heterogeneous), phi = 0 liquid, phi = 1 crystal, two dimensions:
#
#     F = int { eps^2/2 |grad phi|^2 + w g(phi) + p(phi) Delta_f } dA
#
# Model II prescribes the order parameter on the wall, phi = phi_0, which orders the liquid
# next to it: a diffuse wall layer. The contact angle follows from
# cos(Psi) = -1 + 6 phi_0^2 - 4 phi_0^3 (phi_0 = 2/3: Psi = 61.2 deg). The barrier is the
# free energy of the nucleus on the wall minus that of the wall layer alone:
#
#     W = 2 [ F(nucleus deck) - F(layer deck) ]      (2: the mirror half)
#
# The layer exists while w g(phi) + p(phi) Delta_f > 0 for 0 < phi <= phi_0; at
# Delta_f = -w g(2/3)/p(2/3) = -w/64 it becomes unstable and the barrier vanishes: the
# surface spinodal, 92.0 K in the review.
#
# UNITS AND PARAMETERS. Scaled: lengths in delta = sqrt(2 eps^2 / w) = 0.4551 nm, energies
# per unit length in w(T) delta^2; eps^2 = 0.5, w = 1, gamma = 1/12. Ni data of the original
# paper (granasy 2007 heterogeneous, footnote 16; `reported_text`): gamma_SL = 0.364 J/m^2,
# interface thickness 2 nm, L = 2.666e9 J/m^3, T_m = 1728 K, eps^2 and w multiplied by T:
#     D = Delta_f / w(T) = -0.277781 Delta T / (1728 K - Delta T) = -0.0049078 at 30 K.
#
# HOW IT IS SOLVED. Half of the symmetric state on [0, 45] x [0, 25]: phi = 2/3 on the
# wall y = 0, symmetry (natural) at x = 0, liquid on the top side, natural on the right
# side (the wall layer continues there). One backward-Euler step with dt = 1e12 (Newton on the stationary weak form; the layer is a minimum), started from
# the wall layer of the melting point.
#
# EXPECTED RESULT. `free_energy` of postprocessors.csv is F of this deck; the barrier needs
# the partner deck. Classical cap: W_cap = 0.91405. Model II lies well below the cap and
# reaches zero at the surface spinodal (review Fig. 6, dashed line); values in the report.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/walls/b1_wall_model2_dT030_layer_fem
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 45 25'
  divisions = '180 100'       # h = 0.25 delta
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # wall layer of the melting point
  [./cap]
    kind       = parsed
    expression = 1/(1+0.5*exp(y))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = cap
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'
    [./coefficients]
      w  = 1.0
      df = -0.0049078
    []
  []
[]

[Kernels]
  [./dphi_dt]
    type        = time_derivative
    variable    = phi
    coefficient = 1.0
  []
  [./bulk]
    type                 = ac_bulk
    variable             = phi
    material             = f_bulk
    mobility             = 1.0
    gradient_coefficient = 0.5
    spatial_dim          = 2
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.5
    spatial_dim          = 2
  []
[]

[BCs]
  # Model II wall: the order parameter is prescribed
  [./wall]
    variable   = phi
    kind       = dirichlet
    node_set   = ymin
    components = [0]
    values     = [0.6666666666667]
  []
  [./liquid_top]
    variable   = phi
    kind       = dirichlet
    node_set   = ymax
    components = [0]
    values     = [0.0]
  []
[]

[Executioner]
  type            = pf_fem_transient
  scheme          = bdf1
  n_steps         = 1
  dt              = 1.0e12
  max_newton_iter = 60
  newton_tol      = 1.0e-9
  max_cutbacks    = 0
  snapshots_every = 1
[]

[Postprocessors]
  [./free_energy]
    type        = total_free_energy
    bulk_field  = field_phi
    free_energy = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.5
      []
    []
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
