# Part IV, case B1 (PINN+FEM route): heterogeneous crystal nucleus on a foreign wall, Model II,
# wall value phi_0 = 2/3 (contact angle 61.2 deg), undercooling 40 K: the nucleus on the wall
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
#     D = Delta_f / w(T) = -0.277781 Delta T / (1728 K - Delta T) = -0.0065825 at 40 K.
#
# HOW IT IS SOLVED (PINN+FEM route, `pf_pinn_fem_hybrid`). Predictor: the network of the
# PINN route, trained for 40 L-BFGS iterations only. Corrector: one backward-Euler step with a
# very large time step, i.e. PLAIN Newton on the weak stationary equation (no volume
# constraint, unlike the FE deck), started from the network's nodal field. The result is the
# corrector's. The run reports the distance between prediction and corrected field.
#
# EXPECTED RESULT. The state of the FE deck (b1_wall_model2_dT040_nucleus_fem.i, found there by the constrained
# saddle search) to the Newton tolerance. The control deck `b1_wall_model2_dT040_nucleus_newton.i` shows what plain Newton
# does from the circular cap without the network.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 2
  output_dir = results/crystallization/walls/b1_wall_model2_dT040_nucleus_hybrid_seed2
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 40 25'
  divisions = '160 100'       # h = 0.25 delta
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # wall layer of the melting point, with the circular cap of the classical radius on it
  [./cap]
    kind       = parsed
    expression = 1-(1-0.5*(1-tanh((sqrt(x**2+(y+6.0989)**2)-12.6598)/2)))*(1-1/(1+0.5*exp(y)))
  []
  # the wall layer of the melting point: the stationary state without a nucleus (deflated)
  [./layer]
    kind       = parsed
    expression = 1/(1+0.5*exp(y))
  []
  # window of the constrained amount (used by the FE deck only): 1 around the nucleus, 0 beyond 1.5 R along the wall
  [./window]
    kind       = parsed
    expression = 0.5*(1-tanh((x-18.9898)/2))
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
      df = -0.0065825
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
  type = pf_pinn_fem_hybrid
  [./pinn]
    type               = pf_pinn_stationary
    residual_form      = fem         # loss: discrete L2 norm of the assembled weak residual
    hidden_layers      = 4
    neurons            = 48
    pretrain_epochs    = 3000        # regression on the circular cap (the initial condition)
    lbfgs_iters        = 40
    deflation_power    = 1.0         # the state without a nucleus is not a minimum of the loss
    deflation_function = layer
    collapse_tol       = 0.1
    seed               = 2
  []
  [./fem]
    type            = pf_fem_transient
    scheme          = bdf1
    n_steps         = 1
    dt              = 1.0e12        # one backward-Euler step: plain Newton on the stationary equation
    max_newton_iter = 60
    newton_tol      = 1.0e-9
    max_cutbacks    = 0
    snapshots_every = 1
  []
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
