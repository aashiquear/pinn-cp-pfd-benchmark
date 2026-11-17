# Part IV, case A3 (PINN+FEM route): critical nucleus of crystalline Ni in its undercooled melt,
# reduced undercooling Delta T / T_f = 0.275 (T = 1253.02 K).
#
# PHYSICS. Single-field phase-field model of (granasy 2019 nucleation), Eq. 12-15;
# phi = 0 liquid, phi = 1 crystal:
#
#     F[phi] = int { eps^2/2 |grad phi|^2 + w g(phi) + p(phi) Delta_f } dV
#     g = phi^2 (1 - phi)^2 / 4,   p = phi^3 (10 - 15 phi + 6 phi^2)
#
# The critical nucleus is the stationary point (a saddle) of F with phi -> 0 far away:
#
#     eps^2 (phi'' + 2/r phi') = w g'(phi) + p'(phi) Delta_f,   phi'(0) = 0.
#
# Its excess free energy is the nucleation barrier W* = F[phi*] (the liquid has F = 0).
#
# HOW IT IS SOLVED (`pf_pinn_fem_hybrid`). Predictor: a network phi_theta(r) is fitted to a
# ROUGH droplet (the classical radius 2 gamma / |Delta_f| = 0.647 nm, full amplitude, the
# planar interface width) and trained briefly on the squared stationary residual. Corrector:
# one backward-Euler step with a very large time step, i.e. Newton's method on the weak
# stationary equation, started from the network's nodal field. The result is the corrector's
# (FE accuracy). What the network can add is the range of initial guesses from which the
# nucleus is reached. Both sub-solves use the kernels of this deck.
#
# UNITS. Length nm, energy aJ (1e-18 J), energy density aJ/nm^3 = GPa, mobility 1.
#
# PARAMETERS (provenance in docs/14_phase_field_crystallization_plan.md, section 5).
#   Ni, review Table 2 (reported_table): gamma = 0.302 J/m^2, 10-90 % width 1 nm,
#   H_f = 17.29 kJ/mol, T_f = 1728.3 K, rho = 8357 kg/m^3, molar mass 58.6934 g/mol.
#   derived at T_f: eps^2 = 0.412338 aJ/nm, w = 15.92548 GPa, delta = sqrt(2 eps^2/w) = 0.22756 nm
#   derived at T = 0.725 T_f (eps^2 and w proportional to T):
#     eps^2 = 0.298945 aJ/nm,  w = 11.545976 GPa
#     Delta_f = -H_f (T_f - T)/(T_f v_m) = -0.6770002 GPa
#
# EXPECTED RESULT (reference: scripts/cryst_reference.py, shooting solution).
#   W* = 0.33671 aJ = 19.463 kT (kT = 0.0172998 aJ); phi(0) = 0.75187;
#   equimolar radius 0.8539 nm.
#   In results: `energy` of solver_history.csv is W* in aJ.

[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/single_field/a3_ni_nucleus_x0.275_hybrid
[]

[Problem]
  type      = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type              = line_mesh
  domain            = '0.0 12.0'      # nm; ten nucleus radii
  divisions         = 2400            # h = 0.005 nm = delta / 45
  coordinate_system = spherical       # dV = 4 pi r^2 dr; r = 0 is the natural (symmetry) end
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # rough guess: classical droplet, R = 2 gamma(T) / |Delta_f| = 0.6468 nm, amplitude 1
  [./droplet]
    kind       = parsed
    expression = 0.5*(1-tanh((x-0.6468)/(2*0.22756)))   # unquoted, no blanks
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = droplet
  []
[]

[Materials]
  [./f_bulk]
    kind       = parsed_free_energy
    variable   = phi
    expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'
    [./coefficients]
      w  = 11.545976
      df = -0.6770002
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
    gradient_coefficient = 0.298945
    spatial_dim          = 1
  []
  [./interface]
    type                 = ac_interface
    variable             = phi
    mobility             = 1.0
    gradient_coefficient = 0.298945    # eps^2(T)
    spatial_dim          = 1
  []
[]

[BCs]
  # far field: undercooled liquid
  [./liquid]
    variable = phi
    kind     = dirichlet
    node_set = xmax
    components = [0]
    values     = [0.0]
  []
[]

[Executioner]
  type = pf_pinn_fem_hybrid
  [./pinn]
    type                 = pf_pinn_stationary
    hidden_layers        = 3
    neurons              = 32
    collocation_interior = 600
    pretrain_epochs      = 2000
    lbfgs_iters          = 40        # a short training: the corrector finishes the job
    deflation_power      = 1.0       # the liquid phi = 0 is not a minimum of the loss
    collapse_tol         = 0.1
  []
  [./fem]
    type            = pf_fem_transient
    scheme          = bdf1
    n_steps         = 1
    dt              = 1.0e12
    max_newton_iter = 50
    newton_tol      = 1.0e-10
    max_cutbacks    = 0
    snapshots_every = 1
  []
[]

[Postprocessors]
  [./W_star_aJ]
    type          = total_free_energy
    bulk_field    = field_phi
    free_energy   = f_bulk
    [./gradient_fields]
      [./phi]
        field = field_phi
        kappa = 0.298945
      []
    []
  []
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi'
[]
