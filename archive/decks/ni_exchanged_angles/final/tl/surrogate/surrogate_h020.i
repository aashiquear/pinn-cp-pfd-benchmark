# Ni specimen surrogate with the structure of docs/13 P7 (PINNACLE ADR 0030 addendum): a network of
# the orientation (three Bunge angles, entered as the 12 sorted Schmid factors of the loading axis x)
# and the five gage station thicknesses to the engineering stress-strain curve (learned as the ratio
# to the calibrated CP set's own material-point curve at the sample's orientation), with the failure
# strain, the fracture location and the peak stress on their own heads; a five-seed ensemble flags
# queries outside the training cloud. Trained on the 30 synthetic h020 specimens of
# decks/ni/final/tl/surrogate_design_h020.i (targets from decks/ni/final/tl/surrogate/targets_h020.i);
# the 15 real specimens are never training data: they are the `evaluate` rows (real_specimens.csv).
# No multi-fidelity: h040 is 6-20% above h020 on the curve (docs/12 §5.2), so every sample is h020.
#
# The mesh, variable and kernel exist only so the problem builds; nothing is solved here. The
# thicknesses are data-only parameters (path = none): they enter through each sample's own mesh.
# missing_samples = skip: a synthetic run that has not finished is left out, never re-run here.
# Restart of Part III (2026-11-15): single parameter set, total-Lagrangian FE runs.
# Run: pinnacle run decks/ni/final/tl/surrogate/surrogate_h020.i
[Global]
  backend    = torch
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = ../../../../../results/ni_study/tl/final/surrogate/model_h020
[]
[Problem]
  type          = composed
  dimension     = 3
  finite_strain = true
[]
[Mesh]
  type   = file
  path   = ../../../../../results/ni_study/mesh/h020/L.1/mesh.e
  format = exodus
  scale  = 1.0
[]
[Variables]
  [./u]
    family     = fem
    components = 3
    order      = 2
  []
[]
[Materials]
  [./cp]                         # the single CP set (decks/ni/final/tl/fem_base_h020.i)
    type                  = dislocation_density_cp
    C11                   = 253000
    C12                   = 152000
    C44                   = 124000
    mu                    = 26000
    gamma0                = 1e-05
    rate_exponent_n       = 8.33333333333
    tau_c0                = 0.26
    q_lat                 = 1.4
    alpha0                = 0.4
    burgers               = 2.5e-07
    rho0                  = 2.0257e6
    k1                    = 0.0243375
    annihilation_diameter = 4.6689e-07
    [./crystal]
      structure    = fcc
      euler_angles = '0 0 0'     # set per sample by [Parameters]
    []
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
[Parameters]
  [./phi1]
    path  = materials/cp/crystal/euler_angles/0
    range = '0 360'
  []
  [./phi]                        # Phi (the parser lowercases block names)
    path  = materials/cp/crystal/euler_angles/1
    range = '0 180'
  []
  [./phi2]
    path  = materials/cp/crystal/euler_angles/2
    range = '0 360'
  []
  [./t1]
    path  = none
    range = '0.285 0.46'           # design: mean 0.325-0.400, taper -0.04 to +0.06 over the gage
  []
  [./t2]
    path  = none
    range = '0.285 0.46'
  []
  [./t3]
    path  = none
    range = '0.285 0.46'
  []
  [./t4]
    path  = none
    range = '0.285 0.46'
  []
  [./t5]
    path  = none
    range = '0.285 0.46'
  []
[]
[Executioner]
  type            = parametric_surrogate
  sources         = data
  target          = reduced
  variable        = u
  columns         = 'eng_stress'
  abscissa        = gauge_strain
  samples_file    = ../../../../../results/ni_study/tl/final/surrogate/targets_h020/samples.csv
  data_dir        = ../../../../../results/ni_study/tl/final/surrogate/targets_h020
  missing_samples = skip
  holdout         = 0.2              # 6 of the 30 synthetic specimens
  seed            = 0
  # P7.1 invariant orientation inputs
  orientation_features   = schmid
  orientation_parameters = 'phi1 phi phi2'
  loading_axis           = x
  crystal_structure      = fcc
  # P7.2 physics-embedded output
  baseline               = material_point
  baseline_material      = cp
  baseline_form          = ratio
  baseline_strain_column = gauge_strain   # e of target.npz
  baseline_strain_rate   = 1.0e-3         # gauge strain rate of the runs (bar 40 % in 800 s)
  baseline_lateral       = free
  baseline_update        = implicit
  baseline_increments    = 200
  # P7.3 separate heads on a shared trunk
  scalars                = 'failure_gauge_strain location_fraction peak_eng_stress'
  architecture           = trunk_heads
  trunk_layers           = 2
  head_layers            = 1
  neurons                = 32
  drop_fraction          = 0.03           # the curve head's own failure estimate (3% load drop)
  # ensemble
  ensemble               = 5
  ensemble_flag          = 0.1
  epochs_adam            = 3000
  epochs_lbfgs           = 50
  lr                     = 2.0e-3
  certify                = false
  evaluate               = ../../../../../results/ni_study/tl/final/surrogate/targets_h020/real_specimens.csv
[]
[Output]
  formats = 'csv json'
  fields  = 'displacement'
[]
