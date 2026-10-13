# Neural-route check at h40 against the FE reference (the calibrated run):
# curves and equivalent-plastic-strain fields at four strains, element by
# element on the same mesh. Read by scripts/field_compare.py.
[Compare]
  figures    = ../../figures/report
  prefix     = neural_h40
  strains    = '2 8 15 20'
  strain_max = 0.20
  vmax       = 0.6
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
[]
[Runs]
  [./fe]
    dir   = ../../results/cp_study/cal/gn2/base
    label = 'FE (cold Newton)'
  []
  [./linear]
    dir   = ../../results/cp_study/hybrid/h40_linear
    label = 'PINN+FEM, linear-extrapolation predictor + Newton'
  []
  [./hybrid]
    dir   = ../../results/cp_study/hybrid/h40_physics
    label = 'PINN+FEM, physics-trained predictor + Newton'
  []
  [./pinn]
    dir   = ../../results/cp_study/hybrid/h40_pinn_solver
    label = 'PINN as solver (verified, corrector = auto)'
  []
[]
