# Neural routes on the calibrated CP+PFD plate at h40 (docs/09 §5.1 step 6):
# the mechanics sub-solve of every staggered step with (a) the linear
# extrapolation predictor (PINN+FEM, Part I §8) and (b) the grain-aware
# FE-residual-trained network as the verified solver (corrector = auto,
# PINN as solver, Part I §8-9); the base is the FE route. Same deck
# otherwise (Gc 1106, l0 5, 256 steps). Compared on the fracture strain,
# the curve, the fields at matching steps, Newton iterations and time.
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/cp/neural_h40
  concurrent = 2
[]
[Cases]
  [./linear]
    predictor = linear
  []
  [./solver_grain]
    predictor = fem_solver_grain
    corrector = auto
    predictor_device = cuda
    predictor_gpu = true
  []
[]
