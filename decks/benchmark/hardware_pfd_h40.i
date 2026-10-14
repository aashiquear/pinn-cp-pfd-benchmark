# Hardware and time benchmark of the CP+PFD routes on the calibrated h40
# perforated plate (docs/09 §5.1 step 7): the cold FE route on numpy and on
# torch CPU (2 and 8 threads) and the torch GPU; the PINN+FEM route (linear
# extrapolation); the verified grain-aware PINN-as-solver route with its
# network on the GPU and on the CPU. Each run is the first 16 staggered steps
# of the full deck (strain per step unchanged, early stop off), alone on the
# host. Read by scripts/hardware_benchmark.py (deck_file = pfd.i).
[Benchmark]
  out_dir   = ../../results/pfd_study/cp/hardware
  figures   = ../../figures/pfd
  prefix    = hardware_pfd_h40
  deck_file = pfd.i
  n_steps   = 16
  repeats   = 1
[]
[Cases]
  [./fe_numpy_t8]
    deck      = ../../results/pfd_study/cp/neural_h40/base/pipeline.i
    backend   = numpy
    threads   = 8
    predictor = none
  []
  [./fe_torch_cpu_t2]
    deck      = ../../results/pfd_study/cp/neural_h40/base/pipeline.i
    backend   = torch
    use_gpu   = false
    threads   = 2
    predictor = none
  []
  [./fe_torch_cpu_t8]
    deck      = ../../results/pfd_study/cp/neural_h40/base/pipeline.i
    backend   = torch
    use_gpu   = false
    threads   = 8
    predictor = none
  []
  [./fe_torch_gpu]
    deck      = ../../results/pfd_study/cp/neural_h40/base/pipeline.i
    backend   = torch
    use_gpu   = true
    threads   = 8
    predictor = none
  []
  [./hybrid_linear_t8]
    deck      = ../../results/pfd_study/cp/neural_h40/linear/pipeline.i
    backend   = torch
    use_gpu   = false
    threads   = 8
  []
  [./solver_grain_gpu]
    deck             = ../../results/pfd_study/cp/neural_h40/solver_grain/pipeline.i
    backend          = torch
    use_gpu          = false
    threads          = 8
    predictor_device = cuda
  []
  [./solver_grain_cpu]
    deck             = ../../results/pfd_study/cp/neural_h40/solver_grain/pipeline.i
    backend          = torch
    use_gpu          = false
    threads          = 8
    predictor_device = cpu
  []
[]
