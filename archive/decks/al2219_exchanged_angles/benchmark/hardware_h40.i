# Hardware and time benchmark of the routes on the calibrated h40 plate (author's
# plan, 2026-10-12): the cold FE arm on numpy (1, 2, 4, 8 threads), on torch CPU
# and on the torch GPU (torch CPU also at 1, 2, 4 threads, added 2026-10-13 to
# size the concurrency of the PFD campaigns); the PINN+FEM route (linear extrapolation); and the
# verified grain-aware PINN-as-solver route with its network on the GPU and on
# the CPU. Every run is the first 16 load steps of the full deck (strain per
# step unchanged), run alone on the host. Read by scripts/hardware_benchmark.py.
[Benchmark]
  out_dir = ../../results/cp_study/hardware
  figures = ../../figures/report
  prefix  = hardware_h40
  n_steps = 16
  repeats = 1
[]
[Cases]
  [./fe_numpy_t1]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = numpy
    threads   = 1
    predictor = none
  []
  [./fe_numpy_t2]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = numpy
    threads   = 2
    predictor = none
  []
  [./fe_numpy_t4]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = numpy
    threads   = 4
    predictor = none
  []
  [./fe_numpy_t8]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = numpy
    threads   = 8
    predictor = none
  []
  [./fe_torch_cpu_t1]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = torch
    use_gpu   = false
    threads   = 1
    predictor = none
  []
  [./fe_torch_cpu_t2]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = torch
    use_gpu   = false
    threads   = 2
    predictor = none
  []
  [./fe_torch_cpu_t4]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = torch
    use_gpu   = false
    threads   = 4
    predictor = none
  []
  [./fe_torch_cpu_t8]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = torch
    use_gpu   = false
    threads   = 8
    predictor = none
  []
  [./fe_torch_gpu]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = torch
    use_gpu   = true
    threads   = 8
    predictor = none
  []
  [./hybrid_linear_t8]
    deck      = ../cp_study/hybrid/h40_linear.i
    backend   = numpy
    threads   = 8
  []
  [./solver_grain_gpu]
    deck             = ../cp_study/hybrid/h40_fem_solver_grain.i
    backend          = numpy
    threads          = 8
    predictor_device = cuda
  []
  [./solver_grain_cpu]
    deck             = ../cp_study/hybrid/h40_fem_solver_grain.i
    backend          = numpy
    threads          = 8
    predictor_device = cpu
  []
[]
