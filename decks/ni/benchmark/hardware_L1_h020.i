# Hardware and time benchmark of the routes on the Ni single crystal L.1, h020, with the carried
# CP+PFD set (chi 0, W_c 75 MPa, G_c 0.1 MPa mm, l0 0.6 mm) and the implicit (backward-Euler) CP
# update (author 2026-10-25, before the 15-specimen stage): the cold FE route on torch CPU at
# 2 / 4 / 8 threads, on the torch GPU and on numpy at 8 threads (reference); the PINN+FEM route
# (linear extrapolation); the PINN-as-solver route (fem_residual network) on the GPU and on the
# CPU. Every run is the first 16 load steps of the full 300-step deck (strain per step
# unchanged), one at a time on the host. Read by scripts/hardware_benchmark.py (deck_file pfd.i).
[Benchmark]
  out_dir   = ../../../results/ni_study/benchmark/hardware_L1_h020
  figures   = ../../../figures/ni/benchmark
  prefix    = hardware_L1_h020
  deck_file = pfd.i
  n_steps   = 16
  repeats   = 1
[]
[Cases]
  [./fe_torch_cpu_t2]
    deck      = fe_L1_h020.i
    backend   = torch
    use_gpu   = false
    threads   = 2
    predictor = none
  []
  [./fe_torch_cpu_t4]
    deck      = fe_L1_h020.i
    backend   = torch
    use_gpu   = false
    threads   = 4
    predictor = none
  []
  [./fe_torch_cpu_t8]
    deck      = fe_L1_h020.i
    backend   = torch
    use_gpu   = false
    threads   = 8
    predictor = none
  []
  [./fe_torch_gpu]
    deck      = fe_L1_h020.i
    backend   = torch
    use_gpu   = true
    threads   = 8
    predictor = none
  []
  [./fe_numpy_t8]
    deck      = fe_L1_h020.i
    backend   = numpy
    threads   = 8
    predictor = none
  []
  [./hybrid_linear_t2]
    deck      = linear_L1_h020.i
    backend   = torch
    use_gpu   = false
    threads   = 2
  []
  [./solver_gpu_t2]
    deck             = solver_L1_h020.i
    backend          = torch
    use_gpu          = false
    threads          = 2
    predictor_device = cuda
  []
  [./solver_cpu_t2]
    deck             = solver_L1_h020.i
    backend          = torch
    use_gpu          = false
    threads          = 2
    predictor_device = cpu
  []
[]
