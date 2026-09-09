# Backend and device envelope — research question R2

Phase 0, Framework track. Answers *"under what conditions does GPU execution
accelerate each solver route, and where are the crossover points?"*

Host: `ssh:rocm-server`, AMD Ryzen AI MAX+ 395 (32 threads) with an integrated
Radeon 8060S (`gfx1151`, 40 CU, unified memory), PyTorch 2.13.0+rocm7.2,
`float64` throughout, load average 0.02 at submission. Harness:
`scripts/bench_backend_envelope.py`; raw records in `results/backend_envelope.csv`.
Protocol: capability check before timing, one configuration per subprocess,
3 s time-budgeted warmup with a minimum of 8 iterations, median of 11
individually synchronized reps with min–max spread reported, all loops bounded,
`TORCH_BLAS_PREFER_HIPBLASLT=0`, `OMP_NUM_THREADS=32`. All 46 configurations
completed; none capped or timed out.

---

## 1. Result

**The two halves of the coupled problem want opposite hardware.**

| Workload | NumPy CPU | Torch CPU | Torch GPU |
|---|---|---|---|
| CP constitutive update, µs per quadrature point | **10.7** | 45.3 (4.2× slower) | 106.0 (9.9× slower) |
| PINN training step, µs per collocation point | not applicable | 4.46 | **1.51** (2.9× faster) |

Both rates are the asymptotic values at the largest size measured (65,536
quadrature points; 200,000 collocation points).

**Crystal plasticity: no crossover exists.** NumPy is flat at 9.8–10.7 µs per
quadrature point from 64 to 65,536 points, with run-to-run spread below 1 %.
Torch CPU and Torch GPU are also flat once past their small-size noise, so the
gap is a *constant multiplicative penalty per quadrature point*, not a fixed
overhead that amortizes with batch size. Extrapolation is therefore
unnecessary: the GPU does not become competitive at any batch size, because the
work is a large number of independent 3×3 tensor operations with no dense
parallel structure to exploit, and every one pays dispatch latency. This
reproduces ADR 0015's observation from the opposite direction — the ADR recorded
a bundled bicrystal deck running in 17 s on NumPy and not finishing in 50
minutes on Torch — and quantifies the per-point cause.

**PINN training: the GPU wins everywhere.** The GPU is faster at every size
tested, from 1.4× at 1,000 collocation points to 2.9× at 200,000, and the
advantage is stable rather than growing, indicating both configurations are in
their asymptotic regime by 20,000 points. The GPU also wins on *reproducibility*,
which matters for the study's repeated-run timing statistics: spread is
0.17–2.4 % on the GPU against 6.9–85 % on the CPU, where the double-backward
graph competes with 32 threads for memory bandwidth. The 80,000-point CPU
measurement (8.65 µs/point, 85 % spread) is the clearest instance and should be
read as contention, not as a genuine scaling feature.

## 2. State marshalling is not the bottleneck

`compute_batch` receives per-quadrature-point state as a Python list of dicts
and stacks it on every call, which looked like a plausible cause of the Torch
penalty. Measured separately, it is not: marshalling is **3.8 %** of NumPy CP
time, **0.4 %** of Torch CPU time, and **0.2 %** of Torch GPU time. Torch is
actually *faster* than NumPy at the stacking itself (0.18 vs 0.41 µs per point).

This is a useful negative result. A struct-of-arrays refactor of the state
layout would not recover the Torch penalty; the cost is in the constitutive
arithmetic dispatch, which `compute_batch` already batches as far as the current
formulation allows. Closing the gap would require changing the formulation, not
the data layout — consistent with the survey's finding that GPU acceleration of
crystal plasticity requires algorithmic restructuring rather than porting
element loops.

## 3. G11 — the framework cannot give each half its preferred backend

`backend` and `device` are declared only under `[Global]` in
`schemas/pinnacle_input.schema.json` (`/properties/global/properties/backend`)
and documented only as `[Global]` keys in `docs/user_guide.md` §6. There is one
backend and one device per run, with no per-sub-executioner override.

The all-FEM route is unaffected: it runs entirely on NumPy, its best backend.
The hybrid routes are affected directly. NumPy has no autodiff, so a hybrid run
*must* select Torch, and its CP half then pays the penalty above. The only
choice left is Torch CPU versus Torch GPU, which trades 60.7 µs per quadrature
point against 2.95 µs per collocation point. Equating the two gives a
break-even:

> Torch GPU beats Torch CPU for a hybrid run only when
> `N_epochs × N_collocation > 20.6 × N_cp_evaluations × N_quadrature_points`.

For a representative 3-D void problem — 20,000 hexahedra at 8 quadrature points
(160,000 points), 5 load steps × 10 Newton iterations, 20,000 collocation
points — the PINN would need roughly **8,200 training epochs** before the GPU
pays for itself. The bundled decks use 2,000. At these settings the hybrid route
is faster on the CPU, and putting it on the GPU because "PINNs belong on GPUs"
would make it slower.

This has a direct bearing on the study's headline cost result. Without
per-solver backend scoping, the hybrid route's measured cost contains a backend
penalty of 15.2 s per constitutive batch evaluation on the GPU, or 5.5 s on the
CPU, at 160,000 quadrature points — neither of which is an intrinsic property of
the hybrid *method*. Reporting it as such would attribute to "PINN overhead"
what is a framework scoping limitation.

**Two ways to handle it, and they are not equivalent.** Either (i) add
per-sub-executioner `backend`/`device` overrides on `feat/pinn_benchmark`, so
the CP half runs on NumPy and the damage half on the Torch GPU and the hybrid is
measured at its true cost; or (ii) keep the single-backend constraint and report
the penalty as a separately itemized line in the cost accounting. Option (i)
changes the study's conclusion about whether hybrid routes are ever
cost-competitive, so the choice is a scope decision rather than an
implementation detail.

## 4. Backend policy adopted for the study

1. **All-FEM reference route: NumPy, CPU, `float64`.** Fastest by a factor of
   4.2 over Torch CPU and 9.9 over the GPU, and the most reproducible.
2. **PINN training, collocation-batch evaluation, autodiff, surrogate training:
   Torch, GPU, `float64`.** 2.9× faster and an order of magnitude more
   reproducible. `float64` is required by the stiff CP kinetics and by the
   double-backward residual; the fp64 matrix-multiply rate on this GPU
   (0.448 TFLOP/s) exceeds the CPU's, and `torch.autograd.grad` with
   `create_graph=True` was verified working on this device.
3. **Hybrid routes: decided per configuration by the break-even above**, with
   the measurement recorded, pending the scope decision in §3.
4. **Never** `autocast` or `GradScaler`: they target fp16/bf16 matrix multiplies
   and interact badly with `create_graph=True`. Never fp16 on this CPU, where it
   is software-emulated and catastrophically slow.
5. Backend and device are recorded per run in the frozen configuration, and any
   hardware effect is reported separately from any algorithmic effect.

The general instruction to prefer Torch on GPU therefore holds for every
neural workload in this study and is overridden, on measured evidence, for the
crystal-plasticity field solve.
