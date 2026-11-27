# Runs on several GPUs with PINNACLE: development, verification, timing and cost of decomposed runs and of the library `gpumpi`

Status 2026-11-26. This report belongs to Part V of the study (laser powder bed fusion
process chain, `docs/reports/lpbf_process_pinnacle_benchmark_report.md`), whose decks are
the test cases here. It covers the hardware side only: how a PINNACLE run uses several
GPUs, what was measured, what it cost, and what is not tested. Every run is a row of
`docs/run_ledger.md` (rows `lpbf-hardware-*`, `lpbf-decomposed-chain-checks`,
`lpbf-column-wide-smoke`, `lpbf-dwell-*`); the raw results are in
`results/lpbf/hardware/` (not tracked). Framework decisions: PINNACLE ADR 0050, 0051,
0053 (device pool) and ADR 0057 with its addenda (decomposed runs, `gpumpi`); library
guide: PINNACLE `docs/gpumpi.md`. PINNACLE revision: branch
`feat/multi-gpu-and-am-process`, commit 512717f.

## Abstract

Two capabilities were wanted from a machine with several GPUs: a shorter run, and a mesh
larger than the memory of one card. PINNACLE's first multi-GPU path, a device pool that
splits the element loop of one Python process over the cards, gives neither: on four
NVIDIA A10G cards it changes the time per step of a melt-pool flow run by a factor of
about 1.05 to 1.15, and its first card holds the whole mesh. Four separate processes on
the same cards, each with a quarter of that mesh, run at 0.303 s per step each against
1.255 s for the whole mesh on one card and do not disturb each other: the limit is the
single Python process, not the hardware. The design that follows is the one of the MPI
codes: the mesh is cut into parts, every part is an ordinary PINNACLE problem in its own
process on its own card, and the processes exchange the nodal sums at the cuts. The
exchange layer is a separate Python library, `gpumpi`, that depends on NumPy only and is
usable by other array codes (NumPy, PyTorch, CuPy, JAX): a launcher that gives each
process one device, a communicator close to that of `mpi4py`, mesh partitioning, the sum
of nodal values across cuts, and Krylov solvers for linear systems that live in parts.
With it the whole process chain runs on a cut mesh with the result of the single run:
heat and melt-pool flow to round-off, grains with nucleation to $2\times10^{-14}$, the
crystal-plasticity stress stage to the tolerance of its linear solve. On four NVIDIA
A100 cards a flow track of 588,000 elements advances 3.56 times faster than on one card,
a build column with flow and 30 grain fields 3.77 times, and the stress stage 3.7 times;
the device memory per card is that of its part (36 kB per element with flow and
grains). The cloud measurements cost 19.66 US dollars of a limit of 30, of which about
8.8 were lost to two runs that failed or were stopped for faults of the test scripts and
of one set-up routine, not of the solver. Not tested: the MPI transport, the CuPy and
JAX array paths, more than four cards, a mesh that exceeds one card, and a coarse-space
correction for the linear solve of large stress runs.

## 1. Goals and starting point

The author's two goals:

1. **Speed.** Several GPUs should shorten a run.
2. **Memory.** The number of finite elements of a run should not be bound by the memory
   of one card, the limit of every single-device code.

The machine the runs are meant for has eight NVIDIA H200 cards of 141 GB and is shared:
a run can count on at most four cards, and which ones are free is known only when it
starts.

**Starting point: the device pool.** `[Global] devices = 'cuda:0 cuda:1 ...'` (or
`num_gpus = N`) splits work of one process over several devices of one node (PINNACLE
ADR 0050, 0051, 0053): the element loop of the explicit executioner `pf_fem_explicit`
(each device holds a contiguous range of elements and returns nodal contributions, which
are summed on the first device), the forward pass of a network over its batch of points,
and the batched evaluation of the crystal-plasticity law over its quadrature points.
`cuda:N` addresses NVIDIA cards and, in a ROCm build of PyTorch, AMD cards (`rocm:N`,
`hip:N` are accepted). The pool was verified on the one GPU of the study host by naming
that GPU twice and with the CPU:

| path | check on the host (one AMD Radeon 8060S) | time, single / pooled |
|---|---|---|
| explicit element loop | field equal to $4\times10^{-15}$ | 13.1 s / 14.7 s and 18.6 s |
| network forward pass (three network classes; 200,000 points, second derivatives) | outputs, input derivatives, parameter gradients and ten optimizer steps equal; loss equal to 16 digits | 23.8 s / 24.1 s and 34.8 s for 20 steps |
| crystal-plasticity points (implicit and explicit state update) | stress, displacement and plastic strain equal to $2\times10^{-13}$ | 16.3 s / 17.9 s and 26.1 s |

On six example decks trained with and without a two-slot pool, three agree to $10^{-11}$
or better; the other three differ by as much as the same decks differ between one and
eight threads on one device (round-off amplified by the optimizer). Not on the pool:
sparse linear solves, the host-side part of the implicit slip update, the implicit block
executioner `pf_fem_transient`, the J2 and `kalidindi` materials, and any evaluation
inside a differentiation graph of the constitutive law. On one physical device a pool
costs 5 to 60 % more than the single device, so a statement about speed needed a node
with several cards.

## 2. The device pool on four cards

One cloud container with four NVIDIA A10G cards (24 GB each) and 8 CPU cores, PyTorch
2.15 with CUDA, float64 (`scripts/modal_multigpu.py::main`). The case is the single
track with melt-pool flow of Part V (section 6.5 there) on finer meshes, run for a fixed
number of steps with `devices` naming the first one, two or four cards.

**Table 1.** Explicit executioner with the device pool on one, two and four GPUs.

| case | 1 GPU | 2 GPUs | 4 GPUs |
|---|---|---|---|
| 5 um, 588,000 elements: wall time of 200 steps (s) | 306 | 358 | 279 |
| the same, 100 steps (s), two series | 169, 169 | 192, 191 | 149, 146 |
| per step, from the difference of 200 and 100 steps (s) | 1.37 | 1.66 | 1.30 |
| largest difference of the final temperature to one GPU (K) | | 0.0 | 0.0 |
| 3.5 um, 1,720,000 elements: per step in the solver, split part (s) | 2.1 | | 1.7 |

**Table 2.** Where a run of 30 steps at 5 um spends its time (profile of the main
thread).

| part | 1 GPU | 4 GPUs |
|---|---|---|
| set-up on the host, once (element geometry, lumped operators) | 35 s | 35 s |
| element loop on the devices, per step | 0.71 s | 0.49 s |
| services and bookkeeping on the host, per step (the remainder of 1.37 s and 1.30 s) | about 0.7 s | about 0.8 s |

- **The pool is correct on separate cards.** The tensors of every slot are on its card,
  and one, two and four cards give the same field.
- **It does not shorten this solver.** Four cards change the wall time by $-9$ to
  $-14$ %, two cards by $+13$ to $+17$ %. The part of a step that is split gains a
  factor 1.4 on four cards at 588,000 elements and 1.2 at 1.7 million: the gain does not
  grow with the size of the mesh.
- **Causes ruled out.** A larger share of CPU for the container (8 cores in place of the
  default share) changes nothing. Assembling the nodal residual on each card and sending
  only nodal arrays, in place of element arrays assembled on the first card, changes
  nothing (kept in PINNACLE, since it moves eight times less data).
- **About half of a step is serial work on the host** (copies of the fields for the
  services, a nodal volume recomputed in every step, the lumped update), which no number
  of cards shortens. The profile attributes 0.35 s per step of it to named functions; the
  rest was not resolved.
- **Mechanics.** The polycrystal tensile cube of Part V at 4,096 elements (32,768
  points, three load steps) takes 141, 143 and 142 s on one, two and four cards: no
  change. The sparse solve and the host-side part of the slip update are not on the pool.
- **Memory.** The first device keeps the element arrays of the whole mesh, every device
  holds full nodal arrays, and the services of a deck hold the whole mesh on the host:
  the pool does not raise the mesh size a machine can hold.

A plain loop of large tensor operations on four cards from four threads of the same
container reaches a factor 2.6 over the same work done one card after another
(`scripts/modal_multigpu.py::check`), so the cards can work at the same time; the
kernels of the solver are sequences of several hundred small operations issued from
Python by one thread per card.

## 3. One process per card: the experiment that decided the design

The same case as separate processes on the same four A10G cards: the whole mesh on one
card, a quarter of the elements (the block cut along the scan direction) on one card
alone, and four such quarters at the same time, one process per card
(`scripts/modal_multigpu.py::split`). The four processes exchange nothing, so the last
row is the floor of a domain-decomposed run, before the cost of the exchange. Two step
counts (40 and 120) separate the set-up from the time per step.

**Table 3.** Separate processes on four A10G cards, 5 um mesh.

| run | elements per process | wall time of 40 / 120 steps (s) | per step (s) |
|---|---|---|---|
| whole mesh, one process, one card | 588,000 | 107.3 / 207.7 | 1.255 |
| a quarter, one process, one card, alone | 147,000 | 27.2 / 51.7 | 0.306 |
| four quarters at the same time, one process per card (slowest) | 147,000 | 27.8 / 52.0 | 0.303 |

- **Four processes on four cards do not disturb each other.** All elements advance one
  step in 0.30 s in place of 1.26 s, a factor 4.1, and the set-up shrinks by the same
  factor.
- **The limit of the pool is the single Python process** (its threads and the serial
  work on the host), not the cards.
- **Both goals need the same change:** each process owns a part of the mesh on its own
  card, and neighbors exchange the values of the nodes they share.

## 4. How compiled codes use several GPUs

The codes that moved to Kokkos do not get several GPUs from Kokkos. Kokkos is a
portability layer: one loop body compiled for CUDA, HIP, SYCL or OpenMP, on one device
per process. Several GPUs come from MPI with domain decomposition: one MPI rank per GPU,
each rank owns a subdomain with the nodes it shares with its neighbors, and the ranks
exchange only those values (through GPU-aware MPI where available). LAMMPS, PETSc and
Trilinos applications run this way, and MOOSE lists the assembly across several GPUs
with GPU-aware MPI among the requirements of its Kokkos path. Kokkos added several GPUs
inside one process only recently (version 4.3 for CUDA, 4.6 for HIP), and it is manual;
one rank per GPU remains the standard. JAX has the same structure: it can address
several devices, but a finite-element code still has to decompose its mesh itself.

For PINNACLE this means that PyTorch already plays the role of Kokkos (the same kernels
run on CUDA, ROCm and the CPU); the missing half was the MPI half.

Sources (web pages and preprints read for this section on 2026-11-26; they are not
entries of `bib/references.bib`, and nothing in the results depends on them):

- Kokkos documentation, "Multi-GPU Support":
  <https://kokkos.org/kokkos-core-wiki/API/core/MultiGPUSupport.html>
- "VLCs: Managing Parallelism with Virtualized Libraries" (states the versions and that
  the multi-process MPI design is the standard of Kokkos applications):
  <https://arxiv.org/pdf/2512.04320>
- LAMMPS forum, one MPI rank per GPU with the KOKKOS package:
  <https://matsci.org/t/cores-and-gpus-with-kokkos-and-cuda/26452>
- MOOSE framework requirements (partial assembly across several GPUs with and without
  GPU-aware MPI): <https://mooseframework.inl.gov/sqa/framework_srs.html>
- A PETSc algebraic multigrid on GPUs, one MPI rank per card:
  <https://arxiv.org/pdf/2606.24748>; petsc-users on the mapping of ranks to cards:
  <https://lists.mcs.anl.gov/pipermail/petsc-users/2025-November/052459.html>
- PyTorch forum, one process per GPU:
  <https://discuss.pytorch.org/t/cuda-context-size-per-process/79976>
- JAX on several GPUs:
  <https://codelabs.developers.google.com/jax-gpu-06-multi-gpu-training>

## 5. Design

### 5.1 The model

Every element belongs to exactly one part. A part holds its elements and all nodes they
touch, numbered locally; a node on a cut exists in every part that touches it. An
element loop on part $p$ gives, at such a node, the contribution of that part's elements
only, $R_i^{(p)}$ and $M_i^{(p)}$. The parts exchange these values and add them,

$$R_i = \sum_p R_i^{(p)}, \qquad M_i = \sum_p M_i^{(p)}, \qquad
u_i^{n+1} = u_i^n - \Delta t\,R_i/M_i ,$$

in the order of the part numbers on every part, so that every copy of a shared node
holds the same number bit for bit and receives the same update. Only the nodes on the
cuts travel: a surface against a volume.

### 5.2 The library `gpumpi`

A second package of the PINNACLE repository (`src/gpumpi/`) that depends on NumPy only
and imports nothing from PINNACLE.

| piece | what it does |
|---|---|
| `run(function, devices=...)` | starts the ranks, gives each a communicator and a device, returns their results; a rank that fails stops the run with its traceback |
| transports | `process`: one Python process per rank on one machine, joined by pipes (nothing to install); `thread`: ranks as threads, for tests; `mpi`: the ranks of an `mpiexec` run through `mpi4py`, for several nodes |
| device isolation | `CUDA_VISIBLE_DEVICES`, `HIP_VISIBLE_DEVICES` and `ROCR_VISIBLE_DEVICES` are set before the code of a rank loads: each process sees one card, as `cuda:0` |
| `Comm` | `send`, `recv`, `exchange`, `bcast`, `gather`, `scatter`, `allreduce`, `barrier`, close to the lower-case interface of `mpi4py`; sends never block (a sender thread on pipes, `isend` on MPI), so neighbors cannot deadlock; reductions add in rank order |
| `partition_slabs`, `partition_rcb` | element parts as slices along one axis or by recursive coordinate bisection, balanced by count or weight, cut between layers of elements |
| `decompose`, `Part` | parts with local numbering, the owner of every node, the shared nodes per neighbor in a common order |
| `HaloExchange` | `add` (the sum of section 5.1), `copy_from_owner`, sums and extremes over the whole mesh, for NumPy, PyTorch, CuPy and JAX arrays on any device; only the shared rows are moved to the host and sent |
| `cg`, `bicgstab`, `gmres` | Krylov solvers for a linear system that lives in parts: the caller gives the action of its part's operator and a right-hand side that is complete at the shared nodes; inner products count every node once |
| `Jacobi`, `InteriorSolve` | preconditioners: the assembled diagonal; an exact solve of the nodes only one part holds (SciPy sparse LU on the host) with the diagonal on the cuts |

### 5.3 How PINNACLE uses it

`[Global] rank_devices = 'cuda:0 cuda:1 ...'` (or `ranks = N`; `partition = slabs | rcb`),
or `pinnacle.parallel.run_decomposed(config, devices=[...])`.

1. Rank 0 builds the mesh of the deck, cuts it and sends every rank its part as a mesh
   that keeps the node, side and element sets of the whole mesh.
2. Every rank builds the unchanged deck on its part as an ordinary problem on its own
   device. Kernels, materials, boundary conditions and the backend do not know about
   the decomposition.
3. **One hook.** Components that sum element contributions into nodal arrays carry an
   attribute `interface` and call `interface.add(array)` after the sum: the explicit
   executioner for the lumped mass, the nodal volume, the residual of every row and
   sub-step (in sorted order of the rows), the loads and the projected gradients; the
   process state for the lumped mass of the deposited material. A measure of the whole
   mesh that a kernel needs is taken over all parts.
4. **Random draws are properties of the whole mesh.** The nucleation service draws the
   site, its critical undercooling and its order parameter for every node of the whole
   mesh in the order of the node numbers; a part takes the rows of its nodes. A site
   fires in the part that owns its node; the parts share the fired sites every step and
   each applies the seeds that reach into it. Random plate grains are placed in the box
   of the whole mesh.
5. **Services gather at the end.** The process state sends its element records and
   series frames to the first rank, which writes the files of a single run; the defect
   and cracking reports are made there. During the run the pool measures are reduced
   over the parts every step.
6. **Stress stage.** The Newton solver completes the residual across the cuts, takes
   its norms over the whole mesh, and solves the linear system by GMRES over the parts
   with `InteriorSolve`; the line search probes the completed residual. A constitutive
   update that gives up at one quadrature point is raised on all parts, so that they
   cut the load step back together. Reaction forces and the mean stress and strain of a
   step are reduced over the parts.
7. The fields of all ranks are merged on rank 0 in the node and element numbering of
   the whole mesh, in the array shapes of a single run. A function can reduce every
   stored frame on the rank that made it (`snapshot_filter`).

## 6. Verification

Tier 1 (regression) and tier 2 (independent numerical agreement: the cut run against
the single run of the same deck; the solvers against a direct solution).

**Table 4.** Unit tests (PINNACLE `tests/unit/gpumpi`, `tests/unit/parallel`).

| subject | what is checked |
|---|---|
| communicator, on threads and on processes | tags, collectives, exchange of messages larger than a pipe buffer in rings of two and three ranks, a failing rank, device assignment and isolation |
| partitions and decomposition | balance, cuts between layers, local numbering, every node owned once, shared sets in a common order |
| interface sum | equal to the assembly of the whole mesh with nodes shared by four parts; all copies bit-equal; NumPy and PyTorch arrays |
| solvers | `cg`, `bicgstab`, `gmres` with and without preconditioner against the direct solution of the assembled system on a clamped elastic block (and a non-symmetric variant), two to four parts, NumPy and PyTorch |
| decomposed decks against the single run | conduction track on 2 and 3 slabs and 4 boxes; track with flow on 4 boxes; heat of two layers; grains with nucleation over two layers; record, series, defect and cracking reports file by file; stress stage on two ranks; processes against threads; frame filter; refusals |

Full PINNACLE suite at the revision of this report: 4,320 passed, 58 skipped.

**Table 5.** The cut run against the single run of the same deck.

| case | host | parts | largest difference |
|---|---|---|---|
| conduction track, track with flow, heat of a layered build | workstation CPU, float64 | 2 to 4 | relative $3\times10^{-16}$ to $10^{-14}$ |
| two layers with grains and nucleation | workstation CPU | 3 | $2\times10^{-14}$; 35 nucleation events in both |
| stress after the track (240 elements, line search, 4 cut-backs, 12 increments) | workstation CPU | 2 | displacement $2\times10^{-10}$, stress $6\times10^{-9}$ (relative); the same Newton iterations and cut-backs |
| polycrystal tensile cube (216 elements, implicit slip update, 19 cut-backs, 55 increments) | workstation CPU | 3 | $2\times10^{-13}$ / $10^{-12}$; the same increments |
| the same cube at 4,096 elements, 2 load steps | workstation CPU | 4 processes | $10^{-11}$ / $2\times10^{-11}$ |
| track with flow, 73,500 elements | AMD host, one Radeon card named twice | 2 | run completes, 3.1 GB of host memory against 1.6 GB |
| build column with flow and 30 grain fields, 129,600 elements, 60 steps | AMD host, one card named twice | 2 | temperature and grain fields identical (as stored, float32) |
| track with flow, 588,000 elements, 60 steps | 4 x A100 | 2 and 4 | temperature identical, velocity $1.5\times10^{-16}$ m/s of 0.65 m/s, pressure $10^{-19}$ Pa (as stored, float32) |
| wide build column with flow and 30 grain fields, 518,400 elements, 60 steps | 4 x A100 | 4 | temperature $8\times10^{-12}$ K, grain fields $2\times10^{-15}$, velocity $2\times10^{-14}$ m/s, pressure $2\times10^{-9}$ Pa (float64) |
| tensile cube, 4,096 elements, 2 load steps | 4 x A100 | 4 | stress 0.04 Pa of 1.1 GPa |

The differences of the stress stage are at the tolerance of the linear solve ($10^{-8}$
on the residual). The agreement is implementation consistency, not physical validity:
both runs are the same model.

## 7. Timing

All cloud timings: one container with four cards and 8 CPU cores, float64. The time per
step is the difference of two runs with different step counts, which removes the start
of the processes and the set-up.

**Table 6.** The track with flow at 5 um (588,000 elements) on four NVIDIA A100 cards
(40 GB), slabs along the scan direction, runs of 20 and 60 steps
(`scripts/modal_multigpu.py::cut`).

| run | elements per process | per step (s) | factor | device memory per process (GB) | whole process, 60 steps (s) | host memory, all processes (GB) |
|---|---|---|---|---|---|---|
| one process, one card | 588,000 | 0.730 | 1 | not recorded | 97 | 7 |
| 2 ranks, 2 cards | 294,000 | 0.380 | 1.92 | 2.9 | 56 | 11 |
| 4 ranks, 4 cards | 147,000 | 0.205 | 3.56 | 1.45 | 36 | 17 |

**Table 7.** The grain field and the stress stage on four A100 cards
(`scripts/modal_multigpu.py::chain`). Grains: the wide build column of Part V (section
8.3 there) at 10 um, with flow and 30 grain fields, 518,400 elements, runs of 20 and 60
steps. Stress: the polycrystal tensile cube, two load steps.

| case | one process, one card | four parts, four cards | factor |
|---|---|---|---|
| column with flow and grains, per step (s) | 1.561 | 0.414 | 3.77 |
| the same, set-up and 20 steps (s) | 310 | 82 | 3.8 |
| tensile cube, 4,096 elements, two load steps (s) | 83.6 | 22.4 | 3.7 |
| tensile cube, 13,824 elements, two load steps (s) | not run | 105 | |

**Table 8.** Other timings.

| case | host | single | parts |
|---|---|---|---|
| track with flow, 73,500 elements, 40 steps | AMD host, one card named twice | 16 s | 14 s (2 ranks) |
| column with flow and grains, 129,600 elements, 60 steps | AMD host, one card named twice | 69 s | 62 s (2 ranks) |
| stress after the track, 240 elements | workstation CPU | 214 s | 158 s (2 ranks, threads) |
| tensile cube, 216 elements, 55 increments | workstation CPU | 134 s | 126 s (3 ranks, threads) |
| tensile cube, 4,096 elements, 2 load steps | workstation CPU | 72 s (4 threads, direct solve) | 32 s (4 processes of 4 threads) |

- **Both goals are met.** Four cards advance the mesh 3.6 to 3.8 times faster than one,
  and each card holds its part only.
- **The exchange is small.** In the flow track a part with two neighbors shares 8,662
  of its 155,916 nodes and sends 130 kB per exchange, nine exchanges per step. Of the
  factor 4.1 without any exchange (Table 3), 3.56 remain.
- **The fixed cost shrinks less.** Building the operators takes 46 s for the whole flow
  track on one card and 13 s per rank on four, but starting the processes and loading
  PyTorch costs about 10 s whatever the number of ranks, so a run of 60 steps gains a
  factor 2.7. A production run of thousands of steps sees the per-step factor.
- **The workstation rows are not scaling measurements**: the single run of the cube
  used 4 threads and the four processes 4 threads each. They show that the iterative
  solve of the stress stage is not slower than the direct one at these sizes.
- **Rate for planning:** 3.2 us per element, step and card on an A100 for the column
  with flow and grains at 129,600 elements per card. It is assumed, not measured, that
  the rate holds for larger parts.

## 8. Memory

| quantity | value |
|---|---|
| device memory with flow and 30 grain fields, float64 | 36 kB per element: 4.67 GB per part of 129,600 elements on A100; 4.70 GB for 129,600 elements on one part and 2.37 GB each on two parts on the AMD host; 9.37 GB per part of 259,200 elements |
| device memory of the flow track without grains | 2.9 GB per part of 294,000 elements, 1.45 GB per part of 147,000 (about 10 kB per element) |
| host memory of all processes, flow track at 588,000 elements | 7 GB (one process), 11 GB (two ranks), 17 GB (four ranks) |
| host memory, wide column at 518,400 elements | 22 GB at the peak of a single run followed by four ranks; 5.8 GB for two ranks on the AMD host |
| host memory of the mesh cut | rank 0 holds the whole mesh on the host while it cuts it |
| frames and series | gathered on rank 0 at the end of a run; with the frame filter of the column script about 3 GB at 518,400 elements, 25 GB at 4.1 million, 50 GB at 8.1 million (estimates) |

The device memory of a part scales with the part: the mesh a machine can hold grows
with the number of cards. A mesh that exceeds one card was not run.

## 9. The linear solve of the stress stage

The stress stage no longer factorizes its tangent as a whole. Each Newton step is a
GMRES solve over the parts (restart 80, relative tolerance $10^{-8}$) preconditioned by
`InteriorSolve`.

**Table 9.** Iterations.

| system | unknowns | parts | iterations |
|---|---|---|---|
| clamped elastic block, conjugate gradients to $10^{-10}$, `InteriorSolve` | 4,131 | 2 / 4 / 8 | 58 / 90 / 103 |
| the same, `Jacobi` | 4,131 | 2 / 4 / 8 | 161 |
| clamped elastic block, `InteriorSolve` | 28,611 | 2 / 4 / 8 | 85 / 131 / 153 |
| the same, `Jacobi` | 28,611 | 2 / 4 / 8 | 323 |
| tensile cube, crystal plasticity, GMRES per Newton step | 14,739 (4,096 elements) | 4 | 97 (workstation), 130 (A100 run) |
| the same | 46,875 (13,824 elements) | 4 | 164 |

- The count grows with the mesh and with the number of parts. It is slow enough for
  meshes of this size.
- **No coarse space.** For stress runs of millions of unknowns a coarse-space (two-level)
  correction will be needed; it is not built.
- The preconditioner factorizes the interior of each part on the host with SciPy; the
  factorization and the triangular solves are host work, as the direct solve was.

## 10. Cost of the cloud measurements

Account limit set by the author: 30 US dollars, for checking the capability and its
usefulness. All amounts are metered cost covered by account credits; nothing was billed.
Prices at the time: A10G 1.10, L4 0.80, A100 (40 GB) 2.10 US dollars per card and hour,
CPU about 0.047 per core and hour.

**Table 10.** Cloud runs and their cost (account summary before and after).

| measurement | node | account summary (USD) | cost (USD) | of which lost |
|---|---|---|---|---|
| device pool (Tables 1, 2; seven runs, among them the repeats with a larger CPU share, with assembly per card and with the profiler that make up the diagnosis of section 2) | 4 x A10G | to 4.28 | 4.28 | none |
| one process per card (Table 3) | 4 x A10G | 4.28 to 5.41 | 1.13 | none |
| decomposed flow track (Table 6) | 4 x A100 | 5.41 to 15.63 | 10.2 | about 6.8: the first run was killed for memory in the comparison step of the test script and restarted by the platform; a second run was stopped at the same step after five of its six runs; the third run is the result |
| grain field and stress stage (Table 7) | 4 x A100 | 15.63 to 19.66 | 4.0 | about 2: a first attempt ran a memory-hungry grain set-up and was stopped |
| **total** | | | **19.66** | **about 8.8** |

- Waiting cost nothing: a container with four A10G cards waited 35 minutes and one with
  four L4 cards 20 minutes without being scheduled; the measurements moved to A100
  cards, which were available and cost twice as much per card.
- A complete measurement of six runs on four A100 cards takes 7 to 19 minutes and costs
  1.5 to 2 US dollars.
- The account summary is the figure to go by; an earlier reading of 4.66 for the first
  measurement is not reproduced by it.

## 11. Faults found at size

Three faults did not show in the unit tests and were found by the first runs at size.
All three are fixed and have a test or a guard.

| fault | effect | remedy |
|---|---|---|
| comparison of fields of two runs with different array shapes: a scalar field of the single run is `(n,)`, the decomposed run returned `(n, 1)`, and their difference broadcast to n x n | a cloud container killed for memory, and the GPU host of the study unresponsive for 45 minutes; the solver was not involved | the decomposed run returns the shapes of the single run; the tests compare shapes before values; `scripts/lpbf_parts_check.py` does the same |
| order of the exchanges: rows of fields that are advanced together were completed in the order of a Python set of names, which differs from process to process | a two-process run of the column failed with a shape error; thread ranks share one hash seed and could not show it | the rows are completed in sorted order; a test runs the grain deck on separate processes |
| memory of the initial grains: the distance of every node to every grain seed formed at once | 38 GB of host memory within 17 s for 518,400 elements and 4,000 grains | evaluated in pieces of about 250 MB (1.5 GB at that size, the same result) |

Practice kept from them: every run on a shared or small machine goes through
`scripts/memory_guard.py` (it stops a run before memory runs out and reports what was
left behind), and cloud functions run with restarts switched off.

## 12. Limits and what is not tested

| item | state |
|---|---|
| `mpi` transport of `gpumpi` (several nodes) | written against `mpi4py`; not tested (no MPI on the machines of this study) |
| CuPy and JAX array paths | written against the documented interfaces; not tested (not installed) |
| more than four cards; a mesh that exceeds one card | not run |
| the target machine (H200 cards) | not measured; planning numbers are from A100 cards |
| coarse space for the linear solve of large stress runs | not built |
| periodic boundaries, `pinned_fraction`, noise terms of the explicit executioner | refused or drawn per process on a cut mesh |
| boundary terms that find the boundary from the mesh topology | a cut would count as a boundary; side sets are inherited and work |
| stress stage: predictor and corrector services, gauge strains, element-ordered orientation files, the phase-field fracture coupling, snapshots | not on a cut mesh (orientations given by point file, as in the benchmark decks, work) |
| implicit block executioner `pf_fem_transient`, network executioners | not on a cut mesh; the device pool serves the network forward pass |
| Exodus and VTK output of a cut run | not built; `decomposed_fields.npz` and the files of the gathering services are written |
| services other than the scan path, the process state, the nucleation, and the defect and cracking reports | report per part; their files get `.rankN` |
| the mesh cut | on one rank, which holds the whole mesh on the host |
| collectives | one message per rank (no tree): meant for the ranks of one to a few nodes |

## 13. Reproduction

```bash
# device pool on four cards (Tables 1, 2), threads on four cards, processes without exchange (Table 3)
.venv-modal/bin/modal run scripts/modal_multigpu.py --gpu A10G --counts 1,2,4 --steps 200 --mesh 5
.venv-modal/bin/modal run scripts/modal_multigpu.py::check
.venv-modal/bin/modal run scripts/modal_multigpu.py::split
# decomposed flow track (Table 6); grain field and stress stage (Table 7)
.venv-modal/bin/modal run scripts/modal_multigpu.py::cut --gpu A100 --steps 20,60
.venv-modal/bin/modal run scripts/modal_multigpu.py::chain
# a single run against its cut run on any machine: times, differences per field, per-rank statistics
python3 scripts/lpbf_parts_check.py <deck.i> "cuda:0 cuda:1 cuda:2 cuda:3" [key=value ...]
python3 scripts/lpbf_parts_check.py <deck.i> 4 global.use_gpu=false          # four CPU ranks
# any thermal, flow or grain deck on several cards; the cards with enough free memory
python3 scripts/gpu_free.py <GB needed per card> 4
python3 scripts/lpbf_run.py <deck.i> global.rank_devices="cuda:0 cuda:1"
# under a memory guard
python3 scripts/memory_guard.py <min available GB> <max seconds> <command ...>
# tests (PINNACLE checkout)
python3 -m pytest tests/unit/gpumpi tests/unit/parallel
```

The cloud scripts send the PINNACLE source and the decks of this repository to the
provider and run them there; nothing else leaves the machine. The workstation checks of
Table 5 are kept in reduced form as PINNACLE `tests/unit/parallel/test_decomposed.py`.

## 14. Conclusions

- A single Python process that drives several GPUs does not shorten the explicit solver
  and does not raise the mesh size; one process per card with its own part of the mesh
  does both.
- With the mesh cut into parts the whole process chain of Part V gives the result of
  the single run and gains 3.6 to 3.8 times on four cards, at the device memory of a
  part per card.
- The exchange layer is a library of its own and does not depend on PINNACLE; its MPI
  transport and its CuPy and JAX paths wait for a machine on which they can be tested.
- The next limit is the linear solve of large stress runs, which needs a coarse space.
