# PINN and PINN+FEM routes for the CP plate: what PINNACLE has, what is missing, and the plan

Assessment of 2026-10-11 (PINNACLE `feat/pinn_benchmark` at 6517a31), by
reading the code, decks, goldens and ADRs; nothing was run. Purpose: decide
what must be developed before the CP-only plate of `docs/reports/` can be
benchmarked on a PINN and on a PINN+FEM hybrid against the FE results of
the report (h40, calibrated set, stress-strain curve AND field values).

## 1. What exists

| route | executioner | state of the CP support |
|---|---|---|
| pure PINN, CP field | `pinn_transient` on `stress_divergence_cp` | the only reachable residual is **displacement-only at the virgin state** (`CPDisplacementOnlyResidual`: F = I + grad u, material evaluated with Fp = I, i.e. anisotropic hyperelasticity); the per-element `[OrientationField]` is **ignored** (material's global Euler angles only); divergence of the Cauchy stress taken in reference coordinates (approximate at finite strain); one steady solve, no load path, no CP state; stress/strain in the result are zeros for the CP kernel; one `material.compute` per collocation point inside the graph |
| pure PINN, full CP state | `_pinn_strong_residual_full_state` + `TransientCPResidual` | unreachable from a deck (`requires_time = False`, raises on finite strain), returns a scalar loss, no orientations |
| PINN+FEM hybrid | `pinn_fem_hybrid` (`warm_start`, `predictor_corrector`, `preconditioner`, `increment_predictor`) | the FEM corrector **never updates the CP state** (virgin-state response at every step), builds the assembler with `finite_strain = False` regardless of the kernel, does not ramp Dirichlet values over the load steps, returns no CP fields and does not keep the PINN prediction for a vector variable; the network has no load input, so the "increment predictor" delta is identically zero; the 1-D CP golden is consistent with an elastic solve (sigma = C11 eps) |
| hybrid with damage as PINN | `pinn_fem_hybrid_coupled_physics mode = cp_pfd` | CP stays on the mesh in `cp_quasi_static`; only damage is neural. Not a CP-only neural route |
| 0-D neural constitutive | `cp_constitutive` | trajectory network N(t) -> state along one homogeneous path; no spatial field; cannot be used at a quadrature point (ADR 0016) |
| backend | `[Global]` only | one backend per problem; PINN training needs torch, so a hybrid forces the CP FEM arm onto torch, where the CP field solve is 4.5x (CPU) to 17x (GPU) slower than numpy (`06` §6); the device-resident batched CP path is scoped, not started |
| PINN field export | `ComposedModel.payload_for` | exports nodal displacement only, evaluating the raw network without the training-time scaling, Fourier features or hard-BC transform; no stress, strain or state |

ADR 0017 states why there is no pure PINN for CP: a coordinates-in network
has no memory, and the CP state is a functional of the load path. Its three
prerequisites for a future route are a differentiable spectral split (damage
only), a mesh-based state update with the interface (F, state, dt) -> state,
and a coupling loop.

**Conclusion.** Neither the pure PINN nor the hybrid can produce a CP plate
result today: the hybrid would return an elastic solution and the PINN an
elastic single-crystal one. Both need development in PINNACLE before any
benchmark run, and the comparison tooling (field metrics against the FE
series) does not exist either.

## 2. What the benchmark needs

The author's criteria: same plate (h40, 30 grains, calibrated set), reach
20% strain, compare the stress-strain curve AND the fields (displacement,
stress, equivalent plastic strain) against the FE run, train on the torch
GPU, FEM arm on whatever is fastest (numpy), and measure cost (wall, device,
epochs, Newton iterations) against accuracy.

Three framework capabilities are implied:

**C1. A load-path hybrid for CP (PINN predictor, FEM corrector, mesh-held
state).** Per load step: the network proposes the displacement increment,
the FEM Newton solve corrects it to the solver tolerance with the finite-
strain CP tangent, the CP state is advanced by the classical integrator, the
step's history row, snapshot and CP fields are recorded exactly as
`cp_quasi_static` does. Measurable quantity: Newton iterations and wall time
saved per step versus the cold FE solve, at identical accuracy (the
corrector converges to the same tolerance). The honest bound: a cold step at
h40 takes 6.4 Newton iterations; a good predictor removes at most 3-4 of
them, so the gain is bounded by about a factor of two on the FE part minus
the training cost. The predictor needs a load input (x, lambda) and must be
trained on something that knows the plastic state: the cheapest consistent
choice is data (the corrected nodal fields of the previous steps, i.e. a
learned extrapolator) plus the strong-form residual evaluated with the
mesh-held state frozen at the quadrature points.

**C2. A per-step "PINN as solver" for the CP field (the pure PINN column).**
The displacement network replaces the FE Newton solve at each load step:
it is trained on the strong-form residual of the finite-strain balance of
momentum with the CP state held at the quadrature points (the quadrature
points are the collocation points, so the state needs no interpolation),
with the per-point orientations and the first Piola-Kirchhoff divergence;
after training, the state is advanced by the classical integrator at the
network's displacement gradient, and the next step continues from the
network. This is the route that answers "can the network produce the FE
fields cheaply": accuracy is not guaranteed by a corrector, so it is
measured against the FE fields; cost is epochs on the GPU per step. It is
the route ADR 0017 calls for in its requirement 2 and the one the capability
audit lists as R5 "to be developed". It is a research item, not a plumbing
item: whether a network can be trained per step fast enough, and whether the
error stays bounded over 256 steps, is exactly what the benchmark would
establish.

**C3. GPU-resident batched CP evaluation.** Both C1's residual and C2's
training evaluate the CP law at every quadrature point inside the autograd
graph thousands of times per step. The measured per-point torch path (one
Python call per point, 15 host-device crossings) is 17x slower than numpy on
this APU; the training would be dominated by it. The device-resident path
(stacked per-QP state, `compute_batch` with per-point orientations
everywhere, no `to_numpy` inside the loop), scoped in
`docs/reports/cp_field_backend_measurement_2026-10-09.md` in PINNACLE, is a
prerequisite for both routes, and it also makes the FE arm usable on torch
if per-arm backends are not wanted.

Plus two smaller items: **per-arm backends** (G11: PINN arm torch GPU, FEM
arm numpy) or, if C3 lands, running the whole hybrid on torch; and the
**comparison tooling** in this repository (field error metrics against the
FE VTK series: L2 and max of displacement, Cauchy stress components,
equivalent plastic strain on the same mesh; curve RMSE; cost ledger).

## 3. Parameters

The constitutive parameters of the PINN routes are the FE calibrated set
(report Table 2) with the fixed set of Table 1; no re-calibration is
involved, the law is the same. What the neural routes add is the training
parameter set: network depth and width, activation, Fourier features,
collocation counts, epochs and learning rates (Adam then L-BFGS), loss
weights (physics, Dirichlet, data), BC mode, resampling, seed. These are the
"extra parameters required for PINN" of the later sensitivity stage; their
baseline values come from a convergence check on the first step (loss and
field error against the FE step), not from the literature.

## 4. Proposed plan, in order

1. **C3 first** (framework, ~1-2 days): device-resident batched CP
   evaluation with per-point orientation and state on torch; verified
   against numpy to 1e-12 on the h40 plate; measured on the host GPU. Gate:
   the 4-step h40 plate on torch-GPU within 2x of numpy (today 17x).
2. **C1, the load-path hybrid** (framework, ~2 days): a `cp_pinn_fem_hybrid`
   executioner built on `cp_quasi_static` (state, snapshots, history,
   cutbacks unchanged) with a per-step predictor; modes `warm_start` (cold
   step 1, network fitted online) and `physics` (residual with frozen
   state); per-arm backend; records per step the PINN prediction, its
   discrepancy from the corrected field, Newton iterations, training epochs
   and wall time. Tests on the 1-D bar (must reproduce `cp_quasi_static` to
   solver tolerance) and a 2-D plate; ADR.
3. **First check** on h40 with the calibrated set (study, 1 run): curve and
   fields versus the FE run; report with figures. Decision point.
4. **C2, PINN as solver** (framework, ~3 days plus tuning): per-step
   residual training with mesh-held state, per-point orientations and the
   first Piola-Kirchhoff form; same outputs as above. First check on h40
   at a few load steps, then to 20% if the error stays bounded. Decision
   point.
5. **Benchmarks** (study): PINN-parameter sensitivity (network size,
   collocation, epochs, weights) on the chosen route(s); hardware and time
   ledger (numpy CPU FE, torch GPU PINN, hybrid) against accuracy (curve
   RMSE, field L2) for the report's §8-9.

Alternatives the author may prefer: skip C2 and benchmark only the hybrid
(then the "PINN" column is a predictor, not a solver, and the claim is
cost at FE accuracy); or do C2 on a reduced problem first (single crystal,
h20, 2% strain) to settle feasibility before the plate.

## 5. Status log

**2026-10-11, C3 done.** PINNACLE ADR 0027 (`batched_constitutive`,
`compute_stacked`, stacked quadrature-point state): the 4-step h40 plate
runs in 6.0 s on numpy, 4.6 s on torch-CPU and 4.55 s on torch-GPU (was 32
/ 145 / 549 s), identical stresses. Gate met. Side effect: the batched
sub-stepper had been evolving the state under the material-global elastic
stiffness on orientation-field decks; fixed, and the pre-fix FE results of
the report are superseded (docs/07 §9; recalibration awaits the author).

**2026-10-11, C1 built.** PINNACLE ADR 0028: `cp_quasi_static` calls
`predict_load_step` / `observe_load_step` on `[Utilities]`; the
`load_step_predictor` utility offers `mode = linear` (extrapolation) and
`mode = network` (a network u(x, lambda) fitted online to the corrected
fields on its own torch backend/device). Per step the history carries
`newton_time_s`, `predictor_time_s`, `predictor_rel_error` and the training
cost. Pipeline: a `[Predictor]` block of the pipeline deck is rendered into
the generated deck. First check launched on h40 (cold / linear / network,
`decks/cp_study/hybrid/`), compared by `scripts/hybrid_compare.py`.

**C2 design (next): the per-step "PINN as solver".** The same hooks carry
it; what is new is a physics-trained predictor and the executioner's
option to accept its proposal without a Newton correction:

* `load_step_predictor mode = physics`: the network u(x, lambda) is trained,
  per load step, on the strong-form balance of momentum evaluated at the
  quadrature points of the mesh with the crystal-plasticity state FROZEN
  there (the quadrature points are the collocation points, so the state
  needs no interpolation), with the per-point orientations and rotated
  stiffness, plus the Dirichlet term at the nodes and, optionally, the data
  term of the corrected previous steps. The residual needs a
  `compute_stacked` on the predictor's torch backend (the material and the
  orientation stacks rebuilt there from the same deck config; the state
  stacks copied from the FE callback each step) and the divergence of the
  first Piola-Kirchhoff stress by autograd with respect to the quadrature
  coordinates. The executioner hands the predictor what it needs through
  `LoadStepContext.extra`: quadrature-point coordinates, the material
  config, the current state stacks, the per-point orientations and the
  finite-strain flag.
* `cp_quasi_static corrector = newton | none`: with `none` the executioner
  accepts the proposal as the step's solution, evaluates the FE residual
  norm of it (one assembly) and records it beside the predictor's own
  loss, then advances the state exactly as for a Newton solution. That is
  the pure-PINN column: accuracy measured, not enforced. With `newton` the
  same physics-trained predictor is the hybrid's better initial guess.
* Outputs and bookkeeping are those of C1, so the three columns of the
  benchmark (cold FE, hybrid, PINN as solver) come from the same executioner
  and the same decks, differing only in `[Utilities]` and `corrector`.

**2026-10-12, C1 first check (docs/07 §10.5).** Data-trained network
predictor: proposal worse than the previous step (median 1.6x the
increment), 19% MORE Newton iterations than cold, 4.7 h of training over
256 steps (quadratic growth). Linear extrapolation: proposal within 0.03%
of the increment at no cost, about one Newton iteration per step saved
under an absolute convergence criterion (none under the default relative
test, which penalises warm starts; protocol fixed in docs/07 §10.0). The
hybrid column of the benchmark therefore uses the physics-trained
predictor, with linear extrapolation as the reference; the data-trained
mode stays available as a documented negative baseline.

**2026-10-12, C2 first check (docs/07 §10.6): negative.** The physics-trained
network with a soft Dirichlet term and no corrector produces a field that is
not in equilibrium once the prescribed values are pinned (residual equal to
the reaction, meaningless stresses), although its own residual loss is
small. Required before any benchmark of this route: hard Dirichlet
enforcement in the ansatz and the residual evaluated on the accepted
field. Deferred until the FE regeneration is complete.

**2026-10-12, C2 on the plate: strong-form collocation is not viable on a
30-grain polycrystal (docs/07 §10.7).** A network regressed to within 0.5%
(L2) of the FE displacement has an FE residual of 0.73 of the reaction and
a strong-form loss of 5.8 on the mu/L scale; the physics-trained predictor
plateaus at 1.2 and every proposal is rejected. The displacement kinks at
grain boundaries carry the equilibrium and a global smooth network cannot
represent them. Options for the author:
(a) **FE-residual training** (`mode = fem_residual`): train the network on
the discrete weak-form residual of its own nodal values, assembled by the
torch twin of the FE assembler at the frozen state -- the loss Newton
minimises, so a converged network IS an FE solution; the network's value
is then amortisation across the load factor (a better initial guess, or a
solution when `corrector = auto` accepts it), not a new discretisation;
(b) a per-grain (domain-decomposed) ansatz with interface conditions, the
PINN literature's answer to heterogeneous media, considerably more work;
(c) stop at the hybrid with linear extrapolation (-37% Newton iterations)
as the measured benefit of a predictor on this problem, and report the
strong-form PINN as a negative result with the diagnostic above.
