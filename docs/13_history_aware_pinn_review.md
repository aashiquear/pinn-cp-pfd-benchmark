# 13. History-aware neural routes for crystal plasticity: what the literature offers PINNACLE, and at what cost

Written 2026-10-28 for the author, before any framework change. Sources: the two
literature syntheses the author placed in `data/literature/` (`ScholarGPT_history_CP`,
2019-2024 foundations; `SciSpace_history_CP`, 2020-2026 with the recent CP work), read
against the PINNACLE branch `feat/pinn_benchmark` (HEAD 9a706ae; ADR 0028, 0029, 0030,
0033) and the measured costs of the Al2219 plate (report §8, §10, §11, §17, §19) and of the
Ni specimen (report §8, Table 10). Nothing in the framework was changed by this document
when it was written; §8 (added 2026-10-29) records what was implemented afterwards.

Both syntheses are machine-written reviews, not primary sources. Their equations and model
classes are standard; their numerical claims about individual papers are flagged by the
syntheses themselves as unverified, and I treat them that way. The primary references that
matter for us, in the order I would read them: Bonatti, Berisha and Mohr 2022 (CP-RNN
surrogate of CP-FFT, discretization-robust); Liu et al. 2023 (recurrent neural operator,
two to five learned internal variables reproduce a CP polycrystal's history); Maia et al.
2024 (physically recurrent network: real constitutive points embedded in the network carry
the history); He and Semnani 2024 (incremental neural CDE: stability of a learned state
update inside implicit FEM); Masi and Stefanou 2023 (evolution TANN: learn the rate
$\dot{\mathbf z}$, not the increment, so the network is not tied to one step size);
Flaschel et al. 2025 (convex generalized-standard-material networks keep the implicit
Newton integration); Wang et al. 2026 (Plasolver: a neural operator as a warm start that
halves the classical solver's iterations); Vlassis and Sun 2021 (Sobolev training: fit the
derivative, not only the value, when the network must feed a Newton solver).

## 1. The one thing the literature agrees on

A history-dependent material is not a map from the current strain to the stress. It is a
map from the current strain and a *state* to the stress, plus an evolution law for the
state. Every method that works keeps that structure: the state is either the physical
internal variables ($\mathbf F^p$, slip resistances, dislocation densities, backstress) or
a learned low-dimensional stand-in for them, and it is committed only when the global step
is accepted. A network that is asked for $\sigma(\mathbf x, t)$ or $\mathbf u(\mathbf x, \lambda)$
with no state is the "standard PINN without state evolution" that both reviews put at the
bottom of their tables.

PINNACLE already satisfies this: the state lives at the quadrature points, is advanced by
the backward-Euler update (ADR 0033), is evaluated from the committed state at every trial
(the "state-commit discipline" of the SciSpace review §7.5), and the network of the PINN
routes solves one step's equilibrium with the state frozen on the mesh (ADR 0028, point 4).
So the framework is on the right side of the literature's main line. The question the
author asks is a different one: *why does our pure PINN route cost 5 to 7 times the
finite-element route, and what, if anything, from this literature brings it down.*

## 2. Where the time goes today (measured)

The hardware check of the Ni report (§8, Table 10: L.1 h020, 792 elements, 6,336 quadrature
points, first 16 load steps of the carried deck, two threads) gives the anatomy.

| route | Newton iterations (16 steps) | Newton time per iteration | predictor time per step | time per step |
|---|---|---|---|---|
| cold finite element (torch CPU) | 113 | 7.1 s | -- | 50 s |
| PINN+FEM, linear extrapolation | 63 | 8.6 s | 0.0 s | 37 s |
| PINN as solver (`fem_residual`), network on GPU | 63 | 8.7 s | 199 s | 236 s |
| PINN as solver (`fem_residual`), network on CPU | 63 | 8.6 s | 302 s | 339 s |

Three facts come out of the per-step rows of those runs
(`results/ni_study/benchmark/hardware_L1_h020/solver_*_t2/solver/solver_history.csv`):

1. **The predictor's time is assembly time.** Of the 199 s (GPU) and 302 s (CPU) per step,
   98-99% is `predictor_assembly_time_s`: 35 (GPU) to 53 (CPU) assemblies of the FE residual
   and tangent per step, at 5.7 s each. The network's own forward and backward passes are
   negligible. The "GPU advantage" of this route is only that L-BFGS converged in fewer inner
   evaluations there.
2. **Each assembly is a full constitutive integration, ten times over.** With
   `state_update = implicit` the stress at a trial displacement is the backward-Euler
   solution at every quadrature point (a damped local Newton on the 12 slip increments), and
   the tangent is a finite-difference stencil of that stress: nine perturbations of the
   displacement gradient (`cp_field_callback.tangent`, `full_gradient`), each a complete
   re-integration of all 6,336 points. One residual-and-tangent assembly is therefore ten
   backward-Euler integrations of the whole specimen plus the element loops; the linear
   solve is small beside it. This is also the anatomy of a Newton iteration of the
   finite-element route, which is why the Newton iteration costs 7 to 9 s. (A profile of
   four load steps is in §2.1.)
3. **The network does not improve on its baseline.** The equilibrium error ratio of the
   linear-extrapolation baseline at the start of the path is 1.3 to 3.0 (the residual
   exceeds the reaction: at 0.1% strain per step the first steps are far from linear), and
   after 4 to 5 L-BFGS epochs it is 1.3 to 2.96: a change in the third digit. Not one of the
   16 steps was solved by the network; the Newton iterations (63) are exactly those of the
   linear baseline. On the Al2219 plate (report §8.2, Table 8) the same route did better,
   132 of 256 steps solved without Newton with grain-aware inputs, but at 5 to 6 times the
   finite-element cost, because the per-step training cost 10 to 20 s against 1 to 3 s of
   Newton.

So "the predictor stage takes all the cost" has a precise cause. The loss that makes the
network verifiable (the assembled residual, ADR 0029) is evaluated by the most expensive
object in the framework, the implicit crystal-plasticity integration with a finite-difference
tangent, once per epoch; and the per-step field network has nothing cheap to learn beyond
the linear extrapolation it starts from. Neither of these is the history problem the
literature discusses. The history is already handled correctly; it is the *price of each
look at the physics* that is wrong, and the *object the network is asked to learn* that is
poorly chosen.

### 2.1 Profile of a Newton iteration

cProfile of four load steps of the carried L.1 h020 deck (1e-3 grip strain per step, torch
CPU, two threads, local workstation; 28 global Newton iterations, 32 tangent assemblies, 74
stress assemblies; the profiler roughly doubles the wall time, 129 s per step against 50 s
unprofiled, but the shares are what matter):

| part of the mechanics solve | time | share | what it is |
|---|---|---|---|
| backward-Euler integrations, all | 512 s | 99% | 1,322 whole-specimen integrations (`_backward_euler_stacks`) |
| of which the tangent stencil | 405 s | 79% | 32 tangents x 10 integrations (one unperturbed, nine perturbed) |
| of which the stress at a trial | 98 s | 19% | 74 residual assemblies, 62 s of them in the global line search's backtracking |
| element loops, linear solve, rest | 6 s | 1% | |

Inside one whole-specimen integration (0.39 s each): the local Newton of the 12 slip
increments runs about 11 iterations and 15 residual evaluations per integration (14,377
batched 12 x 12 solves and 19,380 residual evaluations over 1,322 integrations), the
batch iterating until its slowest point converges; the per-point `numpy.linalg.solve`
alone is 25% of the run, the residual evaluation (`end_state` + `compute_stacked`) 31%.

Two conclusions. Four fifths of a Newton iteration is the finite-difference tangent,
which P1 (§4.1) removes; and the local solve that remains spends eleven iterations where
a good start would need two or three, which P3 (§4.3) targets. The global line search's
backtracking (12%) is a third, smaller item: a predictor that lands inside the Newton
basin removes it.

## 3. The literature's five families, mapped onto PINNACLE

| family (both reviews) | what is learned | PINNACLE today | relevance to our bottleneck |
|---|---|---|---|
| Recurrent constitutive surrogate (CP-RNN, GRU/LSTM) | stress from a hidden state | none | replaces the local integration; needs its own tangent; drift risk |
| Thermodynamics-informed state models (TANN, eTANN, GSM networks, NN-EVP) | free energy and state evolution | none | correct structure; the GSM networks keep an implicit Newton, so no speed gain by themselves |
| Physics-informed state / hybrid implicit (Strategy B/C of the ScholarGPT review) | the slip increments or a correction; the implicit solve kept, tangent by implicit differentiation | the backward-Euler solve exists; the tangent is finite-difference | **the direct answer to item 2 of §2** |
| Full-field solver / neural operator (Plasolver, PINN as solver) | the field of a step or of a family of problems | `fem_residual` per step (ADR 0029); surrogate (ADR 0030) | the per-step version is what we measured; the operator version is our surrogate route |
| Physically recurrent / embedded material points (PRNN) | an encoder-decoder around real constitutive points | none | the right idea for the specimen surrogate's output structure (§6.3) |

Two statements of the reviews deserve emphasis because they decide what is worth doing:

- *Amdahl.* "Even a 1,000x local speedup produces only about a 10x full-solver speedup if
  90% of the baseline was spent in the accelerated computation." In our Newton iteration
  the constitutive integration (with its stencil) is the 90%; the lever that matters is the
  one that cuts the integrations per iteration, not the one that makes a network faster.
- *Break-even.* A trained model pays when it is queried many times. A network trained per
  load step and per trajectory (our PINN routes) never reaches that regime; the surrogate
  does (Al2219 report §11: 14 queries). The literature does not contain a per-step PINN
  that beats an implicit finite-element step on a plasticity problem; the gains it reports
  are at the material point (replacing an expensive local solve) or across many problems.

## 4. What can be done, with cost and gain

Each item names the code it touches, the development effort I estimate, the measured or
derived gain, and the risk. "Dev" is my own time including tests, ADR, syntax docs and
CHANGELOG under the PINNACLE rules.

### 4.1 P1. Consistent backward-Euler tangent by implicit differentiation (all routes)

What: replace the nine-perturbation finite-difference stencil of `cp_field_callback.tangent`
(implicit branch) by the algorithmic tangent of the backward-Euler update. The local solve
already has the semi-analytic Jacobian $\mathbf J = \mathbf I + \Delta t\,(n|\dot\gamma|/|\tau|)\,\mathbf K + \Delta t\,(n\dot\gamma/\tau_c)\,h\,\mathrm{sign}\,\Delta\gamma$
(`_backward_euler_stacks_graphfree`). The consistent tangent is

$$
\frac{d\boldsymbol\sigma}{d\mathbf H} = \frac{\partial\boldsymbol\sigma}{\partial\mathbf H}\Big|_{\Delta\boldsymbol\gamma}
+ \frac{\partial\boldsymbol\sigma}{\partial\Delta\boldsymbol\gamma}\,\frac{d\Delta\boldsymbol\gamma}{d\mathbf H},
\qquad
\frac{d\Delta\boldsymbol\gamma}{d\mathbf H} = -\,\mathbf J^{-1}\frac{\partial\mathbf R}{\partial\mathbf H},
$$

one $12 \times 12$ solve with nine right-hand sides per point, batched over the points
(Strategy C, "hybrid implicit model", of the ScholarGPT review §7; the same implicit-function
construction the GSM networks of Flaschel et al. use). The ingredients $\partial\mathbf R/\partial\mathbf H$
(through the elastic Schmid coupling) and $\partial\boldsymbol\sigma/\partial\Delta\boldsymbol\gamma$
(through $\mathbf F^p$) are the quantities $\mathbf K$ already encodes; what is new is
assembling them into the chain rule and the finite-strain bookkeeping of the exponential-free
$\mathbf F^p$ update.

Gain: a Newton iteration goes from ten constitutive integrations to one plus a batched
linear-algebra pass; from the ADR 0033 timings (0.15 s per stress and 1.7 s per tangent on
1,760 points) and the profile of §2.1, that is a 3 to 5 times cheaper iteration on every
route (cold FE, PINN+FEM, the solver route's assemblies, the surrogate's data runs). The
tangent becomes exact instead of a 1e-6 difference of a nonsmooth local solve, which is
the cause of the Newton floor at 0.16 N on the multi-slip orientation (group 9 of the
campaign, report §9.3, 210 cutbacks and 46 h on 9.1): quadratic convergence to the
tolerance, no floor, fewer cutbacks. For the campaign that is the difference between 15 to
46 h per specimen and an estimated 4 to 12 h, and between group 9 solving or not.

Cost: 3 to 5 days. Risk: moderate, in the finite-strain terms (the stress is the Cauchy
stress of $\mathbf F^e = \mathbf F\,(\mathbf F^p)^{-1}$ with the lattice rotation; every
term must be differentiated consistently or the quadratic convergence is lost, which the
existing finite-difference stencil makes easy to test against). Validation: the tangent
against central differences at a converged point (test), the Newton iteration count on
the L.1 deck against the stencil's, and the campaign's own cutback counts. No change to
any result: the tangent only drives the iteration.

### 4.2 P2. Cheaper assemblies in the residual bridge (PINN as solver only)

What: `FEMResidualBridge.residual` assembles residual *and* tangent at every epoch because
the backward pass needs $\mathbf K^\mathsf T$. Add `tangent_refresh = n` (reuse the tangent
of epoch 0 for $n$ epochs; the gradient is then that of a modified Newton, inexact but a
descent direction), and prefer Adam to L-BFGS when the tangent is stale. With P1 in place
the refresh matters less; without it, it turns a ten-integration epoch into a
one-integration epoch.

Gain: 5 to 8 times cheaper training per step for the solver route; combined with P1, the
solver route's 236 s per step would become roughly 30 to 50 s, i.e. parity with the cold
finite-element step rather than 5 times it. Cost: half a day. Risk: low; the verification
(the accepted field's residual below `corrector_tol`, ADR 0028 guard) is unchanged.

What it does not do: make the network solve steps. On L.1 the network lowered the
baseline's residual ratio by 0.1% in 4 epochs; on the plate it needed grain-aware inputs
and 50 to 60 assemblies per step to solve half the steps. Cheaper epochs make the route
cheaper, not better.

### 4.3 P3. Neural warm start of the local backward-Euler solve (the material point)

What: the literature's Strategy B in its safest form. A small network per material
(inputs: the trial elastic strain relative to the committed state, the committed slip
resistances and densities, $\Delta t$; outputs: the 12 slip increments) predicts the start
of the local Newton in `_backward_euler_stacks_graphfree`; the Newton then converges from
the prediction, so the result is the exact backward-Euler solution whatever the network
does, and a bad prediction costs only the iterations it would have cost anyway (Plasolver's
"warm start" idea, at the material point instead of the field). Training data are free: every
backward-Euler solve of every campaign run is a labelled sample (trial in, converged
increments out); a writer that logs them and an offline trainer are the new pieces. Trial
evaluations never commit anything, so the state discipline is untouched.

Gain: the local Newton's iteration count (not instrumented today; the profile of §2.1 gives
the share of the stress evaluation in an iteration). If the local solve takes five to eight
iterations and the prediction brings it to one or two, the stress evaluation becomes 3 to
4 times cheaper. After P1 the stress evaluation is the dominant term of a Newton iteration,
so P3 is worth 2 to 3 times on top of P1; before P1 it is worth little (the stencil
dominates). Cost: 3 to 5 days (sample writer, dataset format, a trainer reusing the
surrogate's network code, a deck switch `local_warm_start = <weights>`, tests). Risk: low;
the fallback is the current solver. The Jacobi phase and the sub-increment ladder of ADR
0033 stay as they are.

### 4.4 P4. Replace the local solve by a learned state update (CP-RNN / RNO / eTANN)

What: the family with the largest reported local speedups (10 to 1,000 times at the
material point, Maia et al. and Bonatti et al.; the ScholarGPT review's recommended
architecture §7.2). A state network predicts the end-of-step stress and state directly, with
the physical variables kept as the state (Strategy B) or a learned low-dimensional state
(RNO), with a fallback to the exact solve when an error indicator fires, and a tangent by
automatic differentiation of the network (or Sobolev-trained).

Gain: bounded by Amdahl. After P1 and P3 the constitutive share of a Newton iteration is
perhaps half; replacing it entirely gains at most another factor of two on the iteration,
and the element loops and the linear solve remain. Cost: 2 to 4 weeks (training over
orientations, rates and paths; drift tests over 400-step histories; the tangent; the
fallback logic; ADR), plus the loss of exactness the campaign currently has. Risk: high
(the reviews' own warning: hidden-state drift, inconsistent tangents, "a 2% stress error
at the material point may still fail inside implicit FEM"). Not recommended for this study;
it is a study of its own.

### 4.5 P5. Should the pure-PINN cycle be revamped?

The honest answer from the measurements and the literature: no. The per-step field network
is the one object in the taxonomy that cannot reach break-even, because it is trained once
per step and per trajectory, and the finest instrument we have (the assembled residual) is
the most expensive loss there is. The literature's gains come from the other end: the
material point (P3, P4) or the parameter space (the surrogate). What I would still do for
the one representative specimen of the campaign, because it is cheap and makes the route's
record fair: P1 and P2 (its cost drops from 5 times to about parity with FE), a quadratic
baseline (§4.6), and state-aware node features (the accumulated slip and the phase field $d$
at each node from the callback's state, the single-crystal analogue of the grain-aware inputs
that made the plate network solve half its steps: a network of $(\mathbf x, \lambda)$ cannot
kink where the neck and the crack form, a network that also sees $d$ and $\bar\gamma$ can).
Cost of the features: one day (a `node_features = 'state'` option of `pinn/node_features.py`
reading the callback's committed stacks); gain unknown, to be measured on the representative
specimen; risk none (a deck switch).

### 4.6 P6. PINN+FEM: three cheap improvements to the predictor route

- **Quadratic extrapolation baseline** (`baseline = quadratic`: the parabola through the
  last three accepted fields). The linear extrapolation already halves the Newton
  iterations; on a smooth path the quadratic one should take another 20 to 40%
  (measured on the plate: the linear proposal was within 0.03% of the increment beyond 1%
  strain, so the remaining error is second order). Cost: hours. Risk: none (the guard
  rejects a bad proposal).
- **Tangent reuse across Newton iterations when the proposal is good** (`tangent_every = n`
  on the executioner, a modified Newton for the first $n$ iterations of a step whose
  predictor residual is small). Today the backward-Euler tangent is recomputed at every
  iteration because it is "cheap" relative to the sub-stepper's; it is not cheap in
  absolute terms (nine integrations). With P1 this item becomes moot. Without it, it halves
  the iteration cost at the price of one or two extra iterations. Cost: half a day.
- **Adaptive step growth** when the predictor's relative error is small for several steps
  (the path of the specimen is 400 steps of 1e-3 grip strain because the explicit update
  needed it; ADR 0033 showed backward Euler within 2% of a 2,000-step reference at 1e-3,
  so a 2e-3 step may well be within the same tolerance away from the fillets and the
  crack). Cost: one day plus a convergence check like the one of ADR 0033 (one run).
  Gain: up to half the steps on the smooth part of the path. Risk: the failure strain
  resolution (the calibration tolerance is two load steps); keep 1e-3 near failure.

### 4.7 P7. The specimen surrogate: structure it with the physics it can get for free

The planned surrogate (notebook §8: eight inputs, three Euler angles and five station
thicknesses, to the curve, the failure strain and the location; 30 Latin-hypercube samples,
15 real specimens held out) is a parametric map in the sense of ADR 0030, and the Al2219
numbers say its accuracy will be set by the data (2 to 3% on the curve with 18 samples in a
four-dimensional box; here eight dimensions with 30). The literature suggests four changes
that cost little and attack exactly that:

1. **Invariant orientation inputs.** Euler angles are periodic and, under the cubic
   symmetry, every orientation has 24 equivalent angle triples; a network of raw angles has
   to learn a discontinuous map. Replace them by the Schmid factors of the 12 systems for the
   loading axis (sorted, or the largest three and their count near the maximum) or by the
   symmetric-equivalent angle in the fundamental zone. This is the "equivariant / invariant
   feature construction" of the ScholarGPT architecture §7.2. Cost: hours; the Schmid
   factors are one line from the orientation code.
2. **Physics-embedded output (the PRNN idea).** For any orientation the single-crystal
   material-point curve $\sigma_{\mathrm{MP}}(\varepsilon)$ of the calibrated CP set costs
   seconds (`scripts/taylor_surrogate.py` carries the one-point integrator). Make the
   network learn the *ratio* of the specimen curve to the material-point curve, or feed the
   material-point curve as an input feature. The network then learns only the specimen
   effects (the gage taper, necking, the fillets, the damage), which are smooth in the
   inputs, and the hardening of the orientation, which is not, comes from the physics. Cost:
   one day (a feature hook in the surrogate trainer). Gain: on the Al2219 plate the curve
   surrogate's error came from the elastic knee and the hardening amplitude, both of which
   this removes.
3. **Separate heads, shared trunk.** The failure strain and the location are scalars with
   their own targets; predict them with their own output layers from the shared trunk, with
   the failure strain also recoverable from the curve head (the first load drop) as a
   consistency check. The Al2219 coupled surrogate (report §19) missed one held-out failure
   because the curve head alone placed the drop 0.5% early; a scalar head does not have
   that failure mode.
4. **Multi-fidelity data.** The data cost dominates (30 h020 runs at 15 to 46 h each at
   present; 4 to 12 h after P1). ~~The curve to failure is nearly mesh-insensitive (notebook
   §3: h040 within 1% of h020 on the curve) while the location needs h020.~~ **Correction
   (2026-10-29):** this premise was wrong; notebook §5.2 measured h040 6-20% above h020, and
   the CP+PFD check of §8.5 adds a 5.7% early failure and a location 0.09 off on h040. The
   additive correction exists in PINNACLE, but the campaign runs h020 only. Run the 30
   samples at h040 for the curve and failure strain, and a subset of 8 to 10 at h020 for the
   location head (or correct the h040 location by the h020 subset). Cost: one day of tooling;
   gain: the data campaign shrinks from weeks to days.

Also worth doing because training takes seconds: an ensemble of five networks (different
seeds) to flag a held-out specimen whose inputs lie outside the training cloud (the
reviews' "out-of-domain detector").

## 5. Recommended order and what it buys the campaign

| order | item | dev | gain | affects |
|---|---|---|---|---|
| 1 | P1 consistent tangent | 3-5 d | 3-5x per Newton iteration; no Newton floor; group 9 solvable | every route, the surrogate data, the W_c calibrations |
| 2 | P6a quadratic baseline | hours | -20 to -40% Newton iterations | PINN+FEM |
| 3 | P7.1, P7.2, P7.3 surrogate structure | 2 d | accuracy at fixed data | surrogate |
| 4 | P7.4 multi-fidelity design | 1 d | data cost /3 to /5 | surrogate |
| 5 | P2 bridge tangent refresh + P5 state features | 1.5 d | solver route at ~FE parity; a fair record | PINN as solver |
| 6 | P3 local warm start | 3-5 d | 2-3x on the stress evaluation after P1 | every route |
| -- | P4 learned state update | weeks | <= 2x after 1 and 6; exactness lost | not this study |
| -- | P6c adaptive steps | 1 d + 1 run | up to -50% steps | after a convergence check |

With items 1 and 2 the finite-element specimen goes from 15-46 h to an estimated 3-10 h,
the PINN+FEM specimen to 2-7 h, and the solver route to about the finite-element cost; the
30-sample surrogate data from weeks to days. Items 3 and 4 decide whether the surrogate
reaches the few-percent accuracy the parity plot needs. Item 6 is the first step toward
the literature's material-point acceleration without giving up exactness; item P4 is the
full version and belongs to a later study.

What the literature does *not* promise and I would not claim: a per-step PINN that beats
implicit FEM on this problem; a learned constitutive update that is exact; a surrogate
accurate outside its design box.

## 6. If approved, what each change needs under the PINNACLE rules

P1: ADR 0033 addendum (or ADR 0034) with the derivation; `tests/unit/materials/
test_cp_implicit_update.py` gains the tangent-versus-central-difference test and a
Newton-iteration-count test on the single element; `state_update` docs name the tangent;
CHANGELOG. P2/P5/P6: parameters on `load_step_predictor` and `cp_quasi_static` with syntax
regeneration (`pinnacle syntax md`), schema, tests in `test_load_step_predictor.py`. P3:
new component (`materials/cp_local_warm_start.py` or a parameter of the callback), ADR,
tests, an example deck. P7: benchmark-side scripts only (no framework change), unless the
Schmid-factor feature and the multi-head target go into `parametric_surrogate` (then ADR
0030 addendum).

## 7. Appendix: the profile run

`python3 -m cProfile -o profile.out -m pinnacle run pfd.i` on a copy of the hardware
check's cold-FE deck (`results/ni_study/benchmark/hardware_L1_h020/fe_torch_cpu_t2/pfd.i`)
with `n_steps = 4`, `t_end = 0.08`, the pull expression scaled to the same 1e-3 grip strain
per step, no checkpoint; PINNACLE 9a706ae, torch 2.12.1 CPU, `OMP_NUM_THREADS=2`, local
workstation, 2026-10-28 10:09-10:17. Steps accepted in 11 / 5 / 6 / 6 Newton iterations
(235 / 83 / 100 / 98 s under the profiler). Not a ledger run (no result is kept; the
profile output is in the session scratchpad and the 45-line summary is reproduced in §2.1).

## 8. Status after implementation (2026-10-29)

The author approved items 1 to 6 of the §5 order on 2026-10-28 (a grouped version of the
list: P1; P6a; P7.1-P7.4; P2 + P5 state features; P3). Implemented in PINNACLE
`feat/pinn_benchmark` at 1aa1ff5 (ADR 0034, ADR 0028 and 0029 addenda, CHANGELOG,
regenerated syntax docs, tests); the specimen-level measurement is the notebook's §8.4
and the ledger rows `pinnacle-p1-L1-stencil` / `pinnacle-p1-L1-consistent`.

| item | status | where | what was measured |
|---|---|---|---|
| P1 consistent tangent | **done** | `implicit_tangent = consistent` (default under torch); `cp_field_callback._consistent_tangent`; ADR 0034 | 1e-7 of the largest entry against central differences (stencil 1e-4); 0.04 s against 0.12 s per point; L.1 four-step deck 706 s to 93 s with the warm start below |
| P2 bridge tangent refresh | **done** | `tangent_every` on `load_step_predictor` (`FEMResidualBridge(tangent_every=...)`), reported as `_tangent_assemblies`; ADR 0029 addendum | unit test only (the residual still falls with a stale tangent); not yet run on a specimen, the solver route is paused |
| P3 warm start of the local Newton | **done in its exact form, not the neural one** | the local Newton starts from the point's last converged increments (or the summed ladder sub-increments); ADR 0034 addendum | essential, not optional: without it the one-increment solve fails on nearly every trial, the ladder serves and P1 cannot serve. The learned predictor of §4.3 (sample writer, trainer, deck switch) is not built; the exact warm start removes most of its case |
| P4 learned state update | not built (author 2026-10-29); the development path is documented | PINNACLE `docs/research/learned_state_update_path.md` (what it is, the Amdahl bound after ADR 0034, six steps, gating tests, decision rule) | |
| P5 state-aware node features | **done** | `node_features = 'plastic_strain damage'` (`fem_residual` only; qp to element mean to node mean, scaled); ADR 0029 addendum | unit test only; gain on the representative specimen not measured (route paused) |
| P6a quadratic baseline | **done** | `baseline = quadratic`, `mode = quadratic`; ADR 0028 addendum | exact on a quadratic path; on the Ni single element worse than linear on the step that straddles the yield knee, then better, 20x by the tenth step. The 20-40% Newton-iteration gain on a specimen is not yet measured |
| P6b tangent reuse across Newton iterations | not done | | moot after P1 (the tangent costs one integration) |
| P6c adaptive step growth | **done** (2026-10-29) | `step_growth_max` (+ `_after`, `_factor`, `_tol`, `_max_newton`, `_damage`) on `staggered` / `monolithic`; windows of nominal steps with a damage look-ahead and a drop back near failure; ADR 0028 addendum | on the 1-D bar with the implicit update: a 48-step path in 15 windows within 3 % of the nominal; 2x windows of a 12-step path differ by 2.4 %, 4x by 7 % (backward Euler is first order), so it stays a deck option (default off). The specimen convergence check (one run at `step_growth_max = 2`) is the author's and awaits the campaign |
| P7.1-P7.4 surrogate structure and multi-fidelity data | **done** (2026-10-29, framework side) | ADR 0030 addendum on `parametric_surrogate`: `orientation_features = schmid` (`materials/schmid.py`, invariant under the cubic group), `baseline = material_point` (`materials/material_point.py`, the implicit one-point update of the FE route), `scalars` with heads on a shared trunk (`pinn/networks/trunk_heads.py`, `drop_fraction` consistency check), `hf_scalars` / `hf_overrides` / `hf_samples` (additive correction), `ensemble` / `ensemble_flag`; external `target.npz` with `e`, `s`, `s_hf` | mechanism check on the 1-D CP bar over the three Euler angles and rho0 (24 samples, 6 held out): held-out curve 5.4 % mean / 13.6 % max, scalars 6 %, corrected high-fidelity scalar 10 %, spread 6 %, no flag; 26 s for 24 baseline curves, 16 s training. The specimen gain is measured when the surrogate stage runs (on hold); the benchmark-side writer of the specimens' `target.npz` is the next script |

Two things the implementation taught that the review did not know: the warm start is a
precondition of P1 on this material (§4.3 described it as a 2-3x extra); and the first
step of a run differs by 1.2% in the reaction between the old and the new solver (the old
commits were served by the sub-increment ladder, a different discretization of the same
step), after which they agree to 0.05%. Any rerun of the campaign with the new solver
therefore reproduces the old results, not bit for bit, but within the calibration tolerance.

Update 2026-10-29 (evening): the author approved P7 (all four parts) and P6c for the
framework and asked for P4 to be documented, not built. Implemented in PINNACLE
`feat/pinn_benchmark` (ADR 0030 and 0028 addenda, `docs/research/learned_state_update_path.md`,
CHANGELOG); the table above carries the status. No benchmark run was made (the author's hold).

Correction 2026-10-29 (later): the P7.4 premise (h040 within 1% of h020) does not hold. Measured
on L.1 with the carried set (notebook §8.5): h040 is +5.8% to +18.5% on the curve, fails at 0.438
against 0.464 and 0.09 of the gage further along. The mechanism stays in the framework; the Ni
surrogate's data are h020 runs.
