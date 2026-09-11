# PINNACLE capability audit for the CP–PFD benchmark study

Phase 0, Framework track. Audited revision: `feat/pinn_benchmark` branched from
`master` at `16ff0e6`. Component inventory captured at
`env/pinnacle_syntax_dump.json` (249 registered components).

Method: `pinnacle syntax dump --format json`; `docs/user_guide.md` §1.1 (physics ×
solver matrix), §1.2 (capability index), §1.3 (coupled hybrid decks); the ADR
series (0006, 0015, 0016, 0017, 0018, 0019); and direct reading of
`materials/crystal_plasticity.py`, `materials/hardening.py`,
`materials/dislocation_density.py`, `physics/slip_kinetics.py`,
`system/damage_kernels.py`, `mesh/generators.py`, and `postprocess/`.

---

## 0. Gap status at a glance

The authoritative status of every gap, with the commit that closed it. Severity
codes: **A** blocks a faithful reproduction, **B** blocks part of the study,
**R** required for a reported quantity, **D** documentation only.

| Gap | Sev | Status | Commit | Subject |
|---|---|---|---|---|
| G1a | A | **fixed** | `aac832b` | Driving force was evaluated on the *total* strain, not the elastic strain |
| G1b | A | **fixed** | `aac832b`, `1f00312` | No stored plastic work in the driving force; hybrid route needed a second fix |
| G2 | B | **fixed** | `1c85ab2` | Taylor prefactor α₀ absent, fixed at 1.0 against the papers' 0.4 |
| G3 | B | **fixed** | `1c85ab2` | Declared units of `k1`/`k2` contradicted the implemented law |
| G4 | B | **fixed** | `84bebc7` | Isochoric loading: needed no new component, but exposed G13 |
| G13 | A | **fixed** | `84bebc7` | Per-BC time functions were silently ignored by `cp_quasi_static` |
| G5 | B | **fixed** | study `scripts/make_void_meshes.py` | No conforming spherical-void mesh, no local grading — met from hexaform, no framework change |
| G14 | A | **fixed** | `4e48301` | Exodus reader silently discarded every file-provided node set |
| G15 | A | **fixed** | `43dd6fd`, ADR 0023 | `staggered` never incremented the load and reset CP state each sub-solve |
| G19 | A | **fixed** | (this session) | `cp_quasi_static`'s reported stress was evaluated after the state advance, double-counting plasticity |
| G24 | A | **fixed** | `981e2dd` | Damage never fed back into the mechanics stress: three wiring faults, both CP and J2 routes |
| G27 | A | **fixed** | `5216640` | J2 degradation multiplied the whole stress instead of splitting tension from compression |
| G26 | A | **fixed** | `e74e5fa` | Degradation applied to the Cauchy stress post-hoc instead of to the 2nd Piola-Kirchhoff stress per Eq. (29) |
| G25 | A | open, BLOCKING | — | A load step abandoned at the cutback floor returns a field violating its own Dirichlet BC; study runs now refuse to write results |
| G22 | A | **fixed** | `aef2302` | `_snapshot_material_state` returned None on every call: the alternation state restore was a silent no-op |
| G23 | B | open | — | Coupled restart implemented but NOT bit-exact against an uninterrupted run |
| G21 | R | **fixed** | (this session) | Sparse DIRECT solve scales n^2.91 and made fine meshes unreachable; CG is exact and ~linear |
| G20 | A | **fixed** | `1f0834d` | Load-path resolution: threshold measured, enforced in the harness, guarded in the framework |
| G17 | A | **fixed** | `7156159` | Newton hardcoded a DENSE tangent; 3 of 5 mesh levels could not allocate |
| G18 | A | **fixed** | `7156159` | CP tangent finite-differenced one quadrature point at a time |
| G16 | A | **fixed** | `4340932`, ADR 0024 | `staggered` had no alternation loop: `max_iterations`/`tolerance` were dead, so coupling was lagged |
| G6 | B | open | — | No dog-bone mesh generator |
| G7 | R | **fixed** | `1f00312` | No Miller-direction → Bunge-Euler utility; closes C.6 |
| G8 | R | open | — | Missing study metrics (void volume, fracture plane, Hausdorff, …) |
| G9 | D | noted | — | Exponent and reference-slip-rate mapping, no code change |
| G10 | D | noted | — | Degradation specified in two places; confirm they agree |
| G11 | B | open | — | `backend`/`device` are `[Global]`-only, so a hybrid can't split hardware |
| G12 | A | **fixed** | `b49fa70` | Two contraction errors in the batched CP path |

Twenty-one of twenty-eight closed. **Five remain open, and G25 (severity A)
blocks Phase 1.**

RETRACTION. A previous revision of this section claimed G25 "no longer
blocks the fracture study" once G26 was fixed. That was wrong, and the
evidence for it was badly chosen: it rested on a Gc = 2.0 MPa.um case where
the specimen severed within five load steps, so the run reached a terminal
state and *looked* complete. Re-run at Gc = 30 MPa.um, inside the physical
band, with 60 load steps and the coarsest mesh, the solve fails at step 12
while still carrying 144.8 MPa, 60.2% of its 192.4 MPa peak, with only 12%
of the domain damaged -- that is the onset of crack propagation, not the end
of it. G26 removed a common trigger and let runs get much further; it did
not make the coupled solve traverse propagation.

The other four open items do not block: **G6** and **G11** are development
work, **G23** is coupled restart (implemented but not bit-exact), and
**G8** is the study-metrics register. G9 and G10 are documentation notes
needing no code. The one genuinely open research item, R5 (the
history-aware CP state surrogate), is scoped as its own track rather than
as a gap.

**Both severity-A blockers are now closed.** G15 made the coupled solution
path-correct in load and in constitutive state, which is what made the study's
primary metrics *defined* rather than merely imprecise; G16 makes the solution at
each load level *converged* rather than lagged. A benchmark reference needs
both, and the discretization and solver-convergence step is unblocked.

One consequence to carry into every later comparison: the converged scheme
predicts **lower** stress than the lagged one — `j2_pfd_void_in_cube` fell from
80.90 to 74.80 MPa, about 7.5 % — because a lagged coupling solves mechanics
against a damage field that has not yet responded to the current load, so it
under-damages and over-predicts load capacity. Any coupled number produced
before commit `4340932` is an upper bound on stress.

Note that the register has **grown** as work proceeded, from eleven entries to
nineteen: G12 through G18 were each found by verifying something else or by
profiling a real run, and none was visible from reading the documentation or
the component index. That is worth stating plainly, because it bears on how much of the
framework should be trusted without a check. Every defect found so far was
**silent** — none raised an error — and five were in code paths the
documentation describes as working. Four of them (G13, G14, G15, G16) surfaced
only because a verification was written to measure the physical quantity rather
than to confirm that a run completed: the isochoric checker measures volume from
deformed node positions, and the staggered probe compares displacement fields
across step counts. A check that merely asserted "the run converged" would have
passed in every one of those cases.

Three of the defects share one signature, and it is worth naming as a review
heuristic for the rest of the framework: **a capability declared in a schema,
implemented in a sibling component, and unwired in the one under test.**
Per-BC time functions (G13) were honored by four executioners and dropped by
`cp_quasi_static`; file-provided node sets (G14) were preserved by the
derivation step but never read by the reader; and the load-factor window (G15)
was honored by `fem_quasi_static` and ignored by `cp_quasi_static`, which
declared it in its own parameter list. In each case the deck validated, the run
completed, and the declared parameter did nothing.

Each claimed status is checkable against the branch, e.g.
`grep -c "def direction_to_bunge" src/pinnacle/materials/orientations.py` for
G7. The entries below carry the evidence.

## 1. Solver routes actually available for coupled CP–PFD

| Route | Selection | Status | Use in this study |
|---|---|---|---|
| **R1 — all-FEM** | `[Executioner] type = staggered` with `[./mechanics] cp_quasi_static` + `[./damage] pfd_steady` | Available; bundled decks under `examples/fem/phase_field_damage/ductile_anisotropic_cp/` and `examples/applications/rezwan2025_void_cp_pfd/` | Converged reference. Every comparison is measured against this. |
| **R2 — hybrid, damage as a PINN field** | `pinn_fem_hybrid_coupled_physics`, `mode = cp_pfd`; `u` family `fem`, `d` family `pinn` | Available; 1-D/2-D/3-D decks under `examples/pinn_fem/phase_field_damage/ductile_anisotropic_cp/`. The 1-D deck reproduces the all-FEM staggered reference to ~1e-11 | Primary neural route. |
| **R3 — hybrid, PINN predictor → FEM corrector** | `pinn_fem_hybrid`, `mode = warm_start \| predictor_corrector \| preconditioner` | Available; `void_field_hybrid_3d.i` | Same accuracy as R1 by construction (corrector drives to solver tolerance); the measurable quantity is cost, not error. |
| **R4 — neural constitutive surrogate (0-D)** | `[Executioner] type = cp_constitutive` with a `pinn`-family variable; `void_pinn_constitutive_3d.i` | Available | Constitutive-level neural route. |
| **R5 — history-aware CP state surrogate** | not present | **To be developed** (Phase 3, Neural CP track) | Supplies the genuinely neural CP column. |
| **Pure-PINN coupled CP–PFD field solve** | `pinn_two_field` | **Structurally unavailable**, 1-D elasticity only, raises on non-1-D input | Not a study column. Cited as ADR 0017, not reported as a failure. |

ADR 0017 settles the pure-PINN cell: a coordinates-in network has no memory, and
the CP state at a material point is a functional of the entire load path, so two
points at the same `(x, t)` can be in different states. The ADR lists three
prerequisites for a future neural CP route; requirement 2 — a mesh-based state
update exposing `(F, state, Δt) → state_after_increment` — is exactly what the
R5 track builds. Supporting evidence that the slot is anticipated: the
`classical_fallback_threshold` parameter on
`pinn_fem_hybrid_coupled_physics` is documented as *"Until the integrated
surrogate ships, this value is always 1.0 and the warning is suppressed"*
(ADR 0016).

---

## 2. Gap register

Severity: **B** blocks a study deliverable · **R** required for fidelity to the
source papers · **D** desirable.

### G1 — FIXED — The fracture driving force carries no plastic-work contribution

*Original finding below; on tracing the coupling channel this turned out to be
two separable defects, and both are now closed. See **G1a / G1b** later in this
section for the split, the resolution (`aac832b`, `1f00312`), and the measured
effect.*

`system/damage_kernels.py` implements the history field as `H ← max(H, ψ⁺)`,
where `ψ⁺` is the positive part of the **elastic** strain energy from the
mechanics kernel (lines 267, 369–373, 682–704). A tree-wide search for
`plastic_work` returns no hits anywhere in `src/pinnacle/`.

Both source papers drive damage with elastic *and* stored plastic energy,
`H = ψ⁺_el + β W_p` with `β = 1.0`. Without the `W_p` term the coupled model
cannot reproduce ductile CP–PFD behavior: damage is driven only by recoverable
energy, so slip-induced localization does not feed fracture.

This is the most consequential gap in the audit. It is also a designated
sensitivity variable — the mesh-sensitivity plan sweeps `β` over
`0, 0.25, 0.5, 0.75, 1.0` and hypothesizes that `β` controls ductile versus
brittle-like response — so the parameter must exist before Phase 5 can run.

**Action.** Accumulate plastic work in the CP material state; expose a
`plastic_work_prefactor` (`β`) on `pfd_driving_force`; route `W_p` through both
the `staggered` coupler and the `cp_pfd` hybrid strategy. Register, document,
test, and record an ADR.

### G2 — FIXED — The Taylor hardening prefactor α₀ is not implemented (`1c85ab2`)

`materials/hardening.py:60-65`:

```
def critical_stress(self, rho):
    """Return τ_c^alpha = τ_c0 + μ b sqrt(Σ_beta a^{αβ} ρ^β)."""
    sum_a = ops.einsum(self.backend, "ab,b->a", self.a, rho)
    return self.tau_c0 + self.mu * self.b * ops.sqrt(..., sum_a, 0.0)
```

`crystal_plasticity.py:404` inlines the same expression on the batch path. There
is no `α₀` factor, and no such parameter is declared in the registry — the
implementation is fixed at `α₀ ≡ 1`. Both papers use `α₀ = 0.4`.

Consequence for Al2219 at the published parameter set
(`μ = 25.7 GPa`, `b = 0.000286 µm`, `ρ₀_SSD = 220 µm⁻²`): the initial Taylor
term is `1.0 × 25700 × 0.000286 × √220 = 109.0 MPa` instead of
`0.4 × ... = 43.6 MPa` — an initial slip resistance 2.50× too high, before any
calibration. The companion data sheet's warning about `τ_c` is therefore only
half the hazard: the additive decomposition *is* implemented correctly
(`τ_c0` is a separate lattice-friction term, not the total), but the Taylor
term itself is mis-scaled.

**Resolved** (`1c85ab2`). An `alpha0` parameter was added, defaulting to `1.0` —
the value the law was implicitly fixed at before the parameter existed — so
committed goldens and parity tests are unaffected. It is applied in
`critical_stress` **and** in `hardening_matrix` (`h^{αβ} = ∂τ_c^α/∂ρ^β`),
because a prefactor in the stress but not the tangent leaves the consistent
tangent inconsistent with what it differentiates. Study decks set
`alpha0 = 0.4`.

Verified against the hand derivation: `τ_c(0) = 176.7127 MPa` at the Al2219
set under the interaction-weighted reading settled in G13/C.15. Nine unit tests
in `tests/unit/materials/test_cp_taylor_hardening.py`; the decision and its
evidence are recorded in PINNACLE ADR 0022.

Fixing this also exposed a scalar-versus-batch parity defect en route:
`compute_batch` applied the interaction matrix to the SSD alone and added the
GND densities *outside* the contraction, while `compute` added them inside, so
the two paths disagreed whenever GNDs were present. Now identical.

### G3 — FIXED — Declared units of `k1`/`k2` contradict the implementation (`1c85ab2`)

`materials/dislocation_density.py:59-69` implements

```
dρ^α/dt = (1/b) (k1 √ρ^α − k2 ρ^α) |γ̇^α|
```

For `ρ̇` to have units `L⁻² T⁻¹` with `|γ̇|` in `T⁻¹` and `b` a length, the
bracket must carry `L⁻¹`. Since `√ρ ~ L⁻¹` and `ρ ~ L⁻²`, this requires **`k1`
dimensionless** and **`k2` with units of length**. The registry declares the
opposite: `k1` in `1/m`, `k2` dimensionless.

This matters because the papers report a dimensionless storage coefficient `k₀`
and an annihilation diameter `y_c` in length units — which map onto `k1` and
`k2` respectively under the implementation, and onto nothing coherent under the
declared units. The bundled deck's `k1 = 1.0e8`, `k2 = 10.0` are placeholders
that are non-physical under either reading.

The mapping itself (including whether the paper's law carries a factor 2 on the
annihilation term, i.e. `k2 = 2 y_c`) is derived and verified in Phase 1 against
the equation extracted from the papers; the saturation density is
`ρ* = (k1/k2)²`.

**Resolved** (`1c85ab2`). The units metadata now matches the implemented law —
`k1` dimensionless, `k2` a length — with the dimensional argument written into
the docstrings rather than only the conclusion, so a future reader can check it
instead of trusting it. An `annihilation_diameter` key accepts the papers'
parameterization directly (`k2 = 2 y_c`), and supplying both it and a
non-default `k2` is **refused** rather than resolved by precedence: a deck
setting both is asserting two different recovery rates, and silently choosing
one is how a parameter set gets misreported. `ρ* = (k1/k2)²` is documented.

The mapping was derived from the Al2219 form of the evolution law, which gives
`k1 = k_m` and `k2 = 2 y_c` — so `k1 = 0.3432`, `k2 = 0.00536 µm` for Al2219 and
`k1 = 0.0295`, `k2 = 1.58 × 10⁻⁶ mm` for Ni. Note the Ni paper prints the
annihilation term as `2/y_c`, which is dimensionally inconsistent; the Al2219
form is correct and was adopted for both materials. Full derivation in
`docs/03_material_point_verification.md`.

### G19 — FIXED — the reported stress field was not the equilibrium stress once plasticity was active

**Found while calibrating Gc, and it would have corrupted that calibration
silently.** `cp_quasi_static` builds its reported fields AFTER the load step
completes:

```
u_nodes   = u.reshape(mesh.n_nodes, dim)
strain_qp = assembler._compute_strain(u_nodes)
stress_qp = material.stress(strain_qp)      # <- after update_state
```

`material.stress` is a TRIAL stress: read-only with respect to the stored
constitutive state, which is exactly right inside a Newton iteration. But by
this point the state has already been advanced by `update_state`, so applying
the same TOTAL strain again charges the step's plastic increment a second time
and returns a stress that is not in equilibrium with the displacement field
beside it.

**Evidence, three independent checks.**

1. *Elastic regime — exact.* A homogeneous [001] single crystal at 0.02% strain
   gives sigma_zz = 12.73 MPa against the analytic E_001 * eps = 12.72 MPa
   (E_001 = (C11-C12)(C11+2C12)/(C11+C12) = 63,609 MPa), and sigma_xx = 0.000.
   At 0.1% strain: 63.75 against 63.61, sigma_xx = 0.001. With no plastic
   increment there is nothing to double-count, and the field is right.

2. *Plastic regime — a transverse stress that cannot exist.* The same
   homogeneous single crystal with free lateral faces must have sigma_xx = 0
   identically. It reports sigma_xx = 148-217 MPa against sigma_zz = 319-447
   MPa. The per-step ratios are 148.64/319.4 = 0.465, 217.13/330.3 = 0.657,
   152.39/447.0 = 0.341 and 211.42/418.4 = 0.505, i.e. 0.34-0.66. (An earlier
   revision wrote "0.5-0.8", which came from pairing the extremes of the two
   MPa ranges rather than dividing each step's pair.) For the polycrystal the
   volume average is 207 MPa
   against 695 MPa. The average-stress theorem makes the volume average of
   sigma_xx zero for a traction-free lateral surface, so this is not a
   heterogeneity effect.

3. *Non-monotonic axial stress.* Under monotonic loading with monotonic
   hardening, the single crystal reports 319 -> 330 -> 447 -> 418 MPa. A real
   equilibrium path cannot unload here.

**The kinematics are NOT affected**, which is what localises the defect to the
post-processing rather than the solve: the lateral faces contract by 4.40 um
mean against the incompressible expectation of 4.356 um at 20% axial strain, so
the Newton solve is finding the right displacement field. Only the stress
REPORTED beside it is wrong.

**Consequences.** Every stress-strain and load-displacement curve read from
this executioner in the plastic regime is wrong, which is most of what this
study reports. It would have entered the Gc calibration directly: Gc is fitted
so fracture initiates at the experimental strain, and the driving force is
built from the stress. Any parity golden carrying a plastic-regime stress field
encodes it too.

**Fix direction.** Capture the equilibrium stress DURING the solve -- the trial
stress at the converged displacement against the state as it stood at the START
of the step -- rather than re-evaluating after `update_state`. The same applies
to the per-step snapshot stress added for the damage-evidence work.

### G24 — FIXED — damage never fed back into the mechanics stress, on either route

**The single most consequential defect found in this study.** Softening is
the mechanism phase-field damage exists to represent: without it the driving
force never relaxes, damage spreads across the domain instead of localizing
into a band, and the stress-strain curve has no descending branch. Every
coupled run in this project before this fix was **one-way coupled** —
mechanics drove damage, damage never softened the material.

Three independent wiring faults, in series, each silent:

1. The coupled executioner gated the damage push on
   `hasattr(mech_exec, "kernel")`. Neither mechanics sub-executioner has that
   attribute — the kernels live on `problem.kernels` — so the guard was
   always false and **neither** route ever received a damage field. Same root
   cause as G22.
2. The live-callback accessor read only the CP kernel's `_material_cache`.
   The J2 kernel keeps its callback in `_callback`, so fixing (1) restored CP
   and left **J2 still undamaged**. Both names are now checked.
3. `CPFieldCallback` degrades only when the nodal damage can be scattered to
   quadrature points, which needs connectivity supplied by
   `set_gradient_operator` — wired only when GND was enabled. Every coupled
   CP deck with `gnd = false`, which is all of them here, scattered nothing.

**Measured, not inferred.**

| | before | after |
|---|---|---|
| CP: elements above d = 0.85 | 408 of 432 | — |
| CP: average axial stress there | **356 MPa, still rising** | degraded |
| CP: g(d) = (1-d)² + k demands | 0.0023 → about 0.8 MPa | — |
| CP: scatter calls / null returns | — | 802 / 0 |
| CP: min g(d) applied | never applied | **0.0014** |
| J2: min g(d) applied | never applied | **0.34** |

**This explains the diffuse damage seen throughout the study.** It was
repeatedly attributed to Gc being too low. With no softening the driving
force never relaxes, so damage grows everywhere rather than localizing —
the Gc calibration must be re-derived, not carried over.

**Regression protection.** `tests/unit/coupling/test_damage_feeds_back.py`
asserts the end-to-end property on both routes. It sets the load-step and
alternation counts itself, because damage is solved *after* mechanics within
an alternation, so a one-step one-alternation deck cannot show coupling
however correctly it is wired — and several bundled decks are configured
that way. It also asserts the chosen deck declares a degradation material,
so it cannot pass vacuously.

**An unexplained observation, recorded rather than smoothed over.** After
fix (1) alone, three J2 phase-field goldens moved
(`j2_pfd_function_hardening`, `..._finite_strain`, `j2_pfd_void_in_cube`).
After fixes (2) and (3) they pass again, and the suite is back to the
11-item baseline. By code reading, fix (1) should have changed nothing on
the J2 route, because the accessor still returned None for it until fix (2).
I have no verified explanation for the move-and-return and am not
substituting a guess. It should be settled before the J2 route is used for
any reported comparison.

### G27 — FIXED — the J2 route degraded compression along with tension

Closes the asymmetry left open by G26. The user supplied the reference,
MOOSE `ComputeCrackedStress`, which decomposes the uncracked stress
spectrally and degrades only its tensile part:

    sigma = (P+ h + P-) : sigma0,    P- = I4sym - P+

Because the projectors partition the identity this is exactly

    sigma = sigma0 - (1 - h) sigma+,   sigma+ = sum_a <lambda_a>+ n_a (x) n_a

PINNACLE multiplied the entire stress by g. Cracked material must sustain
compression but not tension; a uniform multiply removes the stiffness that
keeps crack faces apart and takes positive definiteness with it.

**Verification.** Matches the reference form `(P+ h + P-):sigma0` to
**7.1e-14** on a general non-diagonal stress. Limits pinned as tests:

| state | result |
|---|---|
| all principal stresses tensile | reduces to `g * sigma0` — the old behavior, exactly where it was right |
| all principal stresses compressive | **unchanged** — crack closure |
| mixed (+200, −150) | +200 → +0.46 degraded, **−150 preserved** |
| g = 1 | bit-exact no-op |

**Tangent.** The reference degrades the tangent with the same operator,
`J = (P+ h + P-) : J0`. A uniform g on C is no longer the derivative of the
split stress, so with damage bound the tangent is central-differenced from
the split stress. Undamaged J2 keeps its analytic tangent — the CP lesson
(G26) was that an inconsistent tangent is what breaks Newton.

**Modularity, as required.** With no damage bound the split path does not
execute, pinned by `test_degradation_leaves_an_undamaged_callback_untouched`
asserting bit-identical stress and tangent. Suite at the 11-item baseline,
3495 passing.

**Blast radius.** `j2_pfd_void_in_cube` regenerated: final stress +0.034%,
peak +0.73% — small because that case is largely tensile, with the
difference confined to the compressive region around the void.
`j2_pfd_function_hardening` and its finite-strain variant did not move at
all, consistent with being in tension throughout. The replaced unit test
asserted the uniform multiply; a tension-only test cannot distinguish the
two formulations, which is how this went unnoticed.

**Normalization difference, recorded not reconciled.** The reference writes
`h = (1-c)^2 (1-k) + k`, so `h(0) = 1` exactly; PINNACLE's `pfd_degradation`
supplies `g = (1-c)^2 + k`, so `g(0) = 1 + k`. They differ by O(k) = 1e-3
and agree at c = 1. The deck's g is used, so the **split** follows the
reference while the **degradation function** stays the deck's pluggable
choice. These are independent modeling decisions, which matters because the
user plans to add further splits — of the four in common use, two are now
implemented: the volumetric split on S for CP (G26) and the spectral
stress split for J2 (this entry).

### G26 — FIXED — degradation belonged on the 2nd Piola-Kirchhoff stress, not the Cauchy stress

Raised by the user, who asked why the source formulation's degradation of
**S** had not been implemented. It had not, and the substitute was wrong in
three ways that together explain G25 and the diffuse damage.

**The source formulation**, Eq. (29), puts the degradation inside the
constitutive law with a volumetric split:

| branch | S |
|---|---|
| Jₑ ≥ 1 | `[(1-c)² + k] C:Eᵉ` |
| Jₑ < 1 | `[(1-c)² + k] (C:Eᵉ − Jₑ^{2/3} K ω Cₑ⁻¹) + Jₑ^{2/3} K ω Cₑ⁻¹` |

with ω = (3/2)(Jₑ^{2/3} − 1) and K = (1/9)C_iijj. PINNACLE instead multiplied
the **Cauchy** stress by g(d) after the material returned.

For Jₑ ≥ 1 the two are algebraically identical, since σ = (1/J)FᵉS(Fᵉ)ᵀ is
linear in S and g is scalar. Everywhere else the post-hoc form fails:

1. **Compression.** The split shields the volumetric part so a broken point
   still resists interpenetration. A scalar multiply on σ degrades
   compression too, removing the stiffness that keeps crack faces apart,
   and positive definiteness with it.
2. **Plastic driving force.** τ = M : S. With S undegraded, slip runs at
   full rate in fully damaged material — which is what exhausted the CP
   sub-step budget in the damaged region and abandoned the load step (G25).
3. **Tangent.** Both tangent paths finite-difference the material directly
   and never saw g(d), so Newton linearized a stress up to 1/g ≈ 400x
   stiffer than the residual it was correcting.

**Implementation.** `degrade_pk2` realizes Eq. (29) branch-free through the
negative part of ω: `S = g(S − Sv) + Sv` with `Sv = Jₑ^{2/3} K ω⁻ Cₑ⁻¹` and
`ω⁻ = min(ω, 0)`. Identical to the two branches, and continuous at Jₑ = 1,
where a hard branch would put a step in the tangent at the state an element
passes through as it unloads. Verified against the paper's branches to
**0.0e+00** in both tension and compression; batched equals scalar; g = 1 is
bit-exact. Degradation is applied in `stress`, both tangent paths, and the
adaptive sub-stepping loop — everywhere τ is formed.

**Modularity, as required.** `degradation_g` defaults to `None` and the path
is then a strict no-op, so CP-only and J2-only runs are untouched. Confirmed
by the full suite: 11 failures, the exact pre-existing baseline, and 3494
passing, up 3 from the new degradation tests. No CP-only or J2-only golden
moved.

**Effect — the physics now behaves.** Perforated plate, Gc = 2.0 MPa·µm,
12 load steps:

| strain | damaged fraction | peak d | avg σ_xx |
|---|---|---|---|
| 0.0025 | 0.000 | 0.459 | 153.9 |
| 0.0050 | 0.003 | 0.865 | **192.0** (peak) |
| 0.0075 | 0.145 | 0.991 | 100.1 |
| 0.0100 | 0.206 | 0.999 | 21.2 |
| 0.0125 | 0.206 | 1.000 | **0.93** (severed) |

A complete fracture curve with a descending branch, the first this coupling
has produced, and damage **localized to 20.6%** of the domain against the
88% diffuse smear before. The diffuse damage attributed for weeks to Gc
being too low was this defect.

**The J2 route still uses the post-hoc Cauchy multiply.** Deliberate and
recorded rather than overlooked: Eq. (29) is the CP formulation, and
imposing a split on J2 without a reference for it would be guessing.

### G25 — OPEN, severity A — a fully softened solve returns zero boundary displacement

Surfaced immediately once softening worked. On a perforated plate at
Gc = 2.0 MPa·µm — far below the calibrated range, so the specimen fails
almost at once — the first load step is normal (boundary strain 0.0100,
average axial stress 294 MPa), and from the second step onward the reported
boundary strain is **0.00000** with an average axial stress of −15 MPa,
while damage freezes at 0.414 of the domain.

The boundary displacement is a **prescribed Dirichlet value**. It must be
honored whatever the material stiffness; a softened structure should show
falling stress at the imposed displacement, not zero displacement. So this
is not "the specimen broke" — it is the solve returning a state that does
not satisfy its own boundary conditions, and continuing without raising.

**Mechanism established.** It is not the linear solver: CG converged on all
2,946 calls (`info = 0`), and the boundary displacement is *exactly* zero at
every node, which is a boundary condition never applied rather than
numerical drift. The CP integrator exhausts its sub-step budget in the
damaged region and `cp_quasi_static` reaches the cutback floor, where it
appends `{"converged": False, "error": "newton_cutback_floor: ..."}` to the
mechanics history, advances `step_idx`, and **continues** — the load
increment is never solved, and the run returns success at the top level.
A silent failure recorded in a structure nothing reads, which is the same
family as G13, G14, G19, G22 and G24.

**Guarded on the study side, not yet fixed in the framework.**
`run_plate.py` now verifies the *physical invariant* after every run: the
measured boundary strain must equal `target_strain x load_factor` at every
step. On violation it writes `INVALID.txt` and refuses to emit
`results.csv`, since a results file is consumed downstream as a solved
result. Guarding on the reported `converged` flag does **not** work and was
tried first: with a single alternation the staggered record is
`converged: False` by construction, so the flag cannot distinguish a lagged
coupling from an abandoned load step. Validated against three saved runs —
it rejects the failed one (step 1: prescribed strain 0.02000, actual
0.00000) and accepts both healthy ones.

The framework-side question — whether the cutback floor should raise, or
sub-cycle the damage field, or accept a reduced load factor and report it —
is a design decision that has not been made. It must close before any run is taken
to full fracture, which is every run the fracture study needs. Note it
would have been invisible without the new `results.csv`: the boundary
strain column is what exposed it.

### G22 — FIXED — the material-state snapshot and restore never did anything

**Found while implementing restart, and it invalidated work I had already
reported as done.** `_mechanics_material` read `sub_exec.kernel`.
`CPQuasiStaticExecutioner` has no such attribute — the kernels live on
`problem.kernels` — so the accessor returned `None` on **every** call, and
every snapshot and restore built on it was a no-op. That includes the
alternation state restore introduced with ADR 0024, whose whole purpose is to
stop each alternation piling plasticity on the previous one.

A bare `except Exception: return None` wrapped the lookup. That is what let a
no-op pass for working code for as long as it did, and it is gone: a kernel
that advertises `make_material` and then fails is a real defect and must
surface.

**Two further defects surfaced the moment the exception stopped being
swallowed**, neither of which could have been seen while the accessor
returned None:

* `make_material(backend, dim)` cannot supply `n_qp`, which is part of the
  material cache key, so the call always MISSED and tried to construct a
  second callback. On an orientation-field deck it fails outright with an
  orientation-shape mismatch; on a deck where it succeeds it would have
  returned a *fresh* object carrying virgin state, so the snapshot would have
  operated on something the sub-solve never uses. The live cached instance is
  read instead.
* The phase-field irreversibility history `H = max(H, psi+)` lives on the
  driving-force material, not in the damage field, so a restart reinstating
  only `d` would resume with `H = 0`.

**Blast radius on committed results is zero, and that is a measurement, not a
hope.** The full suite was re-run after the fix and **no golden moved** (12
failures against the 11-item baseline, the twelfth being stale generated
syntax docs from three new parameters). The bundled coupled cases converge in
a single alternation, and with one alternation there is nothing to restore.

**Where it does bite: this study's own runs.** The Gc calibration sweep was
dispatched with three alternations per load step. Wherever a step took more
than one, plasticity was counted once per alternation. All six runs were
cancelled rather than harvested.

### G23 — OPEN, severity B — coupled restart is implemented but not bit-exact

`checkpoint_every`, `checkpoint_path` and `restart_from` are declared on the
staggered executioner and write the quadrature-point constitutive state, the
damage field, the irreversibility history, the displacement and the load-step
index. Writing a checkpoint is verified to be bitwise non-perturbing: a run
with checkpointing enabled produces results identical to one without.

**Resuming is not yet equivalent.** A run restarted from step 3 of 6 differs
from an uninterrupted run by 9.5e-4 relative in displacement and 7.1% in
stress. The cause is not yet identified — the constitutive state, damage
field and history are all confirmed present in the checkpoint. Until this
closes the parameters exist but must not be used for production runs, and the
campaign remains limited to what fits in one allocation.

The equivalence test itself is the right acceptance criterion and should be
kept: a restart that merely looks plausible is worse than no restart, because
a long campaign would silently accumulate the discrepancy at every resume.

### G21 — FIXED — the direct linear solve, not the constitutive integration, is what made the fine meshes unreachable

**Found by testing an assumption in my own cost model.** I had extrapolated
campaign cost linearly in element count. A sparse direct solve does not scale
linearly, so I measured it on two levels of the void family, one load step
each:

| level | elements | total | linear solve | tangent |
|---|---|---|---|---|
| d=16 | 5,176 | 156.9 s | 22.2 s (14%) | 44.6 s |
| d=24 | 16,768 | 836.7 s | **677.4 s (81%)** | 324.5 s |

Against a 3.24x increase in elements, the direct solve grew **30.6x** — a
scaling exponent of **2.91**, near-cubic — while the run total grew 5.33x
(exponent 1.42). Extrapolated, the direct solve alone reaches roughly 10,000 s
per Newton iteration at d=32 and is hopeless beyond it. The near-cubic term is
what makes the fine meshes unreachable; the constitutive integration, which I
had assumed was the constraint, is not.

**A Jacobi-preconditioned CG path already existed and no deck used it.**
`linear.py` has offered `method="cg"` alongside `"direct"` all along; every
study deck and every bundled example specifies `direct`.

| level | direct | cg | speedup | max relative displacement difference |
|---|---|---|---|---|
| d=16 | 156.9 s | 157.4 s | 1.00x | 8.3e-15 |
| d=24 | 820.2 s | 488.1 s | **1.68x** | 2.7e-14 |

CG returns the *same answer* — the displacement fields agree to 3e-14 and the
volume-averaged axial stress to six significant figures — so this is a pure
cost change with no accuracy cost. Per-element cost is flat under CG (30.3
against 29.1 ms/element at the two levels), i.e. the total scales roughly
linearly, where the direct solve made it superlinear.

**Action taken.** Every study deck now specifies `linear_method = cg`. The
framework default is left alone: at small sizes the two are indistinguishable,
and changing a default would move bundled goldens for no benefit.

**Caveat worth carrying.** CG assumes a symmetric positive-definite system.
The agreement above shows the CP tangent is close enough to symmetric for
Jacobi-CG to converge to the direct answer at these two levels, but that is a
measurement at d=16 and d=24 under a specific load state, not a proof. The
convergence study should re-verify CG against direct at one damaged, strongly
non-symmetric state before the fine-mesh results are trusted; `linear.py`'s own
TODO notes that AMG/ILU preconditioners remain future work if Jacobi stalls.

### G20 — FIXED — load-path resolution is a first-order discretization axis, and the default was far too coarse

**Found by validating the G19 fix, and it is larger than any mesh effect
measured so far.** With the corrected stress in place, a homogeneous [001]
single crystal was compared against a material-point integration of the same
constitutive law — an independent reference that shares no code with the field
solve. They must agree, and whether they do is controlled by how finely the
load path is subdivided:

| mechanics sub-solve `n_steps` | wall (s) | sigma_zz at 2% strain |
|---|---|---|
| 1 | 23.4 | 866.8 |
| 5 | 92.5 | 269.6 |
| 25 | 355.2 | **309.9** |

The material-point reference is **309.9 MPa**. At 25 sub-steps the field solve
reproduces it exactly; at 1 sub-step — the setting every deck in this study
currently uses — it reads 866.8 MPa, a **180% error**. The material point tells
the same story on its own: integrating to 2% strain in 2 increments returns
-39,879 MPa (a negative stress), in 4 increments 475.9, and converges to 309.9
only by about 25 increments, i.e. a strain increment near 8e-4.

**This is not the mesh.** Both the 1-sub-step and 25-sub-step runs use the same
mesh, the same material, and the same total time, so the same mean strain rate;
only the subdivision of the load path differs. The rate-sensitive flow rule
(m = 0.12) makes the constitutive integration error large and one-sided.

**Consequences for the study.**

* Every coupled result produced so far is time-under-resolved, and the error is
  larger than the mesh effects the convergence study was designed to measure. A
  mesh-convergence study run at 1 sub-step would have converged the wrong
  answer.
* The cost model changes, though by less than the sub-step count suggests --
  see the measured 4.7x below, not 8x, because the state update does not
  scale with subdivision.
  Subdividing the MECHANICS sub-solve rather than the staggered loop is the
  cheaper lever, since the damage solve still runs once per staggered step.
* Load-path resolution must be a declared axis of the convergence study
  alongside element size and l0, and must be converged BEFORE any mesh
  conclusion is drawn.

**Mechanism, established rather than inferred.** `update_state` sub-steps
adaptively (shear-increment cap 1e-3, budget 1e5), so the internal STATE
integrates accurately however coarse the load step. `stress()` does not: it
calls `material._evaluate(F, state)` at the FROZEN state with no `dt` and no
sub-stepping, and that trial stress is what drives the Newton equilibrium and
what is reported. Over a large increment the frozen state is stale, so
equilibrium is solved against a stress that is too stiff. This is why
sub-stepping the STATE does not rescue the answer and only a smaller LOAD
increment does.

**Threshold, measured against a 2000-increment reference.**

| target strain | d_eps = 8e-4 | d_eps = 2.4e-3 | d_eps >= 4.8e-3 |
|---|---|---|---|
| 2% | 0.04% | — | — |
| 12% | 0.05% | 11.1% | 3.3e10 % (diverged) |
| 20% | 0.06% | 1.7% (at 2.0e-3) | 7.1e7 % (diverged) |

The requirement is **d_eps <= 8e-4 per constitutive update**, independent of
target strain. The failure is two-sided: above ~2.4e-3 the integration
diverges visibly (1e11 to 1e14 MPa, or a negative stress), but between the
threshold and that point it returns a plausible number that is simply wrong.

**Resolution.** Three parts, because no single one is sufficient:

1. *Requirement.* d_eps <= 8e-4, recorded here with the measurement behind it.
2. *Enforcement.* `run_polycrystal.py` computes the mechanics sub-step count
   from the target strain and load-step count, and records the achieved
   increment in every result record. Subdividing the MECHANICS sub-solve
   rather than the staggered loop keeps the damage solve at one per load step.
3. *Guard.* `CPFieldCallback._check_strain_increment` warns once when the
   increment exceeds the threshold, on both the scalar and batched paths. It
   fires on 14 bundled decks, which is the point -- they are under-resolved
   today and nothing said so.

**Cost — measured, and lower than the naive estimate.** Subdividing 8x costs
**4.7x**, not 8x: 188.3 s at one mechanics sub-step against 888.9 s at eight
(void d=16). A profile explains why. The constitutive state update is 60% of
the runtime and its own adaptive loop already takes the sub-steps the strain
demands -- 342 `compute_batch` calls for a single load step -- so that 60% is
set by the total strain and the shear cap, not by how the load step is
subdivided. Only the Newton work (tangent, residual, linear solve) multiplies.

An earlier revision of this entry quoted ~15x by assuming the whole run scaled
with the sub-step count. It does not, and the measurement above supersedes it.

**Still to establish.** Whether the threshold tightens under damage coupling,
where the degraded stiffness changes the local strain rate. The Gc sweep
records the achieved increment per run so this can be checked against its
results rather than assumed.

### G17 / G18 — FIXED — The coupled solve was unaffordable for reasons unrelated to the physics (`7156159`)

Found by profiling the first production run, after the coarsest level alone ran
41 minutes and a two-load-step probe exceeded 900 s. Both defects are cost, not
correctness — and both were invisible until a real mesh was put through.

**G17 — the Newton solver hardcoded a dense tangent.** `nonlinear.py` asked for
`format="dense"`, and the assembler's COO path, whose own module docstring says
it is "for larger field solves", was reachable from no solver at all
(`hybrid.py` hardcoded dense too). A dense tangent costs `n_dofs² × 8` bytes:

| level | elements | nodes | DOF | dense tangent | allocatable in 60 GB |
|---|---|---|---|---|---|
| d = 16 | 5,176 | 6,096 | 18,288 | 2.7 GB | yes |
| d = 24 | 16,768 | 18,756 | 56,268 | 25.3 GB | yes |
| d = 32 | 43,136 | 46,730 | 140,190 | 157 GB | **no** |
| d = 48 | 137,216 | 145,034 | 435,102 | 1,515 GB | **no** |
| d = 64 | 338,616 | 352,544 | 1,057,632 | 8,949 GB | **no** |

Three of the five refinement levels could not allocate the matrix at all, so the
mesh-convergence study was impossible before it began.
`BoundaryConditions.apply_system_sparse` mirrors the dense elimination exactly,
including the prescribed-increment right-hand-side correction — which is what
makes the two agree rather than merely both converge. Measured on an identical
coupled problem: `max|u_dense − u_sparse| = 1.8e-9` (3.1e-10 relative) and
`max|d_dense − d_sparse| = 1.1e-12`.

The switch-over threshold was **measured, not guessed**: on a banded stand-in the
dense solve wins at 5,000 DOF (0.37 s against 2.28) and at 10,000 (2.46 against
4.91), and loses by 3.8× at 18,288 (35.6 against 9.3). The factorization cost
overtakes long before memory does, so the crossover sits near 12,000 DOF — far
below where the memory argument alone would have put it, which is why the first
threshold I set from memory reasoning (20,000) was wrong and left the coarsest
level on the slow path.

**G18 — the CP tangent finite-differenced one quadrature point at a time.** The
`tangent_method = analytic` parameter is a misnomer: the tangent is a
finite-difference stencil, and it was evaluated in a triple Python loop over
quadrature points × 3 × 3 perturbation components. Profiling one load step at
the coarsest mesh:

- **716.7 s** total for a single load step with a single alternation;
- **449.5 s** in six tangent assemblies, 74.9 s each, 63 % of the run;
- **2,484,480** calls to `_stress_numpy` for 41,408 quadrature points — exactly
  ten per point per tangent (one base plus nine perturbations).

The stencil is unchanged; the points are now evaluated together through
`compute_batch(..., dt=None)`, which returns the response without advancing the
state and whose returned state is discarded, since a tangent must leave the
material as it found it. The batched route declines and falls back to the
per-QP loop under the same conditions `update_state` does — no `compute_batch`,
a per-QP orientation field, or an active grain-boundary multiplier — so the two
cannot diverge in which path they take.

Measured: the tangent alone is **12.3× faster** at 512 quadrature points and
agrees with the loop it replaces to 1.1e-9 relative; the full profiled case goes
**716.7 s → 185.2 s, 3.9×**.

**No parity golden moved for either change**, which is the point: both are
numerically equivalent to what they replace, unlike G15 and G16 which changed
the answer because the answer had been wrong.

**Cost is still the open constraint.** At 185 s for one load step at the
*coarsest* level, a 20-step run is about an hour at d = 16, and d = 64 carries
65× the elements. Sizing the production campaign against measured cost — rather
than assuming the planned matrix is affordable — is the next step, and it is now
a scheduling question rather than a feasibility one.

### G16 — FIXED — `staggered` had no alternation loop, so the coupling was lagged (`4340932`, ADR 0024)

Exposed by fixing G15: it was the remaining blocker on a *converged* coupled
reference. The executioner carried `max_iterations` (default 50) and
`tolerance`, and drives neither: the outer loop runs mechanics once and damage
once per load step, and `"iterations": 1` is **hardcoded** into the history it
reports. So at each load level the mechanics is solved against the *previous*
level's damage field — an explicit, lagged coupling — rather than alternated to
convergence.

**There is a semantic conflict underneath it, and it changes the meaning of
existing decks.** `staggered`'s `n_steps` has been doing double duty. The
`step_factor = (step + 1) / n_steps` expression has always implied load steps,
which is how G15's fix now uses it. But the bundled decks were written to the
de-facto behavior, and say so — `void_in_cube_j2_pfd.i` comments that
"n_steps = 1 would leave the mechanics undegraded and the coupling one-way",
which only makes sense if `n_steps` is an alternation count at fixed full load.
Both readings are in the tree, and the parameters that would separate them
(`max_iterations` for alternations, `n_steps` for load steps) exist with only
one of them wired.

**The semantics decision this forced, recorded as ADR 0024.** `n_steps` is now a
load loop, and `max_iterations`/`tolerance` drive the alternation. That reading
wins because it can express the other — a single load step with several
alternations is `n_steps = 1, max_iterations = k` — while the alternation
reading cannot express a load path at all, and under it the study's primary
metrics do not exist in the output. Introducing a new parameter to preserve the
bundled decks' reading was considered and rejected: G15 had already moved their
answers, so the protection would have been illusory while leaving `n_steps`
permanently ambiguous. The declared default of 50 alternations is kept, because
the cost is bounded by `tolerance` — convergence took 2 — and lowering it to 1
would ship the lagged coupling as the default and make correctness opt-in.

**A further requirement, which is why this was not folded into G15.** Alternating
at a *fixed* load level requires the constitutive state to be restored to the
state converged at the end of the previous load level before each alternation —
otherwise every alternation re-integrates plasticity on top of the last one and
the plastic strain is counted several times per load step. So G16 needs state
*snapshot and restore* at load-step boundaries, not merely the continuation G15
added. `CPFieldCallback` already exposes `load_state_dict`, so the mechanism
exists, and `stress_divergence_j2` already cached its callback so that "the
plastic state survives across load steps and coupling iterations" — the CP cache
added under G15 brought that route into line with a convention already in the
tree. Snapshot and restore duck-type on those two methods and are inert when the
mechanics kernel has neither, since an elastic problem has no path-dependent
state. Note the error this prevents grows *with* `max_iterations`, so without it
a tighter convergence setting would have produced a worse answer.

**Seven parity goldens were regenerated, and they are the coherent blast
radius.** `elastic_pfd_fem_bar_1d`, `elastic_pfd_hybrid_bar_1d`,
`j2_pfd_function_hardening`, `j2_pfd_function_hardening_finite_strain`,
`j2_pfd_void_in_cube`, `rezwan_field_fem_3d`, and `rezwan_field_hybrid_3d` — every
one a coupled `staggered` or hybrid case, i.e. exactly the set that alternates.
The `rezwan_field_*` diffs show the mechanism directly, moving on
`solver_history` columns `iterations` 1 → 2 and `converged` 0 → 1.

**The converged solutions carry lower stress, and that direction is the
physically expected one.** `j2_pfd_void_in_cube` final stress fell from a golden
80.90 MPa to 74.80 MPa, about 7.5 %. A lagged coupling solves the mechanics
against a damage field that has not yet responded to the current load, so it
under-damages and over-predicts load capacity. Any coupled result produced
before this change should be read as an upper bound on stress.

**On the golden that was left failing under G15.** `j2_pfd_void_in_cube` was
deliberately not regenerated then, so a lagged coupling would not be frozen as
the reference; that objection is resolved here and it is regenerated with the
rest. The open question recorded at the time — why that case moved under G15
while the neighboring `j2_pfd_function_hardening` did not, when both drive
mechanics with `fem_quasi_static`, both set `ramp_dirichlet = true`, both run
`staggered` with `n_steps = 3`, and both have their mechanics `n_steps`
overridden to 2 by the harness — is now partly answered: `j2_pfd_function_hardening`
*did* move once the alternation loop landed. So the two cases differ in
sensitivity to load windowing specifically, not in whether they respond to the
coupled-solver changes at all. The root cause of that differing sensitivity is
still not established, and is not pursued further because both goldens now
encode the converged scheme.

**Runtime rises** by roughly the mean alternation count per load step — about 2×
on the cases measured. Convergence took 2 alternations, not the declared cap of
50, so the cost is bounded by `tolerance` rather than by `max_iterations`.

The failure set is back to the 11 pre-existing: 10 parity failures plus
`test_verify_install_passes`, which also fails on the unmodified tree here.
3491 pass.

**Consequence for the study.** G15 made the solution path-correct in load and in
constitutive state, which is what made the metrics *defined*; G16 makes the
coupled solution at each load level *converged*. With both closed, the coupled
reference is fit to benchmark against — which unblocks the discretization and
solver-convergence step, and Phase 2 behind it.

### G15 — FIXED — `staggered` never incremented the load, and reset the CP state every sub-solve (`43dd6fd`, ADR 0023)

**This was the gating blocker for Phase 2 and the most consequential finding of
the audit.** Before the fix the coupled executioner could not produce a
load-stepped solution at all, so the study's primary metric —
fracture-initiation strain, reported as 7.0 % for isochoric [001] in
(Rezwan 2025 intragranular) — was not measurable, and neither was a
load–displacement curve.

Two defects compound:

**(a) The load is never incremented across outer steps.** The outer loop
computes `step_factor = (step + 1) / n_steps`, but `_run_mechanics` is called
*without it*, and `_compute_psi_plus`'s own docstring says `step_factor` is
"unused here and kept for subclasses that ramp a prescribed field instead". So
the value is computed and used nowhere. Every outer step re-runs the mechanics
sub-solve over its **full** load program.

**(b) The CP state is reset to virgin on every sub-solve.**
`cp_quasi_static.run` calls `kernel.make_material(backend, dim)` each time, and
`StressDivergenceCPKernel.make_material` returns a **fresh** `CPFieldCallback`
whose constructor does `self.state = self._initial_state()`. So the plastic
deformation gradient, the dislocation densities, and the accumulated plastic
work are all discarded between outer steps.

Together these make each outer step an independent, full-load, virgin-material
solve, with only the damage history field carrying over.

**Measured, not only read.** Running the production physics on a 2×2×2 box at
outer `n_steps` = 1 and 3:

| outer `n_steps` | axial stretch | peak damage |
|---|---|---|
| 1 | 1.12000000 | 1.0 |
| 3 | 1.12000000 | 1.0 |

`max|u(n=1) − u(n=3)| = 0.000000e+00` — the displacement field is
**bit-identical**. The first outer step is already at the full 12 % strain, and
the additional outer steps change the mechanics not at all. Damage saturates at
1.0 in both, which is what a single full-load application to virgin material
produces.

**Consequences to carry forward.** Any coupled result from this executioner is
a *fully-loaded end state*, not a path. That means: no initiation strain, no
load–displacement response, and — because of (b) — a plastic-work driving-force
contribution (G1a) computed from a single virgin full-load solve rather than
accumulated along the path, which is not the quantity Eq. (31) defines.

**How it was closed.** Three changes in dependency order, with ADR 0023 recording
the ownership decision and the alternative rejected:
1. **`cp_quasi_static` honors its own declared load window.** It declared
   `initial_load_factor` and `target_load_factor` and wired *neither* into its
   `LoadProgram`, though `fem_quasi_static` honors both — so the continuation
   the fix needed did not exist yet. This is the same failure mode as G13 and
   G14: a capability declared in the schema, implemented in a sibling, and
   unwired here. The initial factor now reaches `LoadProgram`; the target
   rescales the span, applied as a wrapper rather than by mutating the program,
   since `AdaptiveLoadStepper` holds a reference and its cutback logic depends
   on the object it was handed.
2. **The CP callback is cached on the kernel.** Keyed on the backend object
   (held by reference, so the `id()` cannot alias a later object that reused the
   address), its `device` and `precision`, `dim`, `n_qp`, and the config flags —
   everything that would invalidate the stored tensors. The device and precision
   keys are load-bearing: `backend/base.py` documents that a tensor built before
   a later `set_device` stays on the old device, so a stale hit would surface as
   a device mismatch deep inside an operator.
3. **The staggered loop sets the window** `[step/n_steps, (step+1)/n_steps]` on
   the mechanics sub-solve, with continuation enabled after the first step.
   `continue_state` **defaults to false**, so every single-solve deck behaves
   exactly as before, and the failure mode of forgetting the flag is a
   *restarted* path — visibly wrong — rather than a *contaminated* one.

**Verification.** With `Gc` raised so damage does not saturate immediately,
per-step damage increment over four steps at load 0.25 → 1.00:

| | per-step damage increment |
|---|---|
| before | `0.1506, 0.0, 0.0, 0.0` |
| continuation on (shipped) | `0.00993, 0.01361, 0.01401, 0.01402` |
| continuation off | `0.00993, 0.02863, 0.04705, 0.06503` |

Before, damage appeared entirely at the first step and never changed as the load
rose — the signature of one full-load solve repeated. With the fix the increments
**converge**, as a path-correct solution should: each increment adds
progressively less new damage as the field equilibrates. With continuation
deliberately disabled they **grow**, because each step re-plasticizes virgin
material at a higher load. The two sequences differ qualitatively, not just
numerically, which is what makes this a correctness fix rather than a refinement.

At the paper's own `Gc` = 22.75 MPa·µm damage still saturates to 1.0 at the first
load step on a homogeneous cube at 3 % strain. That is not a plumbing problem but
direct evidence for gap **C.16**: the reported value cannot be used as a
transferable material property and must be treated as an effective regularized
quantity, to be settled in calibration.

`tests/unit/coupling/test_staggered_load_path.py` (8 tests) covers the window
parameters, the default-off continuation flag, `LoadProgram`'s blend, callback
reuse, state survival across calls, reset-to-virgin, and `n_qp` invalidating the
cache rather than returning a wrong-shaped state.

**What remains: G16.** This fix makes the solution path-correct in load and in
constitutive state, which is what made the study's metrics *defined*. It does
not make the coupled solution at each load level *converged* — that is G16, the
missing alternation loop, recorded above.

No accuracy or cost comparison between solver routes should be reported until
**G16** also closes: a benchmark reference needs both a correct path and a
converged solution at each point on it.

### G14 — FIXED — The Exodus reader silently discarded every file-provided node set (`4e48301`)

Found while giving the void meshes explicit named sets to avoid the hazard in
G15's neighborhood. The reader looked for node sets on `mio_mesh.nodedata`; a
meshio `Mesh` has no such attribute (it exposes `point_data`, `cell_data`,
`point_sets`, `cell_sets`), so the `hasattr` guard never fired and every
file-provided node set was dropped. Only the sets `_derive_boundary_sets` infers
survived — which is why its documented promise that "sets that the file *did*
provide are never overwritten" had never taken effect: none ever arrived for it
to preserve. It uses `setdefault` and was keeping its side of the contract.

**Why it matters physically, and it is the same trap as G15.** The derived
`all_boundary` is *topological*, so on a mesh with an internal void it includes
the void surface, and `pfd_steady` defaults `boundary_node_set` to it. A deck on
the default therefore pins `d = 0` on the void surface — the one place damage
must nucleate — and the run completes and looks converged. A mesher that names
its own outer-boundary set precisely to avoid this had the set silently
discarded.

**Verification carries its own cross-check.** With the fix, the study meshes'
`void_surface` and `outer_boundary` sets import, and at every level their sizes
sum **exactly** to the independently derived topological `all_boundary` —
218 + 1538 = 1756 at d = 16, and 2906 + 24578 = 27484 at d = 64 — with zero
overlap, and every `void_surface` node sits at radius 5.000000 µm. Two
independent constructions agreeing exactly is stronger evidence than either
alone, and it also confirms the indices are zero-based as assumed.

`tests/unit/mesh/test_exodus_node_sets.py` (5 tests) covers import, coexistence
with the derived sets, a file set overriding a colliding derived name,
`name_mapping` applied to file sets, and zero-based indexing.

### G5 — FIXED — No conforming spherical-void mesh with local grading

Met entirely from hexaform plus PINNACLE's existing `[Mesh] type = file` import;
**no framework change was needed**. Generator:
`scripts/make_void_meshes.py` in the study repository; family and per-level
verification record in `meshes/al2219/mesh_family.json`.

**Route.** Hexaform meshes the sphere as an *inclusion* on a cubed-sphere O-grid
that places nodes on the analytic surface; dropping that block leaves a cavity
whose boundary is the conforming O-grid surface. The void is the absence of the
inclusion, so it inherits the inclusion mesh's geometric convergence. This is
why the study does not use the bundled `box_with_void`, which deletes background
elements by centroid test and therefore staircases the surface — precisely where
damage nucleates.

**The family** (50 µm cube, 5 µm centered void, `l0` = 1 µm, refined radius
8 µm = r + 3·l0):

| d | elements | background h (µm) | first layer at void (µm) | l0 / dr₁ | min scaled J | inclusion vol. error |
|---|---|---|---|---|---|---|
| 16 | 5,176 | 3.125 | 0.3797 | 2.63 | 0.5774 | 3.01e−2 |
| 24 | 16,768 | 2.083 | 0.2167 | 4.61 | 0.5774 | 1.71e−2 |
| 32 | 43,136 | 1.562 | 0.1805 | 5.54 | 0.5774 | 7.64e−3 |
| 48 | 137,216 | 1.042 | 0.0993 | 10.07 | 0.5774 | 4.31e−3 |
| 64 | 338,616 | 0.781 | 0.0718 | 13.92 | 0.5774 | 2.28e−3 |

Zero inverted or negative-Jacobian elements at every level; minimum scaled
Jacobian is 0.5774 = 1/√3 throughout, the cubed-sphere corner value, so it is
constant by construction rather than by luck. Exodus round-trip is exact
(maximum coordinate error 0.0) and surface nodes sit on the analytic sphere to
2.7e−15 µm. Inclusion-volume error converges at mean observed order **1.95**
across the family, consistent with second order; the pair-to-pair oscillation
(1.40, 2.80, 1.41, 2.21) is because the O-grid block size advances in even
integers so it does not scale exactly with density.

Verified through the solver, not only in the mesher: PINNACLE imports each level
and **derives** the boundary sets, giving all six bounding-box node sets plus a
*topological* `all_boundary`. At d = 16 that resolves 218 void-surface nodes, all
at radius exactly 5.000000 µm — the surface is conforming, not staircased — which
is what makes the void a free surface the isochoric BC set can be applied around.

**A measurement changed the design, and this is the substantive point.** A first
version set only the refined radius and left `MacroBlockOptions.n_outer` at its
default of 2. Two radial layers then spanned r = 5 → 8 µm *regardless of
density*, so refinement went into the background and not into the region that
matters: across the family the background element size fell 3.125 → 0.781 µm
while the first layer at the void surface moved only 3.15 → 1.93 µm, leaving
**0.3–0.5 elements across l0 at every level**. A convergence study on that family
would have converged the background and reported the result as mesh convergence,
with the damage band unresolved throughout. The controls existed — `n_outer` and
`outer_grading` — and the generator simply had not used them. (Hexaform
documents the grading ratio on the companion key `inner_grading`, ">1 clusters
at the interface (where plastic strain and phase-field gradients concentrate)";
`outer_grading` is documented only as "Same for the matrix-side shell", which is
the side the material surrounding the void lies on. The achieved first-layer
thickness is therefore measured here rather than inferred from that
documentation.) The generator now sets both, measures
the achieved first-layer thickness per level rather than assuming it, and
refuses a family whose **finest** level fails the resolution target.

Two design points that follow, and both are deliberate:

- `n_outer` **scales with the level** (6, 9, 14, 20, 30). Holding it fixed makes
  dr₁ nearly constant (measured 0.164–0.210 µm at `n_outer` = 12), so the local
  fields would be computed on the same radial discretization at every level. A
  refinement family has to refine everywhere; the script warns when dr₁ is not
  strictly decreasing.
- The resolution gate applies to the **finest** level only. The coarsest level
  (d = 16, 2.63 elements across l0) is deliberately under-resolved — a
  convergence study needs levels that are demonstrably too coarse, or
  convergence is not observable — and is flagged `under_resolved_at_l0` in the
  record so its local fields are never quoted without that qualification.

**Deferred audit questions, now answered.** The sphere can be emitted as a
removed void (by dropping the inclusion block, above), and the refined zone *is*
independently controllable: the macroblock extent sets the refined radius and
`n_outer`/`outer_grading` set the radial resolution within it. Neither was a
hexaform gap.

### G4 — FIXED — No isochoric (constant-volume) loading driver (`84bebc7`)

**No new component was needed, and that is the result.** The study plan requires
composition from registered components to be attempted before any framework
work, and here it succeeds: the paper's constant-volume loading is a `parsed`
function on an ordinary `dirichlet` BC.

The paper specifies (Sec. 2.4 / 3.1.1, independently extracted): *"tension
applied along the crystallographic direction on one face; compression on the two
reciprocal directions; three faces loaded; load ratios adjusted dynamically so
that det(F) = 1 globally; remaining faces fixed."*

With three faces held on their own normals — a roller, removing rigid-body
motion — and the other three given prescribed normal displacements, the
macroscopic deformation is a uniform stretch and the constraint is

    det(F) = λ_x λ_y λ_z = 1,   λ_y = λ_z = λ_x^(−1/2).

That is exact at every load level. And because
`u_lat/u_ax = (λ_x^(−1/2) − 1)/(λ_x − 1)` is **not** constant — it drifts from
−0.500 at small strain to −0.459 at 12 % — the displacement ratio must be
recomputed each increment. So the paper's *"adjusted dynamically"* describes a
nonlinear closed-form relation, **not** a solve-time equilibrium iteration. That
reading is what makes the composition sufficient; had it meant an iterated
constraint, a new executioner would have been required.

The general split `λ_y = λ_x^(−p)`, `λ_z = λ_x^(−(1−p))` is exactly isochoric
for **any** `p`, since the exponents sum to 1. `p = 1/2` is the equal-compression
reading of "the two reciprocal directions"; unequal `p` is a sensitivity, not a
different model.

**Verified in a solve, not just in algebra**
(`decks/al2219/iso_bc_verification_cube.i`, `scripts/check_isochoric.py`):
macroscopic `det(F) = 1.000000000000`, deviation `0.0e0`, with lateral stretch
`0.94491118` against the closed-form `1.12^(−1/2)`. The volume is measured from
the **deformed node positions**, not from the prescribed values, so a BC that
failed to apply would show up rather than be assumed away — which is exactly
what happened on the first run, and is how G13 was found. Note that only the
*macroscopic* volume is constrained; the element-level Jacobian is not 1
pointwise, and the checker reports its spread rather than implying otherwise.

Two deck-level points recorded for the production decks: the roller faces must
constrain **one** component each (fully fixing them over-constrains the
transverse contraction the constraint is derived for), and `finite_strain = true`
is required — `det(F) = 1` is a finite-strain statement and 12 % axial strain is
well outside where `1 + tr(ε)` would substitute.

### G13 — FIXED — Per-BC time functions were silently ignored by `cp_quasi_static` (`84bebc7`)

Found by the G4 verification failing: the first run returned lateral stretches of
**+1.02 (expansion)** where the constraint demands 0.945 — the BC's `values` had
been applied with its `function` ignored.

`BoundaryCondition.dirichlet_spec` has always carried a `function` key, the
`dirichlet` BC documents it (*"Optional Function name whose scalar value at the
current time scales the BC value"*), and four executioners honour it —
`pf_fem_transient`, `scalar_steady`, `vpinn_steady`, `pf_pinn`. But
`cp_quasi_static` scaled every Dirichlet value by a single deck-wide factor (the
`load_function` if declared, else the load factor) and never consulted
`spec["function"]`. A deck could name a function, pass `pinnacle validate`, and
have the raw amplitude applied at every step.

**No bundled deck puts `function =` on a BC on any executioner**, which is both
why this went unnoticed and why fixing it cannot regress one.

Severity **A**: the declared capability is the only way to give two prescribed
components different time dependence, and the study needs exactly that. The
`load_function` key is not a substitute — it scales every BC alike, so it can
express the axial ramp or the lateral relation but never both at once.

A second, smaller trap was found alongside it. The `.i` grammar **tokenizes a
quoted value**, so `expression = '50.0 * t'` reached `Function` as
`[50.0, '*', 't']` and was compiled as a *list literal*, failing with
`float() argument must be a string or a real number, not 'list'` from inside
`eval` — a message naming neither the key nor the quotes. An **unquoted**
expression passes through verbatim *even when it contains spaces*, so quoting is
the one spelling that breaks, which is the opposite of the intuition the grammar
creates elsewhere. `_expression()` now raises a diagnostic naming the cause and
printing the unquoted form to use. **Study decks must never quote a `parsed`
expression.**

### G4 — original finding — No isochoric (constant-volume) loading driver

Registered BCs: `convective`, `dirichlet`, `fcci_reactive_sink`,
`gaussian_flux`, `grain_boundary`, `hard_constraint`, `interface_flux`,
`marshak`, `neumann`, `periodic`, `pf_contact_angle`, `radiative`,
`surface_radiation`, `thermal_resistance`. None imposes a volume constraint,
and no load-ratio controller exists.

The Al2219 paper's loading case (a) applies tension along the crystallographic
direction and compression on the two reciprocal directions, with the ratios
**adjusted dynamically during the run** so total volume is conserved. Case (b),
uniaxial with traction-free lateral faces, is expressible today.

**Action.** New driver component implementing the volume constraint with
in-run ratio adjustment. Its mathematical definition must be verified against
the source publication before use, per the study plan.

### G5 — B — No conforming spherical-void mesh, and no local grading

`mesh/generators.box_with_void` documents: *"Elements whose centroid falls inside
the void are dropped and the exposed surface becomes the `void_boundary` side
set"*, with `divisions` a uniform `(nx, ny, nz)`. The void surface is therefore
**staircased**, and there is no local refinement control.

Damage nucleates at the void equator, so a staircased equator perturbs precisely
the quantity the study measures, and a uniform mesh cannot deliver the
`h_min ≤ 0.2–0.5 µm` near the void that the convergence plan calls for without
an unaffordable global element count (0.25 µm uniform in a 50 µm cube is 200³).

**Route.** `hexaform` builds conforming all-hex cubed-sphere O-grid meshes with
nodes placed on the analytic sphere surface and documented second-order volume
convergence (`examples/02_single_sphere_convergence.py`), and writes Exodus
(`src/hexaform/io/exodus.py`). PINNACLE imports it via `[Mesh] type = file`,
`format = exodus`, `verify = true`, which checks Jacobian positivity and
watertightness and derives boundary sets topologically — so the void surface is
included in `all_boundary` (worked reference:
`examples/applications/cppfd_external_mesh/void_block_cppfd_3d.i`).

**Verification items for Phase 1.** (i) hexaform treats the sphere as an
*inclusion*; confirm it can emit the sphere region as a removed void with a free
surface. (ii) Confirm the refinement radius around the void is independently
controllable so the refined zone can span `r_void + 3l₀` to `r_void + 5l₀`.
(iii) Confirm element quality and orientation conventions after import.

### G6 — B — No dog-bone mesh generator

The 15 Ni specimens are dog-bone tensile bars whose thickness follows a spline
through five measured stations. No registered generator produces this; the
closest, `notched_plate_mesh` and `plate_with_hole_mesh`, are 2-D planform
generators with constant thickness.

**Action.** New registered HEX8 generator, spline type and end conditions
documented, geometric uncertainty between stations quantified, unit tests and
quality checks included. Also resolve whether the paper's 9,880 HEX8 elements
refer to the gauge region or the whole specimen.

### G7 — FIXED — No Miller-direction → Bunge-Euler utility

*Original finding below; resolved in `1f00312`. The resolution, the symmetry
verification, and the quantitative closure of gap C.6 are in the **G7 — FIXED**
entry later in this section.*

`[Crystal] euler_angles` takes three Bunge angles. The four extreme cases are
specified as loading-axis **directions** — `[001]`, `[101]`, `[111]`, and
`[-1 4 9]` — which fix only two rotational degrees of freedom. The rotation
about the loading axis is unspecified in the source and is *not* immaterial
under the isochoric condition, which constrains all three axes.

**Action.** Utility converting `[hkl]` to Euler angles under an explicit,
documented convention for the free rotation, plus a sensitivity sweep across the
admissible range. This is gap C.6 of the companion sheet.

### G8 — R — Missing study metrics

Present and usable: `curve_rmse`, `curve_nrmse`, `field_l2_error`,
`damage_field_error`, `extract_crack_path`, `crack_path_overlap`,
`crack_band_width`, `dissipated_energy`, `refinement_guidance`,
`extract_damage_outputs`, `element_volumes`, `side_set_area`, `macro_stress`,
`stress_triaxiality`, `elastic_energy`, `nodal_slip_activity`, `nodal_rho`,
`nodal_gnd_density`, `nodal_nye_tensor`, plus 47 postprocessors and 11
reductions.

Missing, and required by the predefined metric list: void volume and
void-growth ratio (needed for the 12 %-strain definition); fracture-plane fit
and fracture angle under a stated convention; Hausdorff distance between damage
isosurfaces; damage-centroid displacement; a threshold-independent field error
together with a threshold sweep (the `d = 0.65` initiation criterion must be
reported alongside its sensitivity); a fracture-initiation-strain extractor; and
a mathematical definition of propagation rate.

**Action.** Implement as PINNACLE postprocessors where they are general
(void volume, fracture-plane fit, Hausdorff, centroid) and as study-side
scripts where they encode this study's conventions.

### G9 — D — Parameter-mapping notes, no code change

`physics/slip_kinetics.py:15-26` implements
`γ̇^α = γ̇₀ |τ^α/τ_c^α|^n sign(τ^α)`, while the papers write the exponent as
`1/m_s`. Therefore `rate_exponent_n = 1/m_s = 1/0.12 = 8.333…`, not `0.12` and
not the deck default `20`. Likewise `gamma0` defaults to `1.0e-3` against the
papers' `1.0e-5 s⁻¹`. Both must be set explicitly in every study deck and
recorded in the frozen configuration.

### G12 — FIXED — Two contraction errors in the batched CP constitutive path

Found while verifying that the scalar and batched paths agree, which they must,
being two implementations of one law. Both errors are two-character `einsum`
index mistakes, and both vanish identically when the elastic deformation
gradient `F^e` is symmetric — which is why every committed test missed them.

1. **Elastic strain.** `Ee` was formed as `einsum("bji,bjk->bik", Fe_T, Fe)`.
   Since `Fe_T` is *already* the transpose, that index pattern transposes a
   second time and computes `F^e F^e`, not `(F^e)^T F^e`. The two differ by
   twice the skew part of `F^e`, so the elastic strain was contaminated by the
   elastic rotation. Measured at a state with plastic spin: **8.7 % relative
   error in `E^e`**, propagating to 3.6 × 10⁻⁴ in `S`, 6.2 × 10⁻⁴ in `τ`,
   2.3 × 10⁻³ in `γ̇`, and 6.7 × 10⁻⁴ in the Cauchy stress.
2. **Cauchy push-forward.** `einsum("bij,bjk,blk->bil", Fe, S, Fe_T)` computes
   `F^e S F^e` rather than `F^e S (F^e)^T`. With `S` symmetric the correct
   result is symmetric; the computed one was not, so **the batched "Cauchy
   stress" was not a stress**. Measured asymmetry 5.3 × 10⁻² on a stress of
   order 10² MPa.

Consequence before the fix: the scalar and batched paths drifted apart over a
load history to a relative 2 × 10⁻⁵ in `F^p` and 1 × 10⁻³ in `ρ`, and the
difference was *independent of step size* across a 16× range — the signature of
a formulation difference rather than accumulated round-off. After the fix they
agree to machine precision (3 × 10⁻¹⁶ in `F^p`, 2 × 10⁻¹⁴ in `ρ`) over 200
steps.

**Which runs were affected, and it is not a narrow set.** The batched path is
not opt-in. `CPFieldCallback.update_state` selects it whenever

```python
use_batch = self._qp_directions is None and hasattr(self.material, "compute_batch")
```

that is, for **every CP deck with a material-global orientation** — precisely
the single-crystal configuration this study benchmarks. Only decks carrying a
per-quadrature-point orientation *field* (polycrystals, EBSD-driven cases) fell
through to the scalar path and were unaffected. So any single-crystal CP or
CP–PFD result previously produced with PINNACLE carried a wrong elastic strain
and a non-symmetric stress wherever lattice rotation was present, which for
orientation-dependent slip around a void is everywhere. Consistently with this,
the two committed parity goldens that changed when the fix landed
(`rezwan_field_fem_3d`, `rezwan_field_hybrid_3d`) are exactly the two
single-orientation field decks; they encoded the defect and were regenerated.

This concerns PINNACLE results only. Both source papers were produced with
MOOSE, so their published results are outside the reach of this defect.

Fixed on `feat/pinn_benchmark` at `b49fa70`, with
`tests/unit/materials/test_cp_batch_parity.py`. The tests deliberately drive a
**non-symmetric** plastic deformation gradient; a parity test on a symmetric
state has no power against either error, which is the lesson worth carrying
into any future parity test for this material.

### G1 — FIXED — split into two defects, both now closed (`aac832b`)

The audit originally recorded G1 as "no plastic-work term in the fracture
driving force". Tracing the coupling channel shows the situation is worse, and
is two separable defects.

**G1a — the driving force is evaluated on the total strain.** The channel is
`state["strain"] → MieheSpectralSplit → PFDDrivingForceMaterial.update`, in
`cp_pfd_strategy.history_bridge_fn` for the hybrid route and the equivalent path
in `coupled_executioner` for all-FEM. What `cp_quasi_static` puts in
`state["strain"]` is `assembler._compute_strain(u_nodes)` — the **total**
kinematic strain from the displacement field. So the driving force is an
isotropic elastic energy evaluated on a strain that is largely plastic. It is
neither `ψ⁺_el` nor `W_p`: it grows with plastic deformation, which is why the
bundled decks produce plausible-looking damage localization, but its magnitude
has no thermodynamic basis. **Fix:** export the elastic Green–Lagrange strain
`E^e` per quadrature point and split that instead.

**G1b — no stored plastic work.** The papers drive damage with `a⁺ + ζ W_p` at
`ζ = 1`. Without it the ductile part of ductile phase-field damage is absent,
and `ζ` is one of the Group-C parameters the sensitivity study sweeps.

**Both are now closed** (`aac832b`). `W_p` is CP state (`b49fa70`);
`CPFieldCallback` exposes `elastic_strain()` — read from the last constitutive
evaluation, so it is exactly the measure the stress was built from — and
`plastic_work()`; `pfd_driving_force` gains `plastic_work_prefactor` defaulting
to 0.0, so the brittle limit and every existing deck are unchanged; and both
couplers read them.

Two design points are worth recording because each is a trap:

- The quantities travel under a dedicated `result["coupling"]` key rather than
  at the top level. Top-level keys of a mechanics result are treated as
  **exportable output fields**, so promoting them would have grown the
  documented `[Output] fields` catalogue for every CP deck — a user-visible
  change the coupling does not need. Putting them at the top level first was
  what broke four unrelated output-catalogue parity cases.
- The plastic work is projected onto nodes through the **same** recovery as the
  elastic energy. Since the two are summed, projecting them differently would
  give the ratio that `ζ` exists to control a mesh-dependent component that no
  parameter study could separate out.

The wiring was caught mid-error by its own guard: the first attempt set `ζ = 1`
and the run logged *"plastic_work_prefactor = 1 but the coupler supplied no
plastic work"*, because the coupled executioner whitelists which sub-solve
fields pass through and the new ones were not on the list. A silent version of
that path would have reported a ductile result computed from an elastic
driving force.

**The same defect occurred a second time, on the hybrid route, and was silent**
(`1f00312`). `cp_pfd_strategy.fem_solver_fn` builds its state dictionary as
another hand-written whitelist, and it did not name `coupling`, so the
corrected driving force was **inert on the hybrid route**: the hybrid kept
using the elastic energy of the total strain while the all-FEM route used the
corrected form, meaning the two routes were no longer solving the same model.
Nothing failed — the affected parity case was already failing on stale-hash
drift, so its result carried no information either way, and the guard could not
fire because `pfd_material.update` was still receiving a plastic-work argument
from a bridge that was itself never reached with a populated channel.

It was found by comparing the two routes directly on `cp_pfd_hybrid_2d`, which
runs at `classical_fallback_fraction = 1` and therefore **must** reproduce the
all-FEM answer exactly:

| Route | peak damage |
|---|---|
| hybrid, before the fix | 0.26218608 |
| hybrid, after the fix | 0.25491497 |
| all-FEM staggered, same physics | 0.254914971 |

Eight-significant-figure agreement after, against a 2.8 % disagreement before.
Left unfixed, that 2.8 % would have entered this study as *hybrid method
error* — a number attributed to the neural route that was really a missing
dictionary key. **This is the lesson to carry into Phase 2:** before any
hybrid-versus-FEM accuracy claim, the hybrid must first be shown to reproduce
the all-FEM answer on a configuration where the two are mathematically
identical. A route-to-route comparison that is not anchored that way cannot
distinguish method error from wiring error.

### G7 — FIXED — Miller-direction to Euler-angle utility (`1f00312`)

`orientations.direction_to_bunge(direction, loading_axis, free_rotation_deg)`
maps a crystallographic direction onto a sample axis, and decks may now write
`loading_direction` in the crystal block instead of `euler_angles` (supplying
both is refused). The alignment is the minimal rotation; `free_rotation_deg`
makes the third, undetermined degree of freedom an explicit argument.

Verified by symmetry rather than by stored numbers: a 120° free rotation about
⟨111⟩ is a cubic symmetry operation, and it reproduces the identical stress
state to 1 × 10⁻⁹. That single invariant exercises the new utility, the
framework's Euler-to-matrix convention, and the slip-system rotation together —
if any of the three disagreed, it would break.

**This closes gap C.6 quantitatively, and the answer is that the missing angle
matters.** Sweeping `free_rotation_deg` over 0–60° for `[1 1 1]` under an
asymmetric prescribed deformation:

- axial stress `σ₁₁` moves **0.06 %** (597.29 → 597.65 MPa);
- shear `σ₁₂` moves **56 %** (78.9 → 123.0 MPa);
- transverse `σ₂₂` moves 3.2 %.

For `[-1 4 9]`, which has no symmetry axis along the loading direction, `σ₁₁`
itself moves 1.8 %. So a source reporting only a loading direction leaves the
**axial stress–strain response essentially determined** but the **transverse
and shear state undetermined at the tens-of-percent level** — and it is the
transverse state that sets where damage localizes around a void. C.6 is
therefore a real uncertainty for every damage-location and fracture-angle
result, and a negligible one for stress–strain curves. Both halves of that
have to be stated wherever the four extreme orientations are reported.

One methodological note worth keeping: an early version of this check used a
transversely isotropic drive, `diag(λ, λ^(-1/2), λ^(-1/2))`, and measured
exactly zero sensitivity. That was correct — equal lateral stretches give the
loading a continuous rotational symmetry about the loading axis, so the
response *cannot* depend on rotation about it. A sensitivity study of this
degree of freedom run under symmetric loading will report that it does not
matter, which is true of that load case and false of the papers' void and
isochoric problems.

**Measured effect** on `cp_pfd_staggered_plate_2d`, peak damage:

| Configuration | peak damage | change |
|---|---|---|
| total strain, no plastic work (previous behavior) | 0.262186 | — |
| elastic strain, `ζ = 0` | 0.254915 | −2.77 % |
| elastic strain, `ζ = 1` | 0.256737 | +0.71 % |

On the Rezwan void deck the golden's `damage_max` moved from 0.675277 to
0.678039 (+0.41 %). Both are small **because these are coarse demonstration
decks at ~0.1 % strain**, where the plastic fraction of the total strain is
still small and `W_p` is correspondingly minor. At the source papers' 12 %
strain the plastic fraction dominates, so neither correction can be assumed
small there — quantifying that at the study's own load levels is part of the
convergence and sensitivity work, not something to extrapolate from these
numbers.

### G11 — B — `backend` and `device` are global, so a hybrid run cannot give each half its preferred hardware

`backend` and `device` are declared only under `[Global]`
(`schemas/pinnacle_input.schema.json`, `/properties/global/properties/backend`;
`docs/user_guide.md` §6). Measured on `ssh:rocm-server` at `float64`, the CP
constitutive update runs at 10.7 µs per quadrature point on NumPy against
106.0 µs on the GPU, while a PINN training step runs 2.9× faster on the GPU —
so the two halves of the coupled problem want opposite hardware, and a hybrid
run must pick one. NumPy has no autodiff, so hybrids are forced onto Torch and
the CP half pays a 60.7 µs per quadrature point penalty over Torch CPU.

Full derivation, break-even, and policy in `docs/02_backend_envelope.md` §3.

**Decision (user, Phase 0).** Add per-sub-executioner `backend`/`device`
overrides, so the CP half runs on NumPy and the damage half on the Torch GPU and
the hybrid is measured at its method cost rather than its backend cost. This is
a Phase 1 development item on `feat/pinn_benchmark`, registered, documented,
tested, with an ADR. Note the hazard `backend/base.py` documents explicitly: a
tensor created before a later `set_device` stays on the old device, producing a
mismatch deep inside an operator, so device-cache decorators and any
pre-construction tensor need auditing as part of this change.

### G10 — Noted — Duplicate degradation specification

`pinn_fem_hybrid_coupled_physics` takes `damage_gd_func` while the `[Materials]`
block carries `pfd_degradation`. Two places specify `g(d)`. Confirm they agree
in every hybrid deck, or the hybrid and all-FEM routes are not comparing the
same model.

---

## 3. Consequence for the study

Six blocking gaps (G1–G6) stand between the current framework and a faithful
reproduction of either source paper. G1 and G2 are physics defects that would
silently produce wrong answers rather than errors: without G1 the model is not
the papers' model, and without G2 the initial slip resistance is 2.5× too high.
Both must be closed and verified at the material-point level before any
convergence or calibration run, which is why Phase 1 places constitutive
verification ahead of every mesh study.

None of the six is research-grade difficulty. The one genuinely open item is R5,
the history-aware state surrogate, which is scoped as its own track.
