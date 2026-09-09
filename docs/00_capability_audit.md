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
| G15 | A | open | — | `staggered` never increments the load and resets CP state each sub-solve — blocks Phase 2 |
| G6 | B | open | — | No dog-bone mesh generator |
| G7 | R | **fixed** | `1f00312` | No Miller-direction → Bunge-Euler utility; closes C.6 |
| G8 | R | open | — | Missing study metrics (void volume, fracture plane, Hausdorff, …) |
| G9 | D | noted | — | Exponent and reference-slip-rate mapping, no code change |
| G10 | D | noted | — | Degradation specified in two places; confirm they agree |
| G11 | B | open | — | `backend`/`device` are `[Global]`-only, so a hybrid can't split hardware |
| G12 | A | **fixed** | `b49fa70` | Two contraction errors in the batched CP path |

Ten of sixteen closed. **Four remain open**: one **A** item (G15), two **B**
items (G6, G11), and one **R** item (G8, the missing study metrics). G9 and G10
are documentation notes requiring no code.

**G15 gates Phase 2 and is the one open item that is not routine.** The other
three are development work with no unknowns in them. G15 needs load continuation
into the coupled sub-solve *and* quadrature-point state that survives across
sub-solve calls — the second is architectural and warrants an ADR. Until it
closes, no coupled result is a path-correct solution, so no solver-route
accuracy or cost comparison should be reported. The separate genuinely open
research item is R5, the history-aware CP state surrogate, scoped as its own
track rather than as a gap.

Note that the register has **grown** as work proceeded, from eleven entries to
sixteen: G12, G13, G14, and G15 were each found by verifying something else, and
none was visible from reading the documentation or the component index. That is
worth stating plainly, because it bears on how much of the framework should be
trusted without a check. Every defect found so far was **silent** — none raised
an error — and four of them were in code paths the documentation describes as
working. Three (G13, G14, G15) were found only because a verification was
written to measure the physical quantity rather than to confirm that a run
completed: the isochoric checker measures volume from deformed node positions,
and the staggered probe compares displacement fields across step counts. A check
that merely asserted "the run converged" would have passed in all three cases.

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

### G15 — OPEN, severity A — `staggered` never increments the load, and resets the CP state every sub-solve

**This is the gating blocker for Phase 2 and the most consequential finding so
far.** As it stands, the coupled executioner cannot produce a load-stepped
coupled fracture solution, so the study's primary metric — fracture-initiation
strain, reported as 7.0 % for isochoric [001] in (Rezwan 2025 intragranular) —
is not measurable, and neither is a load–displacement curve.

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

**What closing it requires.** Not a small patch:

1. Load continuation into the sub-solve. The machinery exists —
   `cp_quasi_static` already accepts `initial_load_factor` / `target_load_factor`,
   and `LoadProgram.load_factor_at_step` blends from `initial_load_factor`
   explicitly "so continuation runs start at non-zero" — so the outer loop
   passing `(prev_step_factor, step_factor)` is the intended composition. This
   also composes correctly with the per-BC functions closed in G13, since
   `t_phys = λ·t_end` then advances along the path.
2. QP state persistence across sub-solve calls. This is architectural: the CP
   state must be owned somewhere that outlives a single `run`, rather than being
   rebuilt by `make_material`. A cached callback per kernel is the obvious
   route, but it interacts with the device-cache hazard noted under G11 and with
   the assembler's ownership of the material, so it needs an ADR rather than an
   in-place edit.

Until it closes, no accuracy or cost comparison between solver routes should be
reported, because the reference itself would not be a path-correct solution.

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
