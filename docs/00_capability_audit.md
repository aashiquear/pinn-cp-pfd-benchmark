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

### G1 — B — The fracture driving force carries no plastic-work contribution

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

### G2 — B — The Taylor hardening prefactor α₀ is not implemented

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

**Action.** Add an `alpha0` parameter defaulting to `1.0` so committed goldens
and parity tests are unaffected; apply it in `critical_stress` **and** in
`hardening_matrix` (`h^{αβ} = ∂τ_c^α/∂ρ^β`) so the consistent tangent stays
consistent; set `alpha0 = 0.4` in every study deck.

### G3 — B — Declared units of `k1`/`k2` contradict the implementation

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

**Action.** Correct the units metadata; add a parameterization option accepting
`(k0, y_c)` directly with documented conversion; document `ρ*`.

### G4 — B — No isochoric (constant-volume) loading driver

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

### G7 — R — No Miller-direction → Bunge-Euler utility

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
