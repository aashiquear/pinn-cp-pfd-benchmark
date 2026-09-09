# Material-point verification, and the resolution of gap C.15

Phase 1. Harness `scripts/verify_material_point.py`; results in
`results/mp_paths.npz`, `results/mp_numerics.json`, `results/mp_c15.json`,
`results/mp_c15_rate_scan.json`. NumPy backend, `float64`, micrometre–megapascal–second.
Framework revision: `feat/pinn_benchmark` at `1c85ab2`.

Everything here is a single quadrature point: no mesh, no solver, no damage
coupling. The point is to establish that the constitutive law is the papers' law
and that it is integrated correctly, before any of those things can hide a
defect.

---

## 1. Framework corrections made first

Three gaps from the audit had to close before verification was meaningful,
because each changes the answer rather than raising an error. All are on
`feat/pinn_benchmark`, commit `1c85ab2`, with ADR 0022 and nine unit tests.

**G2, Taylor prefactor.** `HardeningLaw` had no `alpha_0` and was fixed at 1.0.
Added, defaulting to 1.0 so committed goldens cannot move, and applied in
`hardening_matrix` as well as `critical_stress` — a prefactor in the stress but
not in `h^{ab} = ∂τ_c^α/∂ρ^β` would leave the consistent tangent inconsistent
with the stress it differentiates. Verified against a central-difference
Jacobian to a relative tolerance of 1e-6.

**G3, Kocks–Mecking units.** The declared units (`k1` in 1/m, `k2`
dimensionless) contradicted the implemented law. The dimensional argument is in
ADR 0022; the conclusion is `k1` dimensionless and `k2` a length. Metadata
corrected, and `annihilation_diameter` accepted as an alternative to `k2`
with the conversion `k2 = 2 y_c` performed internally. Supplying both raises,
rather than resolving by precedence.

**A scalar/batch parity defect found en route.** `compute_batch` applied the
interaction matrix to the SSD density and then added the GND densities
*outside* the contraction, while `compute` added them inside it. The two paths
therefore disagreed whenever GND densities were present — which is the papers'
configuration. Now exact agreement (0.0 maximum difference), pinned by a test.

## 2. Kocks–Mecking mapping

The Al2219 paper's Eq. (12) is

```
d(rho_SSD^a)/dt = ( k_m sqrt(rho) - 2 y_c rho_SSD^a ) gamma_dot^a / b
```

and the framework implements

```
d(rho^a)/dt = (1/b)( k1 sqrt(rho^a) - k2 rho^a ) |gamma_dot^a|
```

giving **`k1 = k_m`** (dimensionless) and **`k2 = 2 y_c`** (length):
`k1 = 0.3432`, `k2 = 0.00536 µm` for Al2219; `k1 = 0.0295`,
`k2 = 1.58 × 10⁻⁶ mm` for Ni. Saturation is `rho* = (k1/k2)² = 4100 µm⁻²`,
verified numerically as the density at which the storage and recovery terms
balance.

Two discrepancies in the sources, both resolved and recorded:

- The Ni paper's Eq. (10) prints `2/y_c` where the Al2219 Eq. (12) prints
  `2 y_c`. With `y_c` a length the Ni form is dimensionally inconsistent
  (`L⁻³` against `L⁻¹`), so it is a typographical error and the Al2219 form is
  used for both materials.
- The storage term's density: taking it as the interaction-weighted sum, as
  written, gives a saturation density of 6.7 × 10¹⁵ m⁻², an order of magnitude
  above what is attainable in aluminum, while the per-system reading gives
  4.1 × 10¹⁴ m⁻², squarely within the observed range. The per-system reading is
  adopted, which is both the physically admissible one and the standard
  Kocks–Mecking form.

The framework accumulates with `|gamma_dot|` where the papers write
`gamma_dot`. For monotonic loading these agree; under reversal the signed form
would *decrease* a statistically-stored density, which is unphysical, so the
absolute-value form is correct and is what the papers must intend. This matters
only for the cyclic and nonproportional paths, which the surrogate training set
will contain, and is recorded as an intentional deviation from the printed
equation. Verified: the mean density is monotone non-decreasing on all four
paths, including unload–reload.

## 3. Gap C.15 resolved: the slip resistance uses the interaction-weighted sum

The Taylor term's density had two readings differing by `sqrt(1 + 11q) = 4.05`
at `q = 1.4`:

| Reading | τ_c(0), Al2219 |
|---|---|
| interaction-weighted sum (Eq. 10 with ρ from Eq. 11) | **176.71 MPa** |
| per-system SSD only (the companion sheet's check) | 43.72 MPa |

**Textual evidence.** Eq. (10) writes `sqrt(rho)` and Eq. (11) defines `rho` as
`sum_b chi^{ab}(rho_SSD^b + |rho_eGND^b| + |rho_sGND^b|)`. Read together these
are unambiguous. Decisively, the implemented SSD evolution uses the per-system
density (§2), so under the self-only reading the reported latent-hardening ratio
`q = 1.4` would appear in **no equation at all** while still occupying a row in
the papers' parameter tables. A parameter that is reported but inert is not a
plausible reading.

**Physical corroboration.** A Taylor (isostrain) polycrystal average under true
uniaxial stress, 40 Haar-uniform orientations:

| ε̇ (s⁻¹) | σ_y(0.2 %), weighted sum | σ_y(0.2 %), self-only | published FE polycrystal |
|---|---|---|---|
| 1 × 10⁻³ | 776.4 MPa | 216.8 MPa | 363.4 MPa |
| 1 × 10⁻⁴ | 598.5 MPa | 165.9 MPa | 363.4 MPa |
| 1 × 10⁻⁵ | 461.1 MPa | 127.3 MPa | 363.4 MPa |

The table is the rate scan (`results/mp_c15_rate_scan.json`, 120 load
increments). An independent repeat at 150 increments
(`results/mp_c15.json`) gives 774.7 MPa and 215.9 MPa at `ε̇ = 10⁻³ s⁻¹`, a
0.2 % difference from the tabulated values and far below the margins the
argument below turns on.

Taylor isostrain is an **upper bound** on aggregate strength, so a
finite-element polycrystal of the same material must yield at or below these
values. The self-only reading puts that upper bound at 0.35–0.60 of the papers'
own finite-element result, which would require the finite-element aggregate to
be stronger than the upper bound of the same material. It is excluded at every
strain rate plausible for a quasi-static calibration. The weighted-sum reading
is admissible at all three rates, and at `ε̇ = γ̇₀ = 10⁻⁵ s⁻¹` the ratio is
1.27 — an ordinary degree of Taylor overprediction.

The bound is one-sided by construction, so this corroborates rather than proves;
the resolution rests on the equations, with the bound test excluding the
alternative. Accuracy tier 2 (independently verified numerical agreement).

Two validations make the bound test trustworthy rather than merely arithmetic. A
prescribed isochoric stretch is *not* a uniaxial-stress test, and an earlier
version of this check used one; its apparent modulus was 53.3 GPa and its
0.2 %-offset yield was not comparable to a measured tensile curve. Solving the
lateral stretch for zero aggregate transverse stress instead recovers
**E = 71.3 GPa and ν = 0.346**, the correct elastic response for aluminum, and
the digitized calibration curve independently gives E = 72–76 GPa. Separately,
the measured rate-sensitivity exponent `d ln σ_y / d ln ε̇` is **0.113–0.116**
against the papers' reported `m_s = 0.12`, which verifies gap G9's mapping
`rate_exponent_n = 1/m_s = 8.333` numerically.

**Consequence for calibration.** Because both readings can be made to fit the
published curve by re-fitting `rho_SSD^0`, we cannot confirm that the published
`(rho_SSD^0, k_m)` were fitted under the convention now adopted. The frozen
configuration must therefore **re-calibrate** those two parameters against the
digitized curve rather than adopt the published values, and report them as
`calibrated` rather than `reported_table`. This is a change of status for two
parameters, not a change of model.

## 4. Strain-path and numerical checks

Four histories at `ε̇ = 10⁻³ s⁻¹` to 2 % strain, single crystal, cube
orientation, 400 increments:

| Path | σ₁₁ final (MPa) | mean ρ (µm⁻²) | max τ_c (MPa) | ρ monotone |
|---|---|---|---|---|
| monotonic | 478.6 | 248.9 | 188.2 | yes |
| unload–reload | 459.8 | 250.0 | 188.6 | yes |
| biaxial | 294.4 | 304.0 | 208.4 | yes |
| nonproportional (L-path) | 301.7 | 230.8 | 181.1 | yes |

All four integrate without intervention. Density accumulates monotonically in
every case, τ_c rises with it, and the biaxial path — which activates more slip
systems — hardens most, as it should.

**Time-step refinement.** The state update is explicit Euler, so first-order
convergence is the expectation. Error in σ₁₁ against the finest solution:
4.11 × 10⁻¹, 1.91 × 10⁻¹, 8.19 × 10⁻², 2.73 × 10⁻² MPa as the increment count
doubles from 200 to 1600, an observed order of 1.10–1.59. The apparent rise
toward the finest level is the usual artifact of measuring against the finest
solution rather than an exact one. At 400 increments over 1 % strain the
discretization error in σ₁₁ is 0.19 MPa, or 0.04 % — negligible against every
tolerance this study reports, but the load-increment sensitivity study still
runs separately at the field level, where path-dependent localization can
amplify it.

**Differentiability.** The central-difference Jacobian `∂σ₁₁/∂F` at a
plastically active state is step-independent to a relative 8.4 × 10⁻¹¹ under a
doubling of the perturbation, on a Jacobian of magnitude 1.1 × 10⁵ MPa. The
update is smooth and differentiable at the working point, which is what the
hybrid and surrogate routes require of it.

## 5. Note on the state of the parity harness

Ten cases in `tests/parity` fail on `master`, and the identical ten fail with
these changes applied, so the changes introduce none of them. All ten are
`fields.vtk` hash drift, including four coupled CP–PFD cases central to this
study: `cp_pfd_hybrid_2d`, `cp_phasefield_void_2d`, `cpd_field_ductile_void_3d`,
and `cpd_field_ductile_void_cylinder_3d`.

The mechanism matters. The harness compares a byte hash of a VTK file, so any
formatting or ordering change fails the case *before* its declared numeric
tolerances are applied — `cp_pfd_hybrid_2d` declares `abs = 0.1, rel = 0.1` and
never reaches them. Those four cases therefore carry no effective tier-1 numeric
protection at present, which reinforces the conclusion in
`docs/01_bundled_deck_assessment.md`: this study's own frozen configurations and
metrics have to be the physics-level contract. Refreshing the goldens, or
comparing fields numerically instead of by hash, is recorded as a separate item
and is not undertaken here.
