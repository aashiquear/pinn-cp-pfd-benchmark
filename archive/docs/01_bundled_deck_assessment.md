# Assessment: the bundled Rezwan 2025 application decks

Phase 0, Framework track. Audited at `feat/pinn_benchmark` (`16ff0e6`), against
`examples/applications/rezwan2025_void_cp_pfd/`.

**Conclusion: the three bundled decks are a working demonstration of the
composition architecture, not a reproduction of the published benchmark.** Every
physical parameter, the geometry, the units, and the loading mode differ from the
source paper. They are the right structural templates and the wrong numbers, so
Phase 1 builds new study decks rather than editing these.

The distinction matters for the accuracy classification the study plan requires:
the bundled decks are pinned at tier 1 (golden-reference agreement) and cannot
support tier 3 (agreement with a previously published simulation) or tier 4
(agreement with experiment) at their committed settings.

---

## 1. Parameter-by-parameter comparison

Paper values from the Al2219 source (Rezwan 2025 intragranular), Tables 1–3,
independently re-extracted in the Sources track; the paper works in
µm–MPa–s.

| Quantity | Paper | `void_field_fem_3d.i` / `..._hybrid_3d.i` | `void_pinn_constitutive_3d.i` | Assessment |
|---|---|---|---|---|
| Domain edge | 50 µm | 2 mm | 1 mm | 40× / 20× oversized |
| Void radius | 5 µm | 0.4 mm | 0.1 mm | — |
| Void-radius ratio r/L | 0.10 | **0.20** | 0.10 | field decks have the **wrong ratio**, so the ligament is not the paper's |
| C11, C12, C44 | 108, 61.3, 28.5 GPa | 108.2e9, 61.3e9, 28.5e9 Pa | 108.2, 61.3, 28.5 (comment says GPa) | values agree; C11 differs in the third digit |
| Shear modulus µ | 25.7 GPa | not set → defaults to C44 = 28.5 GPa | not set | **11 % high**, and silently defaulted |
| Reference slip rate γ̇₀ | 1e-5 s⁻¹ | 1.0e-3 | 1.0e-3 | 100× high |
| Rate exponent | m_s = 0.12, i.e. n = 1/m_s = 8.33 | `rate_exponent_n = 20` | 20 | 2.4× too rate-insensitive |
| Lattice friction τ_c | 0.112 MPa | `tau_c0 = 50.0e6` Pa = 50 MPa | not set → default 50 MPa | 446× high; conflates lattice friction with total CRSS |
| Latent hardening q | 1.4 | not set → `q_lat = 1.0` | not set | isotropic instead of latent; removes slip-system competition |
| Taylor prefactor α₀ | 0.4 | not expressible (gap G2) | — | fixed at 1.0 in code |
| Initial SSD ρ₀ | 220 µm⁻² = 2.2e14 m⁻² | `rho0 = 1.0e10` m⁻² | 1.0e10 | 2.2e4× low |
| Burgers b | 0.000286 µm = 2.86e-10 m | 2.86e-10 | 2.86e-10 | agrees |
| Storage / annihilation | k₀ = 0.3432, y_c = 0.00268 µm | `k1 = 1.0e8`, `k2 = 10.0` | same | placeholders; non-physical under either unit reading (gap G3) |
| Critical energy Gc | 22.75 MPa·µm = 22.75 J/m² | 17.0e3 J/m² | 17.0e3 | 747× high |
| Length scale l₀ | 1 µm | 0.02 mm = 20 µm | 0.02 | 20× high in absolute terms |
| l₀ / domain edge | 0.02 | 0.01 | 0.02 | — |
| Residual stiffness k | 1e-3 | not set | not set | defaults used |
| Kinetic coefficient L | 1e-3 s/µm | not set | not set | defaults used |
| Plastic-work prefactor β | 1.0 | not expressible (gap G1) | — | term absent from the formulation |
| Damage model | — | `at2` | `at2` | AT2 damages from the first increment; paper convention to be confirmed |
| Anisotropy | slip-plane resolved | `slip_plane`, weight 0.85, exponent 4.0 | same | structurally right; the two constants have no stated provenance |
| Orientation | `[001]`, `[101]`, `[111]`, `[-1 4 9]` + 93 measured | `'30 45 0'` | `'30 45 0'` | none of the paper's cases |
| Loading | isochoric (volume-conserving) and uniaxial, displacement-driven | Neumann traction 200 MPa on `xmax` | same | wrong control mode; isochoric not expressible (gap G4) |
| Mesh | conforming, density 1 at boundary / 0.65 at void | `box_with_void`, `6 6 6` uniform, staircased void | `8 8 8` | 216 elements before void removal; h = 0.333 mm ≫ l₀ |
| h vs l₀ | plan requires h ≤ l₀/2 | h/l₀ = 16.7 | h/l₀ = 6.25 | under-resolved by one to two orders |

## 2. Structural observations

**Unit systems are mixed within a single deck.** `void_pinn_constitutive_3d.i`
sets `C11 = 108.2` with a `# GPa` comment and a domain commented `mm`, while
`Gc = 17.0e3` (J/m²) and `burgers = 2.86e-10` (m) are SI. The two field decks are
consistently SI. The study decks will declare one unit system in the header and
state it in the frozen configuration, per the plan's requirement.

**The two field decks disagree with the constitutive deck on geometry.** The
void-to-domain ratio is 0.2 in the field decks and 0.1 in the constitutive deck.
Only the latter matches the paper. Any three-route comparison run on the bundled
decks as committed would be comparing different problems.

**The `mu` default is a silent trap.** `dislocation_density_cp` documents
`mu` as *"Shear modulus used in Taylor hardening. Defaults to C44."* No bundled
deck sets it, so all three use 28.5 GPa where the paper uses 25.7 GPa. Combined
with gap G2 (α₀ fixed at 1.0) the Taylor term is wrong by a factor of
`(1/0.4) × (28.5/25.7) = 2.77` before any other consideration.

**The loading story is inverted relative to the paper.** The decks apply a
Neumann traction and treat the `xmin`/`ymin`/`zmin` faces as symmetry planes.
The paper is displacement-controlled on prescribed faces with the remaining
faces fixed, and the isochoric case additionally requires in-run load-ratio
adjustment. Traction control cannot reproduce a load drop after damage
localization, which is one of the study's failure-strain definitions.

## 3. What CI actually pins

The application README states that
`tests/regression/test_rezwan_benchmark.py` pins curve RMSE, damage-field error,
and crack-path overlap, and calls this *"the reproducibility contract."*

**That file does not exist in the tree.** It is referenced from three places —
the application README, `benchmarks/rezwan2025_validation.md:61`, and
`tests/data/generate_case01_reference.py:5` — and `tests/regression/` contains
only `test_field_accuracy.py`.

What does exist, and does run, is the parity harness: `tests/parity/cases.py`
registers all three decks (`rezwan_field_fem_3d`,
`rezwan_field_hybrid_3d`, `rezwan_pinn_constitutive_3d`) with committed goldens
under `tests/parity/goldens/`, plus deck-translation goldens under
`tests/unit/framework/goldens/applications/`. So the decks are executed and
their numeric output is pinned at coarse settings — genuine tier-1 coverage —
but no physics-metric regression against published results is enforced anywhere.
`benchmarks/rezwan2025_validation.md` is consistent with this: it lists
crack-path overlap and damage-field error as *"requires external reference"* and
records that MOOSE/Roadrunner reference fields are not committed.

**Action.** Correct the stale reference in the application README and
`rezwan2025_validation.md` on `feat/pinn_benchmark`, so the tree does not claim
a contract it does not enforce. The study's own frozen-configuration decks and
metrics become the physics-level contract.

## 4. Consequence for the study

The bundled decks are reused as **composition templates only**: the variable
grouping, kernel-to-material wiring, the `degradation`/`fracture`/
`driving_force` parameter links, and the staggered and hybrid executioner block
structure are all correct and worth copying. The external-mesh deck
(`examples/applications/cppfd_external_mesh/void_block_cppfd_3d.i`) is a better
template still for the study, because it is displacement-controlled and reads an
imported mesh — both of which the study needs.

New study decks live in `decks/al2219/` and `decks/ni/` in this repository, not
in the PINNACLE tree, and carry the independently extracted parameters with
provenance recorded per value.
