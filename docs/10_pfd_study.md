# Part II: phase-field ductile fracture on the perforated plate (J2, then CP-PFD)

Lab notebook of the PFD stages (plan: `09_pfd_benchmark_plan.md`; mesh
workflow: `workflow_pfd_plate.md`). Results enter the report only when the
author agrees; this file records everything, including what was tried and
rejected.

## 1. PINNACLE development for Part II (2026-10-13, ADR 0031)

The survey of 2026-10-13 (docs/09 §7) found six gaps between the shipped
staggered J2+PFD / CP+PFD path and a run on the perforated plate. All six
are closed in PINNACLE `feat/pinn_benchmark` under ADR 0031; what was
measured while closing them:

### 1.1 Dense against sparse damage solve

`scripts/tools/damage_solve_benchmark.py`, local workstation, AT2 system
with a random history field, Dirichlet on the two x faces, one assembly
and one direct solve per path (`results/pfd_study/damage_solve/dense_vs_sparse.csv`):

| mesh | nodes | dense: solve, operator memory | sparse: solve, memory | nonzeros | max abs difference of d |
|---|---|---|---|---|---|
| h40 hole 5 | 1,864 | 0.90 s, 141 MB (3.5 M entries) | 0.35 s (first call, SciPy import included), 6 MB | 30,444 | 3.9e-15 |
| h80 hole 5 | 6,736 | 4.3 s, 1.8 GB (45 M entries) | 0.05 s, 24 MB | 115,244 | 3.0e-15 |
| h160 hole 5 | 26,128 | not run (5.5 GB operator, 683 M entries) | 0.34 s, 94 MB | 458,508 | -- |

The two paths give the same field to round-off (the unit test pins them to
1e-10 on dense and sparse, direct and CG). The dense path scales as n^2 in
memory and n^3 in time: at h80 it costs 4.3 s per damage solve, and a
staggered run solves it up to three times per load step for 256 steps
(55 min of damage solves against 40 s sparsely); at h160 it does not fit.
`pfd_steady` now assembles sparsely at or above `sparse_threshold = 500`
nodes (every plate of this study) and densely below it, so the small
examples and goldens are unchanged.

### 1.2 What else changed, and why it mattered on this plate

- **Damage boundary set.** `pfd_steady` pinned d = 0 on EVERY node when the
  named set was missing or empty, and its default `all_boundary` is
  topological, so on the hole mesh it included the hole surface. Now
  `boundary_node_set = none` is allowed and an unknown or empty set is an
  error. The hexaform plate carries a `grips` set (the two x faces) for the
  study's damage Dirichlet condition.
- **Exodus reader.** Through `meshio` the hexaform plate came back with its
  node sets under the wrong names (`x-` had the 932 nodes of `z-`), the two
  element blocks merged under `hexahedron`, side sets expanded to every face
  of their elements, and no element variables. The reader now reads the
  Exodus schema through netCDF4; block names, exact `(element, face)` sides
  and element variables (a `grain_id` sets `mesh.elem_grain`) survive.
- **Reaction force.** `reaction_node_set` on both quasi-static executioners
  records `reaction_x/y/z` per accepted step; the engineering stress of the
  perforated plate is `reaction_x` over the gross section (50 um^2), which
  the volume-mean stress of Part I is not.
- **Viscosity and threshold.** `pfd_fracture viscosity` (eta) and
  `pfd_driving_force plastic_work_threshold` (W_c) exist, both default 0.
  eta is what the author fixes for the quasi-static study; it had to exist
  to be fixed.
- **Per-step coupled rows.** The `staggered` row of each load step now
  carries the last mechanics alternation's scalars and the damage counters
  (`max_damage`, `nodes_above_threshold`, `elements_above_threshold`,
  `initiated` at `damage_threshold = 0.8`), so the fracture-initiation step
  and the load-displacement curve are columns of `solver_history.csv`. The
  damage sub-solve receives the step-start damage and the step's time
  increment; irreversibility is enforced on d between steps (it rested on H
  alone before).
- **Element tangent assembly.** The generic FEM tangent (used by J2 and
  elasticity; the CP arm has its own batched path) evaluated a four-operand
  `einsum` without a contraction path: 0.32 s per assembly on 288 hexahedra,
  19 of 24 s of a one-step plate solve. NumPy's greedy path is now used for
  three or more operands (two-operand calls are untouched, so the CP arm's
  rounding is unchanged); the same assembly takes a few milliseconds.
- **J2 Newton "stall" was two things, neither a radius.** First, my test
  decks set the absolute Newton tolerance (1e-6, and PINNACLE's 1e-10
  default) far below the force scale of the problem (E x strain x area, 1e4
  MPa um^2 on the plate), so Newton hit the floating-point floor and reported
  non-convergence after 30 iterations; Part I's criterion (1e-3 absolute, no
  relative test) is the right one. Second, the degraded J2 tangent was a
  central difference of the split stress whose diagonal strain columns were
  twice the derivative (the diagonal perturbation was the full step, the
  quotient divided by it as a half step), so Newton on a damaged J2 problem
  converged linearly (12 iterations per sub-step at 0.125% strain) and each
  tangent cost eighteen extra return-mapping integrations (3 s per assembly
  on the h40 plate, 100 of 240 s of a one-step solve). The tangent is now
  the analytic `(I - (1 - g) P+(sigma0)) : C0` with the positive spectral
  projector verified against finite differences. The undamaged algorithmic
  tangent was checked against a finite difference of the stress on elastic
  and plastic states: exact.
- **Newton on a degraded state diverges without a line search.** The first
  coupled plate step (my deliberately brittle check, Gc = 20 MPa um, 1%
  strain) solves its first alternation in 7 iterations; the second starts
  from max d = 0.79 (mean 0.52) and the plain Newton iteration diverges from
  its first step (free-DOF residual 7e3 to 6e7 in fifty iterations) although
  the tangent is consistent. PINNACLE's Newton solver now has a backtracking
  line search (`line_search = backtracking` on both quasi-static
  executioners) and `fem_quasi_static` has load-increment cutbacks like the
  CP executioner; with the line search the same alternation converges in
  41 iterations (max d 0.97: the plate is cut at 1% strain, as that Gc
  implies).
- **Sparse Newton solve for J2.** `fem_quasi_static` now exposes
  `linear_method` and `sparse_threshold`; the h40 plate (5,592 DOFs) was
  solving its tangent densely at 2.2 s per Newton iteration under load.

Tests: `tests/unit/phase_field_damage/test_damage_sparse_solve.py` (8),
`tests/unit/mesh/test_exodus_reader_netcdf.py` (4),
`tests/unit/fem/test_reactions.py` (3); the mesh, FEM, phase-field,
coupling, executioner and solver suites pass (the VPINN steady failure
predates this work). Example: `examples/fem/phase_field_damage/ductile_j2/plate_with_hole_j2_pfd.i`.

## 2. J2 fit (docs/09 §5.1 step 1)

`scripts/j2_fit.py decks/pfd_study/j2_fit.i` → `results/pfd_study/j2_fit/`,
`figures/pfd/j2_fit.png`. Elastic constants from the Voigt average of the
CP study's cubic constants (E = 71.2 GPa, nu = 0.346), so the J2 and CP
plates share one elastic response; the digitized Vaughan curve converted
to true stress and log strain (the J2 plate runs at finite strain); plastic
strain alpha = eps_true - sigma_true / E.

| law | parameters | error at the experiment markers (2 to 17%) |
|---|---|---|
| power, sigma_y0 + K alpha^n | sigma_y0 = 299 MPa, K = 829 MPa, n = 0.44 | +1.0% at 2%, 0.0% at 7%, -3.3% at 11%, -5.0% at 13%, -6.0% at 15%, -5.6% at 17% |
| piecewise table (24 log-spaced bins of alpha, median true stress, end point anchored, non-decreasing) | 26 points | within 0.3% at every marker from 0.6% to 17% |

The power law cannot follow the nearly linear late hardening of Al2219;
the piecewise table (PINNACLE `hardening_mode = piecewise`) is the J2 law
of the study. The marker at 0.14% strain (137 MPa) sits on the digitized
elastic branch above E x strain for either law (as in Part I: that marker
reflects the test record's elastic slope, not the material's modulus).
