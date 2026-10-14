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

## 3. J2+PFD on the h40 perforated plate: probe and Gc(l0) calibration

Host timing probe (`decks/pfd_study/j2/probe_h40.i`, first 16 of 256 steps,
torch CPU, 2 threads): 6.2 s per staggered step with 3 alternations
(Newton 3 to 22 iterations per mechanics solve as damage grows), so a full
256-step run is about 27 minutes and three bisections side by side fit the
host with room for four more cases. At Gc = 300 MPa um, l0 = 5 um, zeta = 1
the maximum damage is already 0.46 at 1.25% strain and the driving force at
the hole edge 106 MPa (plastic work dominated: the strain concentration at
the hole puts the local plastic strain far above the nominal), so the
calibrated Gc sits in the thousands of MPa um; the bisection bracket was
widened to [10, 20000].

Alternation tolerance: the staggered test (max |d_k - d_{k-1}| < tolerance
between alternations) is not met at 1e-4 within 3 alternations once damage
grows (the per-step damage increment is 0.04); the study runs 3
alternations with tolerance 1e-3 and records `converged` per step, and the
sensitivity stage checks 5 alternations against 3 on one case.

Finite-strain against small-strain J2 on the plate (local, the 16-step
probe at Gc = 300, l0 = 5, numpy): identical damage (0.463 against 0.483
maximum at 1.25%), reactions within 0.7%, Newton 2 to 17 iterations per
sub-step for both, growing with the damage. The finite-strain tangent had
lacked the Hencky chain term (fixed, ADR 0031 item 10); after the fix the
growth of the iteration count with damage is the same in both kinematics,
so it is a property of the degraded, softening state (the line search
handles it), not of the strain measure. The study keeps finite strain, as
Part I did for crystal plasticity.

Calibration, first bisection trial (Gc = 447 MPa um, the log midpoint of
[10, 20000]): fracture at 2.7% nominal strain for l0 = 5 and 10 (l0 = 10:
the mechanics stopped converging at max d = 0.78, counted as fracture at
the stop strain; l0 = 5 reached d = 0.99 with the mechanics converging).
So Gc for a 17% fracture strain lies well above 447.

### 3.1 Why the first calibration trials failed, and the fix

Every trial of the first calibration ended in a mechanics failure at
moderate damage (max d 0.5 to 0.6, at 10.8%, 10.7% and 19.2% strain for
l0 = 20, 10 and 5 at Gc = 2990 MPa um), and the per-alternation Newton
count had climbed to 20 to 40 by 5% strain while the damage was still
below 0.1. Diagnostics (local, numpy, h40 hole mesh): the plate without
damage converges in 4 iterations per step to 10% strain (finite or small
strain alike); the hardening law (power, table, smoothed split) changes
nothing; the J2 tangent, degraded or not, matches a finite difference of
the stress at 30% strain. The cause was in the driver: the mechanics
sub-executioner started every alternation from zero displacement, so each
was a cold solve of the whole plastic problem at full load. The executioners
now take the last converged field from the coupled driver (PINNACLE ADR
0031 item 11); the calibration was restarted with it.
With the warm start the 16-step probe (Gc = 300, l0 = 5, local numpy) takes
1 to 5 Newton iterations per alternation instead of 2 to 17, with the
damage and reaction histories unchanged to the printed digits. The
calibration (`decks/pfd_study/j2/cal_h40.i`) was relaunched on the host at
02:40 on 2026-10-14; the first attempt's log is kept as
`results/pfd_study/j2/cal_h40_attempt1.log` on the host.

### 3.2 Calibrated Gc(l0) at h40 (J2)

`results/pfd_study/j2/cal_h40/` (calibration.csv: every trial; calibrated.json;
one directory per trial with its pipeline deck, PINNACLE deck, history,
results.csv, fracture.json and field series). Bisection in log Gc between
10 and 20,000 MPa um, accepted when the fracture strain (first step with a
node at d >= 0.8) is within one load step (0.08%) of 17.2%:

| l0 (um) | l0/h | Gc (MPa um) | Gc l0 (MPa um^2) | fracture strain | trials |
|---|---|---|---|---|---|
| 5 | 2 | 2063 | 10,300 | 17.34% | 6 |
| 10 | 4 | 2696 | 27,000 | 17.27% | 5 |
| 20 | 8 | 2818 | 56,400 | 17.19% | 4 |

Gc grows only weakly with l0 (x1.37 for x4 in l0), so Gc l0 is far from
constant: the driving force is the plastic work at the hole edge, and the
damage band width set by l0 changes the initiation strain less than it
would for an elastic driving force. The base case of the sweeps is l0 = 5
um (2h, the brittle rule's lower bound) with Gc = 2063.

Figures: `figures/pfd/j2_cal_h40_curves.png` (calibrated curves with the
initiation point and the experiment markers), `j2_cal_h40_bisection.png`,
`j2_cal_h40_gc_l0.png`; table `results/pfd_study/j2/cal_h40/j2_cal_h40_table.csv`.

**For the author.** The gross-section engineering stress of the perforated
plate peaks early (415 MPa at 6.1% for l0 = 5, 388 at 3.6% for l0 = 10,
368 at 1.8% for l0 = 20) and softens from there to the initiation point:
with AT2 (no elastic threshold) and the plastic-work driving force, the
damage is diffuse long before any node reaches 0.8 (maximum d grows from
the first plastic step; the mean damage is 0.1 to 0.2 at 10%), and the
larger l0 spreads it wider. The calibration target was the initiation
strain only, so this is the model's behaviour, not a fit error. Two
sensitivity cases are queued to quantify it: `damage_model = at1` (an
elastic threshold before damage starts) and `zeta = 0.5`, besides the
W_c threshold cases of `sens_h40.i`. The experiment markers in the figure
are the plain specimen's and are not a target for the perforated plate.

## 4. J2+PFD sweeps at h40 (docs/09 §5.1 steps 3-4)

### 4.1 Length scale at fixed Gc (the Fig. 8 analogue)

`decks/pfd_study/j2/sweep_l_h40.i`: Gc = 2063 MPa um (calibrated at l0 = 5),
l0 = 2.5, 5, 10, 20, 40 um on the h40 mesh (l0/h = 1 to 16), 256 steps to
20%, torch CPU 2 threads, 13-20 min per case. Figures
`figures/pfd/j2_sweep_l_h40_{curves,panels}.png`; summary
`results/pfd_study/j2/sweep_l_h40/sweep_l_h40_summary.csv`.

| l0 (um) | l0/h | peak stress (MPa) at strain | fracture strain (d >= 0.8) | energy to 20% (MPa) | steps / Newton |
|---|---|---|---|---|---|
| 2.5 | 1 | 485 at 10.8% | none by 20% (max d 0.57) | 90.1 | 256 / 754 |
| 5 | 2 | 415 at 6.1% | 17.3% | 74.6 | 256 / 872 |
| 10 | 4 | 373 at 2.0% | 12.7% | 49.0 | 212 (stopped at d 0.95) / 844 |
| 20 | 8 | 356 at 1.2% | 12.3% | 41.0 | 231 (stopped) / 988 |
| 40 | 16 | 353 at 1.2% | 14.8% | 42.7 | 256 / 1149 |

As in (chakraborty 2016 multiscale) Fig. 8 the peak stress and the fracture
energy fall monotonically with l0 and saturate (356 to 353 MPa, 41 to 43
MPa between l0 = 20 and 40); the fracture strain falls from "none" at
l0 = h to 12.3% at l0 = 8h and then rises again at l0 = 16h, where the band
(2 l0 = 80 um) exceeds the ligament and the damage is no longer localised
at the hole. At l0 = h (the brittle rule violated) the damage cannot
localise on the mesh at all. The study therefore uses l0 >= 2h and the
mesh sweep (4.2) decides between 2h and 4h.

### 4.2 Parameter sensitivity at h40, l0 = 5 um (one at a time)

`decks/pfd_study/j2/sens_h40.i`, 9 cases plus the base, 7 side by side at 2
threads (12-23 min each); figures `figures/pfd/j2_sens_h40_{curves,panels}.png`,
summary `results/pfd_study/j2/sens_h40/sens_h40_summary.csv`.

| case | fracture strain | peak stress (MPa) at strain | energy to 20% (MPa) | note |
|---|---|---|---|---|
| base (Gc 2063, l0 5, eta 0, W_c 0, 3 alternations) | 17.34% | 415 at 6.1% | 74.6 | |
| Gc x 0.5 (1032) | 7.58% | 368 at 1.8% | 28.6 (to 7.6%) | stopped at d 0.95 |
| Gc x 2 (4127) | none by 20% (max d 0.45) | 490 at 10.9% | 90.8 | |
| l0 2.5 (= h) | none by 20% (max d 0.57) | 485 at 10.8% | 90.1 | band unresolved |
| l0 10 (= 4h) | 12.66% | 373 at 2.0% | 49.0 | |
| eta 1 MPa s | 17.34% | 416 | 74.7 | identical to the base |
| eta 10 MPa s | 17.66% | 418 | 75.3 | +0.3% |
| W_c 5 MPa | 17.27% | 432 at 5.9% | 76.7 | less early softening |
| W_c 20 MPa | 17.66% | 486 at 6.2% | 82.9 | peak +17%, initiation +0.3% |
| 5 alternations | 17.34% | 415 | 74.6 | identical to 3 (tolerance met within 3) |

Fracture initiation is governed by Gc and l0 (half Gc halves the fracture
strain; l0 = 4h moves it from 17.3 to 12.7%); the fixed viscosity is
immaterial at the quasi-static rate (eta = 1 MPa s changes nothing, eta =
10 by one load step); the plastic-work threshold W_c moves the initiation
strain by at most 0.3% while it removes most of the early diffuse softening
(peak 486 against 415 MPa at W_c = 20 MPa): it is the lever for the author's
decision on the softening noted in §3.2. The alternation count does not
matter: the 1e-3 test is met within 3 alternations at every step.

### 4.4 Second sensitivity set: AT1 and the plastic-work prefactor

`decks/pfd_study/j2/sens2_h40.i`, at the calibrated Gc = 2063 and l0 = 5:
`damage_model = at1` never damages by 20% (max d = 0.0; the gross-section
curve hardens to 681 MPa like a plain plate, since AT1's elastic threshold
Gc/(l0) sits above the driving force everywhere at this Gc); `zeta = 0.5`
reaches max d = 0.48 by 20% with no initiation and a peak of 485 MPa at
10.8%. Either variant would need its own, lower, Gc to initiate at 17.2%;
both remove the early diffuse softening of the AT2 base (§3.2). Figures
`figures/pfd/j2_sens2_h40_{curves,panels}.png`. The choice between AT2
with a plastic-work threshold (§4.2: W_c = 20 MPa keeps the initiation
strain and restores the peak) and AT1 is the author's.

## 5. CP+PFD on the h40 perforated polycrystal plate

### 5.1 Cost of the coupled CP solve, and a fix

The first CP+PFD calibration trials ran at 28 s per staggered step on the
host (torch CPU, 2 threads), against 0.07 s per Newton iteration for the
CP arm alone (report §10). A local profile of two steps (numpy, 3 threads,
206 s) put 160 s in `load_state_dict`: the coupled driver restores the
step-start quadrature-point state before every alternation, and the CP
callback's snapshot was a nested Python list rebuilt point by point with
a conversion of the whole list per point and key (261,632 `np.asarray`
calls). With arrays in the snapshot and one conversion per key (PINNACLE
50fd7ed, ADR 0031 item 12) the same two steps take 45 s, of which the
mechanics solves are 43 s. The calibration was restarted with the fix.

### 5.2 Calibrated Gc(l0) at h40 (CP): FIRST RESULT, VOID

**Retracted (2026-10-14 08:00).** The CP+PFD plate carried only 150 MPa at
0.3% strain and 308 MPa at 17%, half the J2 plate and half the Part I CP
plate's yield: inside a coupled load window `cp_quasi_static` built its
load program with the path's full t_end (200 s) over the window's single
sub-step, so the rate-dependent slip law saw a 200 s increment for a
0.78 s step (strain rate 256 times too low; flow stress halved, n = 8.33).
Fixed in PINNACLE (ADR 0031 item 13: the window's duration is
t_end (target - initial); history rows carry dt; test). The numbers below
are kept only as the record of the wrong-rate run (`results/pfd_study/cp/
cal_h40_wrong_dt` on the host); the calibration was rerun (§5.3).


`decks/pfd_study/cp/cal_h40.i`, same bisection as the J2 one (log Gc in
[10, 20000], accepted within one load step of 17.2%), 6-7 trials per l0,
4-13 min each with the state-restore fix (3 trials side by side, torch
CPU 2 threads), every step converged:

| l0 (um) | l0/h | Gc (MPa um) | Gc l0 (MPa um^2) | fracture strain | J2 Gc for the same l0 |
|---|---|---|---|---|---|
| 5 | 2 | 519 | 2,590 | 17.27% | 2063 |
| 10 | 4 | 847 | 8,470 | 17.19% | 2696 |
| 20 | 8 | 1302 | 26,000 | 17.11% | 2818 |

The crystal-plasticity plate needs a quarter to a half of the J2 Gc to
initiate at the same strain: the plastic work that drives the damage
concentrates in the favourably oriented grains at the hole (the hole sits
on a four-grain junction, docs/09 §3), so the local driving force at the
initiation site is larger than the J2 continuum's. Gc rises faster with
l0 than for J2 (x2.5 for x4 in l0 against x1.37): a wider band averages
over more grains and loses the intragranular peak. Figures
`figures/pfd/cp_cal_h40_{curves,bisection,gc_l0}.png`, table
`results/pfd_study/cp/cal_h40/cp_cal_h40_table.csv`.

### 4.5 Author's decision on the J2 configuration (2026-10-14)

The J2+PFD plate carries AT2 with the plastic-work threshold W_c = 20 MPa
(§4.2: peak restored, initiation strain moved by 0.3%). Because the
threshold moves the initiation strain, the J2 Gc(l0) calibration, the
length-scale sweep, the mesh sweep and the sensitivity set are repeated
with it (`results/pfd_study/j2/*_wc20`); the W_c = 0 pass above is the
comparison. CP+PFD keeps W_c = 0: its curves harden to initiation.

### 5.3 Calibrated Gc(l0) at h40 (CP), corrected time increment

`decks/pfd_study/cp/cal_h40.i` rerun with PINNACLE f3eeae7 (results
`results/pfd_study/cp/cal_h40/`, 6-9 trials per l0, 4-13 min each):

| l0 (um) | l0/h | Gc (MPa um) | fracture strain | J2 Gc (W_c = 0) |
|---|---|---|---|---|
| 5 | 2 | 1106 | 17.19% | 2063 |
| 10 | 4 | 1832 | 17.11% | 2696 |
| 20 | 8 | 2818 | 17.27% | 2818 |

With the right strain rate the crystal-plasticity plate needs about half
the J2 Gc at l0 = 5, two thirds at l0 = 10 and the same at l0 = 20: the
intragranular plastic work at the hole edge is the stronger driving force
when the band is narrow enough to see one grain, and the difference
vanishes once the band averages over the four grains at the junction.
Figures `figures/pfd/cp_cal_h40_{curves,bisection,gc_l0}.png` (replacing
the void ones), table `results/pfd_study/cp/cal_h40/cp_cal_h40_table.csv`.

Observation on the corrected CP curves (`cp_cal_h40_curves.png`): the
perforated CP plate's gross-section engineering stress follows the plain
specimen's experiment markers within 1 to 2% from 0.6% to 17% (peak 596
MPa at 16.6% for l0 = 5, against 593 at 17% in the experiment) and drops
at initiation; the J2 plate with the same hole peaks at 415 MPa at 6%
(W_c = 0). The agreement with the plain-specimen markers is not a
validation: the perforated plate's gross section is 10% smaller and the
Part I plain CP plate matched the same markers, so the two effects (net
section against strain concentration and the different hardening of the
CP polycrystal under the hole's stress state) happen to cancel on this
geometry; it is recorded as a tier-4 agreement of a different specimen.

### 5.4 CP length scale at fixed Gc (the Fig. 8 analogue)

`decks/pfd_study/cp/sweep_l_h40.i`: Gc = 1106 MPa um (calibrated at
l0 = 5), l0 = 2.5 to 40 um on the h40 mesh, 3 cases side by side.
Figures `figures/pfd/cp_sweep_l_h40_{curves,panels}.png`; summary
`results/pfd_study/cp/sweep_l_h40/sweep_l_h40_summary.csv`.

| l0 (um) | l0/h | peak stress (MPa) | fracture strain | energy to the end (MPa) |
|---|---|---|---|---|
| 2.5 | 1 | 658 | none by 20% | 103 |
| 5 | 2 | 596 | 17.2% | 85.0 |
| 10 | 4 | 486 | 10.5% | 43.7 |
| 20 | 8 | 418 | 7.1% | 26.7 |
| 40 | 16 | 407 | 6.9% | 25.4 |

The crystal-plasticity plate gives the clean form of the reference's
Fig. 8: peak stress, fracture strain and fracture energy all fall
monotonically with l0 and saturate between 8h and 16h (418 to 407 MPa,
7.1 to 6.9%). The J2 plate's rise of the fracture strain at l0 = 16h
(§4.1) does not occur here: with the intragranular driving force the
widest band still localises at the hole. At l0 = h no initiation on either
material: the brittle rule l0 >= 2h holds for both couplings, and the
mesh sweep at fixed l0 (J2 §4.3; CP to follow on the h80 mesh) decides
between 2h and 4h.

### 5.5 CP parameter sensitivity at h40, l0 = 5 um (one at a time)

`decks/pfd_study/cp/sens_h40.i`, 9 cases plus the base (Gc 1106, l0 5,
eta 0, W_c 0, 3 alternations), 3 side by side at 2 threads; figures
`figures/pfd/cp_sens_h40_{curves,panels}.png`, summary
`results/pfd_study/cp/sens_h40/sens_h40_summary.csv`.

| case | fracture strain | peak stress (MPa) at strain | energy (MPa) | max d at the end |
|---|---|---|---|---|
| base | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| Gc x 0.5 (553) | 8.59% | 472 at 8.3% | 35.1 | 0.99 |
| Gc x 2 (2212) | none by 20% | 657 at 20.0% | 103.1 | 0.46 |
| l0 2.5 (= h) | none by 20% | 658 at 20.0% | 103.3 | 0.57 |
| l0 10 (= 4h) | 10.47% | 486 at 9.5% | 43.7 | 0.96 |
| eta 1 MPa s | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| eta 10 MPa s | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| W_c 5 MPa | 17.66% | 602 at 17.1% | 88.2 | 0.98 |
| W_c 20 MPa | 19.06% | 619 at 18.5% | 97.9 | 0.96 |
| 5 alternations | 17.11% | 595 at 16.6% | 84.1 | 0.97 |

As for J2, Gc and l0 govern initiation (half Gc: 8.7%; l0 = 4h: 10.5%;
l0 = h and twice Gc: no initiation by 20%) and the viscosity is immaterial
at the quasi-static rate (eta = 1 and 10 MPa s: identical fracture strain).
Unlike J2, the plastic-work threshold moves the CP initiation strain
(W_c = 5: +0.5%, W_c = 20: +2.0%): with the intragranular driving force
the plastic work at the initiation site is spent in a narrower range, so a
threshold of 20 MPa removes a larger share of it. Five alternations
against three change the fracture strain by one load step (17.11 against
17.19%) and the final damage by 0.01: the 1e-3 alternation test is met
within three. The CP configuration therefore keeps W_c = 0 (§4.5).

## 6. Fields of the calibrated plates (J2 at W_c 0, CP)

`scripts/pfd_fields.py decks/pfd_study/report/fields_cal_h40.i` →
`figures/pfd/cal_h40_damage.png` (nodal damage on the deformed plate at 5,
10, 15, 17.2 and the last snapshot), `cal_h40_crack.png` (elements with
mean d >= 0.8 at the end), `results/pfd_study/cal_h40_fields.csv`.

| strain | CP max / mean d | J2 max / mean d |
|---|---|---|
| 5% | 0.22 / 0.14 | 0.17 / 0.09 |
| 10% | 0.42 / 0.27 | 0.44 / 0.18 |
| 15% | 0.61 / 0.39 | 0.66 / 0.23 |
| 17.2% (initiation) | 0.81 / 0.44 | 0.78 / 0.25 |

The two couplings localise differently. J2 forms a transverse band
through the hole from 10% on (the classic crack path across the
ligament) with the far field at d < 0.2. CP damages the whole plate
(mean 0.44 at initiation, 0.6 in the grains adjacent to the hole) and
initiates in the element row at the hole's equator; the band after
initiation runs along the grain structure rather than straight across.
**For the author:** the CP plate's diffuse damage is the plastic-work
driving force acting in every grain from the first plastic step (as the
J2 softening of §3.2 was); the CP sensitivity (§5.5) showed a threshold
W_c = 20 MPa moves the CP initiation strain by +2%, so a recalibration
with W_c > 0 would localise the CP damage at the price of a new Gc.
The equivalent plastic strain is not in the coupled runs' snapshots
(the staggered writer carries displacement, damage, stress, strain and
the driving force); the plastic-strain row of the figure is empty until
the exporter carries it.

### 4.6 J2 calibration with W_c = 20 MPa (the carried J2 configuration)

`decks/pfd_study/j2/cal_h40.i` (second pass, results
`results/pfd_study/j2/cal_h40_wc20/`), same bisection and acceptance:

| l0 (um) | Gc, W_c = 20 (MPa um) | fracture strain | Gc, W_c = 0 |
|---|---|---|---|
| 5 | 2003 | 17.11% | 2063 |
| 10 | 2466 | 17.27% | 2696 |
| 20 | 2394 | 17.19% | 2818 |

The threshold lowers the calibrated Gc by 3% at l0 = 5 and by 9-15% at
the wider bands, and flattens Gc(l0): with the first 20 MPa of plastic
work withheld, the wide band's driving force loses its early diffuse
part and needs less toughness to reach the same initiation strain. The
W_c = 20 sweeps and sensitivity (`*_wc20`) run from Gc = 2003, l0 = 5.
The peak gross-section stress of the calibrated l0 = 5 plate rises from
415 MPa (W_c 0) to 485 MPa, as the sensitivity case predicted, and the
energy to 20% from 74.6 to 81.6 MPa. Figures
`figures/pfd/j2_cal_h40_wc20_{curves,bisection,gc_l0}.png`.

### 4.3 J2 mesh sensitivity at fixed l0 = 5 um (W_c 0 pass)

`decks/pfd_study/j2/sweep_h_l5.i` (first pass; results
`results/pfd_study/j2/sweep_h_l5/`): the calibrated h40 case (Gc 2063)
against the h80 hole mesh (l0/h = 4), same deck. h160 (26,128 nodes) was
started and killed after 5.5 h without a completed load step (sparse LU of
78k DOFs per Newton iteration, three alternations per step: days per
path); it is left to a scaling study.

| mesh | h (um) | l0/h | fracture strain | peak stress (MPa) | energy to 20% (MPa) |
|---|---|---|---|---|---|
| h40 | 2.5 | 2 | 17.34% | 415 | 74.6 |
| h80 | 1.25 | 4 | 17.42% | 412 | 72.5 |

One load step (0.08%) in the initiation strain, 0.7% in the peak and 3%
in the energy between l0 = 2h and l0 = 4h: the l0 >= 2h rule of brittle
phase-field fracture holds for the J2 coupling on this plate, and h40 is
the production mesh of Part II as it was of Part I. The W_c = 20 pass
repeats the comparison (§4.7).

## 5.6 CP mesh sensitivity at fixed l0 = 5 um

`decks/pfd_study/cp/sweep_h_l5.i`, Gc 1106 (results
`results/pfd_study/cp/sweep_h_l5/`, figures
`figures/pfd/cp_sweep_h_l5_{curves,panels}.png`):

| mesh | h (um) | l0/h | fracture strain | peak stress (MPa) | energy (MPa) |
|---|---|---|---|---|---|
| h40 | 2.5 | 2 | 17.19% | 596 | 85.0 |
| h80 | 1.25 | 4 | 17.58% | 587 | 86.3 |

Five load steps (0.4%) in the initiation strain and 1.4% in the peak:
larger than for J2 (the intragranular driving force is sharper at the
hole, and the finer mesh resolves the hot spot at the grain junction
better) but within the 2% the Part I grain-instance check accepted. The
production mesh stays h40 with l0 = 5 um (l0 = 2h) for both couplings.

## 5.7 CP grain-instance check (Part II analogue of Part I §7)

`decks/pfd_study/cp/grains_h40.i`: the calibrated case (Gc 1106, l0 5,
h40) on the two other Part I grain structures (instances sB and sC: new
Voronoi seedings grown the same way, new orientation draws), everything
else unchanged. Figures `figures/pfd/cp_grains_h40_{curves,panels}.png`.

| instance | fracture strain | peak stress (MPa) | energy (MPa) | max d at the end |
|---|---|---|---|---|
| A (calibration) | 17.19% | 596 | 85.0 | 0.98 |
| sB | none by 20% | 613 | 99.8 | see csv |
| sC | 16.09% | 612 | 83.4 | 0.98 |

The grain configuration at the hole governs the initiation strain: with
the same Gc and l0, instance sC initiates 1.1% earlier and instance sB
not at all within the 20% path, while the gross-section curves of the
three instances agree within 3% up to 15% (as the Part I plain plates
did, within 2%). This is the orientation control the plan anticipated
(docs/09 §3: the plastic zone of the hole is smaller than a grain). A
Gc calibrated on one instance is therefore an instance-specific number;
the report carries the three instances side by side and the calibration
as the reference instance's, and the author may prefer to calibrate
against the mean of several instances (three runs of 13 min per trial).
