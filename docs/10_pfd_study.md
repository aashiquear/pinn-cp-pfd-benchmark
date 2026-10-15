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

> Cleanup of 2026-10-15 (author's request): the figures and figure decks of the
> W_c 0 and W_c 20 J2 passes (§3.2, §4.1-4.9), the AT1/zeta set (§4.4), the
> probe deck and the h160 mesh were removed from the repository; the numbers
> above stay as the calibration record, the carried J2 results are §4.10-4.13.

### 4.10 W_c scan: the carried J2 configuration (author's revision of 2026-10-14)

The author read the W_c = 0 and W_c = 20 calibrations (§3.2, §4.6) as
softening too early (onset at 6% strain while the experiment hardens to
17%) and their Gc (2000 MPa um = 2 kJ/m^2) as too high for an aluminium
alloy, and asked for a threshold that delays the onset past most of the
hardening with a lower Gc. Decks `decks/pfd_study/j2/cal_h40_wc{200,400,600,800}.i`
(bases `base_h40_wc*.i`), one log-bisection of Gc each at l0 = 5 um,
2026-10-14 23:30 - 2026-10-15 01:00 (4-9 trials of 4-16 min); figure
`figures/pfd/j2_wc_scan_h40_{curves,panels}.png` by
`scripts/pfd_wc_scan_figures.py decks/pfd_study/j2/wc_scan_h40_figures.i`,
table `results/pfd_study/j2/j2_wc_scan_h40_table.csv`. The scan values
come from the plastic-work density at the hole of the W_c 20 run (max
driving force 91 MPa at 5%, 428 at 10%, 928 at 15%, 1160 MPa at 17%).

| W_c (MPa) | calibrated Gc (MPa um) | fracture strain | softening onset (strain at peak) | peak stress (MPa) | energy (MPa) |
|---|---|---|---|---|---|
| 0 | 2063 | 17.34% | 6.1% | 415 | 74.6 |
| 20 | 2003 | 17.11% | 6.1% | 485 | 81.6 |
| 200 | 385.5 | 17.19% | 13.9% | 571 | 87.0 |
| 400 | 155.9 | 17.19% | 14.1% | 485 | 75.1 |
| 600 | 144.7 | 17.19% | 13.4% | 472 | 73.8 |
| 800 | 144.7 | 17.19% | 13.4% | 472 | 73.8 |

Findings. (i) The threshold and Gc trade against each other: raising W_c
from 20 to 200 MPa moves the onset from 6.1 to 13.9% and lets Gc fall
5x (2003 -> 386), with the peak rising from 485 to 571 MPa (97% of the
experiment's 590). (ii) Beyond that the trade turns: at W_c = 400 the onset is one step later
(14.1%) but the peak is already down to 485 MPa; at W_c = 600 the
plastic work exceeds the threshold only in the root elements late in the
path, the bisection needs Gc = 145 to still initiate by 17.2%, and the low
Gc degrades the root so fast that the peak falls back to 472 MPa with the
onset at 13.4%; W_c = 800 reproduces W_c = 600 to the digit because the
work above 600 MPa and above 800 MPa accrue at the same elements in the
same steps (the threshold is already above the plastic work everywhere
except the root). (iii) The carried J2 configuration is therefore
W_c = 200 MPa, Gc = 385.5 MPa um (= 0.39 kJ/m^2), l0 = 5 um: the J2
plate then follows the hardening to 13.9% and the markers to within 3%
at the peak, against the CP plate's 16.6% onset and 596 MPa. What the J2
route still cannot do is reproduce the CP plate's late, sharp loss of
load (the J2 curve softens over 3% of strain before initiation while CP
softens over 0.6%): with a homogeneous flow stress the plastic work at
the hole spreads over a wider zone, so the damage that the threshold
delays then grows over a larger volume. The J2 l0 sweep, sensitivity (four
values per parameter) and mesh study (h20, h40, h80) are rerun at this
configuration (§4.11-4.13); the W_c 0 and W_c 20 results of §3.2, §4.1-4.9
are superseded and kept only as the calibration record.

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

`decks/pfd_study/cp/sens_h40.i`, 17 cases plus the base (Gc 1106, l0 5,
eta 0, W_c 0, 3 alternations), four values per parameter (the author's
2026-10-14 request; the x0.25 / x4 Gc, l0 3.75 / 7.5, eta 3 / 30 and
W_c 10 / 40 cases were added on 2026-10-15 at 3 side by side, 2 threads);
figures `figures/pfd/cp_sens_h40_{curves,panels}.png` (one colour per
parameter, one line style per value), summary
`results/pfd_study/cp/sens_h40/sens_h40_summary.csv`.

| case | fracture strain | peak stress (MPa) at strain | energy (MPa) | max d at the end |
|---|---|---|---|---|
| base | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| Gc x 0.25 (277) | 4.22% | 381 at 3.9% | 14.6 | 1.00 |
| Gc x 0.5 (553) | 8.59% | 472 at 8.3% | 35.1 | 0.99 |
| Gc x 2 (2212) | none by 20% | 657 at 20.0% | 103.1 | 0.46 |
| Gc x 4 (4424) | none by 20% | 666 at 20.0% | 104.0 | 0.27 |
| l0 2.5 (= h) | none by 20% | 658 at 20.0% | 103.3 | 0.57 |
| l0 3.75 (= 1.5h) | none by 20% | 641 at 20.0% | 102.0 | 0.74 |
| l0 7.5 (= 3h) | 12.89% | 529 at 12.0% | 57.8 | 1.00 |
| l0 10 (= 4h) | 10.47% | 486 at 9.5% | 43.7 | 0.96 |
| eta 1 MPa s | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| eta 3 MPa s | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| eta 10 MPa s | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| eta 30 MPa s | 17.19% | 596 at 16.6% | 85.0 | 0.98 |
| W_c 5 MPa | 17.66% | 602 at 17.1% | 88.2 | 0.98 |
| W_c 10 MPa | 18.13% | 608 at 17.6% | 91.4 | 0.98 |
| W_c 20 MPa | 19.06% | 619 at 18.5% | 97.9 | 0.96 |
| W_c 40 MPa | none by 20% | 641 at 20.0% | 103.0 | 0.68 |
| 5 alternations | 17.11% | 595 at 16.6% | 84.1 | 0.97 |

Gc and l0 govern initiation and do so monotonically over the four values:
Gc x 0.25 / x 0.5 bring initiation to 4.2 / 8.6%, Gc x 2 / x 4 remove it
from the path (max d 0.46 / 0.25 at 20%); l0 = 3h / 4h bring it to 12.9 /
10.5%, l0 = h / 1.5h remove it. The viscosity is immaterial at the
quasi-static rate (eta = 1 to 30 MPa s: the same fracture strain, the
peak within 0.06 MPa). The plastic-work threshold moves the CP initiation
strain by about +0.1% strain per MPa of W_c (5 / 10 / 20 MPa: +0.5 / +0.9
/ +1.9%) and W_c = 40 MPa removes initiation from the path; with the
intragranular driving force the plastic work at the initiation site is
spent in a narrower range than in J2, so a threshold removes a larger
share of it. Five alternations against three change the fracture strain
by one load step (17.11 against 17.19%) and the final damage by 0.01: the
1e-3 alternation test is met within three. The CP configuration therefore
keeps W_c = 0 (§4.5).

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
`figures/pfd/cp_sweep_h_l5_{curves,panels}.png`; h20 added on 2026-10-15
at the author's request, mesh `decks/pfd_study/mesh/h20_hole5.i`, 308 hex8,
700 nodes, figure `figures/pfd/mesh_h20_hole5.png`):

| mesh | h (um) | l0/h | fracture strain | peak stress (MPa) | energy (MPa) |
|---|---|---|---|---|---|
| h20 | 5 | 1 | 16.09% | 609 | 81.0 |
| h40 | 2.5 | 2 | 17.19% | 596 | 85.0 |
| h80 | 1.25 | 4 | 17.58% | 587 | 86.3 |

h20 (l0 = h, the regularisation length equal to the cell) initiates 1.1%
strain early and 2% too stiff; h40 to h80 differ by five load steps (0.4%)
in the initiation strain and 1.4% in the peak. The sequence converges from
below in the initiation strain and from above in the peak, as the
Chakraborty 2016 Fig. 8 sequence does, and the l0 >= 2h rule is where the
change drops under the 2% the Part I grain-instance check accepted. The
production mesh stays h40 with l0 = 5 um (l0 = 2h) for both couplings.

## 5.7 CP grain-instance check (Part II analogue of Part I §7)

`decks/pfd_study/cp/grains_h40.i`: the calibrated case (Gc 1106, l0 5,
h40) on the two other Part I grain structures (instances sB and sC: new
Voronoi seedings grown the same way, new orientation draws), everything
else unchanged. Figures `figures/pfd/cp_grains_h40_{curves,panels}.png`.

| instance | fracture strain | peak stress (MPa) | energy (MPa) | max d at the end |
|---|---|---|---|---|
| A (calibration) | 17.19% | 596 | 85.0 | 0.98 |
| sB | none by 20% | 613 | 99.8 | 0.73 |
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

### 4.7 J2 length scale at fixed Gc, W_c = 20 (the carried configuration)

`decks/pfd_study/j2/sweep_l_h40.i` second pass (Gc 2003, results
`results/pfd_study/j2/sweep_l_h40_wc20/`, figures
`figures/pfd/j2_sweep_l_h40_wc20_{curves,panels}.png`):

| l0 (um) | l0/h | peak stress (MPa) | fracture strain | energy to the end (MPa) | W_c 0 fracture strain |
|---|---|---|---|---|---|
| 2.5 | 1 | 521 | none by 20% | 95.3 | none |
| 5 | 2 | 485 | 17.1% | 81.6 | 17.3% |
| 10 | 4 | 467 | 14.4% | 62.8 | 12.7% |
| 20 | 8 | 455 | 14.8% | 57.2 | 12.3% |
| 40 | 16 | 453 | 17.6% | 56.7 | 14.8% |

Same form as the W_c = 0 pass (§4.1): no initiation at l0 = h, peak
stress and energy falling and saturating with l0, the fracture strain
minimal at l0 = 4h to 8h and rising again at 16h where the band no longer
localises at the hole. The threshold raises every fracture strain beyond
l0 = 2h by 1.7 to 2.8% and the peak stresses by 70 to 100 MPa; the
saturation of the peak (455 to 453 MPa) and of the energy is unchanged.

### 4.8 J2 parameter sensitivity at W_c = 20 (the carried configuration)

`decks/pfd_study/j2/sens_h40.i` second pass (results
`results/pfd_study/j2/sens_h40_wc20/`, figures
`figures/pfd/j2_sens_h40_wc20_{curves,panels}.png`):

| case | fracture strain | peak stress (MPa) at strain | energy (MPa) |
|---|---|---|---|
| base (Gc 2003, l0 5, W_c 20, eta 0) | 17.11% | 485 at 6.1% | 81.6 |
| Gc x 0.5 (1002) | 9.22% | 460 at 4.6% | 39.7 |
| Gc x 2 (4006) | none by 20% | 527 at 10.6% | 96.2 |
| l0 2.5 (= h) | none by 20% | 521 at 10.1% | 95.3 |
| l0 10 (= 4h) | 14.38% | 467 at 4.9% | 62.8 |
| eta 1 MPa s | 17.11% | 486 at 6.1% | 81.7 |
| eta 10 MPa s | 17.50% | 489 at 6.2% | 82.5 |
| W_c 0 | 16.72% | 413 at 5.9% | 73.1 |
| W_c 5 MPa | 16.72% | 429 at 5.7% | 75.3 |
| 5 alternations | 17.11% | 485 at 6.1% | 81.6 |

The same picture as at W_c = 0 (§4.2): Gc and l0 govern the initiation
strain, the viscosity is immaterial (eta = 1 identical; eta = 10 one
load step later), the threshold moves initiation by at most 0.4% between
W_c = 0 and 20 while setting the peak stress, and five alternations equal
three. The carried J2 configuration (Gc 2003, l0 5 um, W_c 20 MPa,
eta 0, AT2, 3 alternations at 1e-3) is therefore robust in every
parameter except the two calibration targets themselves.

### 4.9 J2 mesh sensitivity at W_c = 20 (the carried configuration)

`decks/pfd_study/j2/sweep_h_l5.i` second pass (Gc 2003, l0 5; results
`results/pfd_study/j2/sweep_h_l5_wc20/`, figures
`figures/pfd/j2_sweep_h_l5_wc20_{curves,panels}.png`):

| mesh | h (um) | l0/h | fracture strain | peak stress (MPa) | energy to 20% (MPa) |
|---|---|---|---|---|---|
| h40 | 2.5 | 2 | 17.11% | 485 | 81.6 |
| h80 | 1.25 | 4 | 17.27% | 483 | 79.9 |

Two load steps in the initiation strain, 0.4% in the peak, 2% in the
energy: the carried J2 configuration is mesh-converged at l0 = 2h, as the
W_c = 0 pass was (§4.3). Together with the CP result (§5.6: five steps,
1.4%) this is the Part II answer to the l0/h question: for both couplings
the brittle rule l0 >= 2h is sufficient on this plate, with the CP plate
the more mesh-sensitive of the two.

## 7. Neural routes and surrogate on the calibrated CP+PFD plate

### 7.1 Surrogate design (data source, PINNACLE ADR 0030)

`decks/pfd_study/cp/surrogate/design_h40.i`, 24 LHS samples of (Gc 553 to
2212 MPa um, l0 4 to 10 um, rho0 24 to 96, tau_c0 2.6 to 10.6 MPa; all log
scale) on the calibrated base deck, each a coupled run with the early stop
at d = 0.95 (`scripts/surrogate_design.py design run collect`, 4 side by
side): 17 of 24 initiate (fracture strain 5.9 to 17.9%), 7 do not by 20%;
5 to 17 min per sample (median 11), 4.5 h of host time at 4 side by side.
Targets stored per sample: the engineering stress and the maximum damage
along the path (`samples/<i>/target.npz`) and the nodal damage field at
the end (`samples_damage/<i>/target.npz`). The surrogate decks
`decks/pfd_study/cp/surrogate/plate_h40_surrogate_{curve,field}.i` train on
them (held-out 25%, 1000 queries) after the hardware benchmark.

### 7.2 Neural routes inside the staggered loop

`decks/pfd_study/cp/neural_h40.i` (results `results/pfd_study/cp/neural_h40/`,
figures `figures/pfd/cp_neural_h40_{curves,panels}.png`): the calibrated
CP+PFD plate with the mechanics sub-solve of every staggered step served by
(a) the linear extrapolation predictor (PINN+FEM) and (b) the grain-aware
FE-residual-trained network as the verified solver (`corrector = auto`,
network on the GPU), against the plain FE route, all three on torch CPU
with 2 threads. What PINNACLE needed for this (ADR 0031 items 11, 14, 15):
the warm start of alternations, a predictor history that survives the
re-solved windows, and an unshared backend instance for the network's
device.

| route | fracture strain | peak stress (MPa) | fracture energy (MPa) | steps | Newton iterations, last alternation of each step | Newton time (s) | predictor time (s) | wall (s) |
|---|---|---|---|---|---|---|---|---|
| FE (cold Newton start) | 17.19% | 595.8 | 85.0 | 224 | 1008 | 147 | -- | 642 |
| PINN+FEM, linear extrapolation | 17.19% | 595.8 | 85.0 | 224 | 1008 | 145 | 0.0 | 609 |
| PINN as solver, grain-aware, verified | 17.11% | 595.6 | 83.6 | 222 | 955 | 122 | 43.5 (training 2 epochs per retrain) | 481 |

The linear extrapolation reproduces the FE run to the printed digits
(fracture strain, peak, every curve point and the Newton count of the
last alternations) with a median proposal error of 0.6% of the step
increment; its saving (5% of the wall time) comes from the first
alternation, which the hardware benchmark (§8) resolves with the per-step
counters on the same code. The verified solver accepts its proposal at
216 of 222 first alternations (the first two steps have no history to
extrapolate from and steps 5-8 were corrected by Newton after a failed
verification: 16 Newton iterations in all first alternations, median
residual of the accepted proposals 0.4% of the step) and reproduces the curve
within one load step of initiation; the staggered scheme then re-solves
the mechanics once or twice per step against the updated damage with
ordinary warm-started Newton (1068 iterations over all alternations, of
which 1052 in the re-solves), so the route's saving is the first
alternation only: wall time 481 s against 642 s (25% less, the predictor's
43.5 s included). A predictor of the re-solves (a network of the damage
increment) would be the next step and is out of this study's scope. The
solver run ended one step after initiation on a mechanics failure at
max d = 0.99 (the last row carries the cut-back window and no reaction),
where the FE and linear runs ended on the d >= 0.95 early stop.

## 8. Hardware and time benchmark of the CP+PFD routes (h40)

`decks/benchmark/hardware_pfd_h40.i`, `scripts/hardware_benchmark.py`
(`deck_file = pfd.i`), results `results/pfd_study/cp/hardware/
hardware_pfd_h40.csv`, figure `figures/pfd/hardware_pfd_h40_time_per_step.png`.
Seven cases, each the first 16 staggered steps of the calibrated CP+PFD
plate (strain per step unchanged, early stop off, max d 0.05 at the end),
one case at a time on the host, 2026-10-14 13:38-13:52. The first attempt
stopped before any run: the FE cases pointed at `cp/base_h40.i`, whose
output directory was never staged (the neural study staged its own base);
they now read the staged `neural_h40/base/pipeline.i`. The per-step
counters of ADR 0031 (item 12) separate the first alternation of every
step, which the predictor serves, from the re-solves against the updated
damage.

| case | backend / threads | wall per step (s) | Newton s per iteration | Newton iterations: first alternations / all (16 steps) | predictor s per step | peak memory (MB) |
|---|---|---|---|---|---|---|
| FE, numpy | numpy / 8 | 2.68 | 0.450 | 49 / 80 | -- | 264 |
| FE, torch CPU | torch / 2 | 1.64 | 0.150 | 49 / 80 | -- | 1,328 |
| FE, torch CPU | torch / 8 | 1.47 | 0.117 | 49 / 80 | -- | 1,328 |
| FE, torch GPU | torch GPU / 8 | 2.50 | 0.141 | 49 / 80 | -- | 1,718 |
| PINN+FEM, linear extrapolation | torch / 8 | 1.46 | 0.117 | 43 / 74 | 0.000 | 1,718 |
| PINN as solver, grain-aware, verified, network on GPU | torch / 8 | 2.75 | 0.116 | 16 / 51 | 1.22 (training 18.3 s in all, 422 assemblies) | 2,124 |
| same, network on CPU | torch / 8 | 2.60 | 0.116 | 16 / 51 | 1.17 (17.7 s, 396 assemblies) | 2,124 |

Findings. (i) The backend ranking of Part I holds for the coupled problem:
torch CPU is 3.0x (2 threads) to 3.9x (8 threads) faster per Newton
iteration than numpy, the GPU no faster than the CPU at this size (860
elements) and 2.5 s per step against 1.5 s because the damage solve and
the staggered bookkeeping run on the CPU either way; the production arm
(torch CPU, 2 threads, seven runs side by side) stays the right choice.
(ii) The coupled step costs 5 Newton iterations in all (3 in the first
alternation, 2 in the re-solves, 1.9 alternations per step at this early
stage); the mechanics Newton is 17% of the wall time of a step on torch (37%
on numpy); the rest is the damage assembly and solve (sparse, 1,864
nodes), the CP state handling and the reaction and field evaluations of
the coupled loop, which the backend does not speed up. (iii) The linear extrapolation saves 6 of the 49
first-alternation iterations (12%) at no cost, less than Part I's 20%
because the coupled step is already warm-started by the re-solves. (iv)
The verified solver removes the first-alternation Newton on 10 of the
16 steps (the 16 iterations are the two cold-start steps without history
and four steps whose proposal failed the 5% residual test and were
corrected by Newton; the same six steps as in the full run, after which
every proposal is accepted) but trains for 1.2 s per step, so its step is 1.8x
the FE step at this size; in the full run (§7.2) the retraining is
cheaper per step (43.5 s over 222 steps, 0.2 s) and the route ends 25%
faster than the FE. The network's device does not matter (the cost is the
numpy assemblies of the FE residual, as in Part I).

### 7.4 Plate surrogates (PINNACLE ADR 0030, data source) on the CP+PFD plate

Decks `decks/pfd_study/cp/surrogate/plate_h40_surrogate_{curve,field}.i`
(data source over the 24 samples of §7.1, 18/6 split, seed 0, 3x64 tanh,
Adam 4000 + L-BFGS 50, CPU, `evaluate` = 1,000 random rows of the design
box), comparison `decks/pfd_study/cp/surrogate/compare_h40.i` by
`scripts/surrogate_compare.py` (now accepts a study summary as the routes'
cost and draws held-out field maps when `mesh_vtk` is given); results
`results/pfd_study/cp/surrogate/{h40_curve,h40_field,cp_surrogate_h40_summary.csv}`,
figures `figures/pfd/cp_surrogate_h40_{curve_holdout,field_holdout_fields,cost_per_query}.png`.
The curve deck first refused to run: with two declared variables (u, d)
`parametric_surrogate` needs `variable` named even for a reduced target
(it only sets the result layout); `variable = u` added. The sample cost
in `samples.csv` is now the wall time of the coupled run (the collector
summed the Newton time, 17-37% of a coupled step, §8).

| surrogate | target | held-out rel. L2 mean / max | training (s) | query (ms) | data cost (s, 18 training samples) | break-even vs FE / linear / verified solver (runs) |
|---|---|---|---|---|---|---|
| curve | eng. stress and max damage vs nominal strain over (Gc, l0, rho0, tau_c0) | 14.7% / 69.8% (one sample; 1.6-4.9% on the other five) | 10.3 | 0.1 | 12,987 | 20 / 21 / 27 |
| field | nodal damage at the end of the path | 5.8% / 7.1% | 198 | 0.4 | 12,987 | 21 / 22 / 27 |

Held-out samples (1, 7, 16, 17, 21, 23), curve surrogate, per column
(stress rel. L2, max-damage rel. L2, fracture strain FE -> surrogate at
the first crossing of d = 0.8):

| sample | Gc | l0 | end of path | stress | max d | fracture strain FE -> surrogate |
|---|---|---|---|---|---|---|
| 1 | 702 | 9.0 | 6.0% | 69.8% | 72.3% | 5.86% -> 5.31% |
| 7 | 828 | 8.6 | 8.7% | 4.9% | 3.1% | 8.52% -> 8.44% |
| 16 | 1205 | 8.1 | 12.9% | 4.9% | 4.8% | 12.66% -> 12.89% |
| 17 | 872 | 8.0 | 10.5% | 2.5% | 2.3% | 10.31% -> 10.31% |
| 21 | 947 | 6.0 | 10.9% | 4.7% | 6.6% | 10.78% -> none (predicted max d below 0.8; the stress drop is predicted at 10.9%) |
| 23 | 1341 | 7.3 | 15.2% | 1.6% | 0.9% | 14.92% -> 14.84% |

Findings. (i) The curve surrogate predicts the fracture strain within
0.25% strain on four of the six held-out samples and the stress drop on
five; its mean error is set by sample 1, the lowest Gc and highest l0 of
the held-out set, where the drop comes 0.55% strain early and the whole
rel. L2 of a curve that ends at 6% is paid by the shifted collapse
(peak 478 against 474 MPa). As in Part I, accuracy is data-limited at
the box's edge: the training set has three samples at l0 > 8 um. (ii)
The damage-field surrogate is uniformly accurate (2.8-7.1%); the
difference maps put the error in the band below the hole where the
crack grows, with the surrogate over-predicting the band's damage on four
samples and under-predicting it on the earliest-fracturing one; it does
not learn the sample-to-sample state at the end of the path (some samples
end at initiation, some at 20% without it), so it is a map of "damage at
the end of the path", not at a fixed strain. (iii) Cost: data 3.6 h of
host time (four samples side by side at two threads) against 10 s and
198 s of training; a query costs 0.1-0.4 ms; the surrogate pays for
itself after about 20 FE runs (27 verified-solver runs), the same
ballpark as Part I (14). The query for the fracture strain of a new
(Gc, l0, rho0, tau_c0) is the use the author asked for in the plan
(docs/09 §5.1 step 8); its accuracy away from the box's edge is a few
tenths of a percent strain.

## 9. Stress-concentrator shape study (CP+PFD, h40)

`decks/pfd_study/cp/shape_h40.i` (meshes `decks/pfd_study/mesh/h40_ell_*.i`,
figures `figures/pfd/mesh_h40_ell_*.png`; the author's go-ahead of
2026-10-14 for all three series of docs/09 §6), the calibrated CP+PFD set
(Gc 1106, l0 5, W_c 0) on the same 30-grain structure, 3 side by side at
2 threads, 2026-10-15 00:18-00:50; figures
`figures/pfd/cp_shape_h40_{curves,panels,damage,equivalent_plastic_strain,crack}.png`
and `cp_shape_h40_evolution.gif`; summary
`results/pfd_study/cp/shape_h40/shape_h40_summary.csv`. All six ellipse
meshes keep the h40 background and the O-grid ring (24 segments; min scaled
Jacobian 0.41-0.71); the sharpest root (5 x 1.25, root radius 0.16 um) is
under-resolved by the 0.45 um hole arc, which the plan foresaw.

| series | case | K_t (Inglis) | fracture strain | peak stress (MPa) | energy (MPa) | cracked elements at the end |
|---|---|---|---|---|---|---|
| -- | circle 5 um | 3 | 17.19% | 596 | 85.0 | 40 |
| A (a = 2.5 um) | 5 x 2.5 | 5 | 17.81% | 602 | 89.1 | 61 |
| A | 5 x 1.25 | 9 | 18.36% | 606 | 92.2 | 7 |
| B (area of the circle) | 7.07 x 3.54 | 5 | 16.33% | 628 (a last-step spike; 590 before it) | 79.0 | 45 |
| B | 10 x 2.5 | 9 | 15.47% | 572 | 74.1 | 11 |
| C (5 x 2.5 rotated) | 45 deg | -- | 18.52% | 606 | 93.0 | 82 |
| C | 90 deg | 2 | 18.91% | 610 | 95.2 | 56 |

Findings. (i) Series A contradicts the elastic expectation: at fixed
transverse extent a sharper root initiates later (17.2 -> 17.8 -> 18.4%),
not earlier. In the CP plate the driving force is the plastic work, which
the sharp root concentrates into a smaller volume: the damage maps show
the 5 x 1.25 case's hot spot confined to the root elements (7 cracked
elements at 18.4%) while the circle's spreads over the grain below the
hole. The under-resolved root (one 0.45 um element on a 0.16 um radius)
also smooths the concentration; the case is reported with that caveat.
(ii) Series B follows the Irwin scaling: a larger transverse extent (r_p
proportional to a) initiates earlier (16.3% at 7.07 um, 15.5% at 10 um)
with a lower energy, and the crack leaves the root along the transverse
direction. (iii) Series C: rotating the 5 x 2.5 ellipse delays initiation
(17.8% at 0 deg, 18.5% at 45 deg, 18.9% at 90 deg) as K_t falls from 5 to
2, and the 45 deg case tilts the crack: it leaves the root along the
grain boundary below the hole and grows into the grain on the lower
right, where the circle's crack also turns (the same grain-boundary path
that the fields of §6 show). (iv) Every case initiates at the hole's root
or the grain boundary adjacent to it, and the peak stress varies by only
+/-2% because the gross section is 95% intact until initiation. The
sharpest-root case would need the h80 ring (0.25 um arc) to be
quantitative.
