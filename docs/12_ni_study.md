# Ni single-crystal benchmark (Part III): lab notebook

Running record of the Part III work (plan: `docs/11_ni_plan.md`). Every
number traces to a deck under `decks/ni/` and a row of `docs/run_ledger.md`.

## 1. Pipeline check and the CP integration gap (2026-10-15)

`decks/ni/cp/check_L1_h040.i` (`scripts/ni_specimen.py`): specimen L.1 on
the h040 dog-bone (220 hex8), the paper's Table 3 parameter set, grip
strain 2% in 16 steps, local numpy. The pipeline runs end to end (mesh ->
deck -> run -> gauge strain and engineering stress), but the curve is far
above the experiment: 228 MPa at 0.18% gauge strain and 507 MPa at 3%,
against about 85 and 92 MPa in Fig. 3.

Diagnosis, in order:

| test | result |
|---|---|
| material point (`scripts/verify_material_point.py` helpers), Table 3 set, L orientation, isochoric uniaxial drive at 1e-3/s | sigma_11 62 MPa at 0.1%, 82 MPa at 3%; tau_c(rho_0) = 16.4 MPa: the constitutive law with the paper's constants gives the paper's stress level |
| the same in um units (b 2.5e-4 um, rho_0 2.35 um^-2, y_c 7.9e-4 um) | identical: the law is unit-invariant |
| the FE check with gamma_0 x 100 | stress changes by 3%: the FE overshoot is not the rate term |
| the FE check with the mesh scaled to um | identical to the mm run |
| the FE check with 16 / 64 / 128 / 256 load steps to 2% grip strain (gauge increment 2.3e-3 / 5.7e-4 / 2.9e-4 / 1.5e-4) | stress at 3% gauge strain 507 / 180 / 133 / 110 MPa: first-order convergence in the step size, limit about 87 MPa |

Cause (PINNACLE, `cp_field_callback.py`): during the Newton solve the
stress is evaluated with the plastic state frozen at the last accepted
step (`compute(..., dt=None)`), and the state advances only after
convergence (`update_state`, sub-stepped). The step's plastic slip never
enters the equilibrium stress of that step, so the flow stress overshoots
by about E x (strain increment); the code warns about it above 8e-4 and the
error is first order below. Parts I and II ran at 7.8e-4 nominal increment
and absorbed the overshoot in the calibrated Al2219 set (the "framework
sensitive" parameters of Part I); the Ni study keeps the paper's constants,
so the integration must be converged at a practical step count (a 50%
path at 1.5e-5 increment would need 30,000 steps).

Reference for the check: the material point under TRUE uniaxial stress
(the lateral stretches iterated to zero lateral stress at every step of a
3,000-step explicit path, `scripts/verify_material_point.py` helpers), L
orientation, paper's set: 92 MPa at 0.1% and 120 MPa at 3% at 1e-3/s,
126 MPa at 3% at 1.5e-3/s (the gage's rate when the grips move at 1e-3/s).
So the paper's constants in PINNACLE's law give a flow stress about 30%
above the Fig. 3 experiment (92 MPa at 3%); the strain rate of the test is
not stated (a tenfold lower rate would lower the stress by 24% with
m = 0.12), and the calibration of (rho_0, k_m, y_c, tau_c0) absorbs the
rest. The FE must reproduce this 126 MPa before any calibration.

### 1.1 Integrator study on the dog-bone (PINNACLE ADR 0033)

Three state-update schemes of `stress_divergence_cp` (PINNACLE ADR 0033),
the same h040 dog-bone and parameter set, grip strain 2% at 1e-3/s, numpy,
3-4 threads; "stress" is the engineering stress from the reaction on the
x+ grip (summed after the state advance in every scheme), compared at 3%
gauge strain with the 126 MPa reference:

| scheme | load steps (gauge increment) | stress at 3% gauge (MPa) | error | wall |
|---|---|---|---|---|
| explicit (Parts I-II) | 16 (2.3e-3) | 507 | +300% | 1 min |
| explicit | 64 (5.7e-4) | 189 | +50% | 1 min |
| explicit | 128 (2.9e-4) | 133 | +6% | 2 min |
| explicit | 256 (1.5e-4) | 111 | -12% | 4 min |
| explicit | 1024 (3.7e-5) | 121 | -4% | 17 min |
| iterated, 2 passes, cap 1e-4 | 16 | 327 | +160% | 3 min |
| iterated, 2 passes, cap 1e-4 | 64 | 178 | +41% | 5 min |
| iterated, 2 passes, cap 1e-5 | 256 (1.5e-4) | 136 | +8% | 15 min |
| iterated, 3 passes, cap 1e-5 | 64 | 144 | +14% | 10 min |
| implicit, sub-stepping integrator, cap 1e-5 (single element, 16 steps of 2.5e-4) | -- | within 0.5% of a 1,024-step explicit path | -- | 12 s; 25 s per step on the dog-bone |
| implicit, backward-Euler integrator (single element, 4 and 16 steps) | 4 (1e-3) and 16 (2.5e-4) | 91.1 against 90.3 of the 1,024-step explicit path | +1% | 5 Newton iterations per step |
| implicit, backward-Euler integrator | 64 (5.7e-4) | 120.7 at 3.0%; 142.3 at 4.75% (gage field 142.6) against about 140 | -4% at 3%, +2% at 4.75% | 20 s per step, 5.7 Newton iterations per step |

The explicit scheme is not monotonic in the step count because its two
errors have opposite signs (the trial stress overshoots, the recorded
relaxed-state reaction undershoots) and cross between 128 and 256 steps;
it needs increments near 1e-4 for 5% and is 4% low at 3.7e-5. The
iterated scheme (the step re-solved with the state integrated along the
converged path frozen) is not better than the explicit one at the same
step count: its last pass still leaves the step's flow out of the
equilibrium. The implicit scheme integrates the state inside every
residual evaluation; with the path-following explicit sub-stepping it is
accurate but costs an integration per evaluation (25 s per step here),
with the backward-Euler slip update (a local Newton on the slip
increments of every point, the material's exact evaluation at the trial
end state, a semi-analytic Jacobian, backtracking, a Jacobi phase for
trials far off equilibrium, backward-Euler sub-increments and the
sub-stepping as fallbacks) it matches the material point within 2% at
increments up to 1e-3 and converges the global Newton in 5-6 iterations
with the full-gradient finite-difference tangent. Part III runs the
implicit backward-Euler update (`state_update = implicit`, the pipeline's
default). What the checks also showed: the reaction is to be recorded
after the state advance in every scheme (with the explicit update the
relaxed state's load is the better first-order estimate of the two, with
the implicit one the committed state is what the step solved; the Part II
CP curves are therefore the explicit scheme's post-advance reaction, which
sits about 7% below its own trial equilibrium at 7.8e-4 increments, an
error the Part I calibration absorbed), and the J2 route is unaffected
(its committed state does not change the stress at fixed strain).

Decision: PINNACLE gets a consistent stress update (ADR 0033): with
`state_update = implicit` the callback's `stress()` integrates the state
from the committed one over the step's time increment (the existing
sub-stepper) at every residual evaluation, and the tangent is the
central-difference derivative of that stress per quadrature point; the
default stays explicit so Parts I-II reproduce. Validation: the implicit
FE at 64 steps must reproduce the explicit limit (87 MPa at 3%) and a
single element under uniaxial stress must match the material point.

## 2. CP calibration on specimen L.1 (docs/11 §5 stage 2)

Deck-driven Gauss-Newton as in Part I: `decks/ni/cp/base_L1_h040.i` (the
paper's Table 3 set, h040, 360 steps of 1e-3 grip strain, implicit update)
is the center, `decks/ni/cp/cal_gnN.i` runs it with one multiplicative
perturbation per free parameter (x1.2; tau_c0 x1.5) through
`scripts/ni_study.py`, and `scripts/ni_gauss_newton.py` takes the damped
step in log space against the Fig. 3 experiment markers
(`data/digitized/ni_fig3_experiment_markers.csv`, engineering stress
against GAUGE strain, markers below 0.3% strain excluded) and writes the next
center and GN deck. Runs on the host, torch CPU, 2 threads, 5 side by side,
about 90 min per iteration (10-15 s per step; the first step 45 s).

### 2.1 Iteration 1: the paper's set at the test's strain rate (2026-10-15 19:40-21:15)

The author set the bar strain rate to 5e-2/s (t_end = strain / 0.05; an
earlier 1e-3/s iteration, `results/ni_study/cal/gn1_rate1e-3`, gave RMSE
182.6 MPa and is kept as the record of the paper's set at that rate).

| gauge strain | experiment (MPa) | FE, paper's set, 5e-2/s (MPa) |
|---|---|---|
| 2% | 87.7 | 169 |
| 5% | 98.4 | 225 |
| 10% | 129.3 | 332 |
| 20% | 180.1 | 521 |
| 30% | 211.8 | 652 |
| 40% | 232.8 | 758 |
| 50% | 238.4 | 841 |

RMSE 373.1 MPa: the paper's set in PINNACLE's law gives 942 MPa at 66%
gauge strain against the experiment's 240 at 51%, i.e. it hardens about
four times too fast at this rate (the rate term alone adds 60% between
1e-3 and 5e-2/s with m = 0.12). Sensitivities dS/dln p at 2 / 10 / 30%
(MPa): rho_0 70 / 58 / 50, k_m 29 / 226 / 562, y_c -3 / -55 / -181,
tau_c0 3 / 3 / 3: at this rate the Peierls term is immaterial, the yield
is set by rho_0 and the rate term, the hardening by k_m and y_c. The step
(clipped to +-0.7 in ln p) goes to rho_0 1.17e6, k_m 0.0146, y_c 0.74e-6,
tau_c0 0.255 (linear prediction 25.3 MPa RMSE); iteration 2 runs from
that center, and the clipping means two or three more iterations.

### 2.2 Iteration 2 (2026-10-15 21:15 - 2026-10-16 00:00)

Center rho_0 1.17e6, k_m 0.0146, y_c 0.74e-6, tau_c0 0.255:

| gauge strain | experiment (MPa) | FE (MPa) |
|---|---|---|
| 2% | 87.7 | 117 |
| 5% | 98.4 | 148 |
| 10% | 129.3 | 203 |
| 20% | 180.1 | 301 |
| 30% | 211.8 | 373 |
| 40% | 232.8 | 432 |
| 50% | 238.4 | 480 |

RMSE 141.2 MPa (from 373): the step halved the hardening but the curve is
still 1.9x too high at 30%. Sensitivities at 2 / 10 / 30%: rho_0 50 / 50
/ 40, k_m 14 / 102 / 281, y_c -2 / -25 / -90, tau_c0 3 / 3 / 3. Next
center (clipped step) rho_0 7.64e5, k_m 0.0108, y_c 1.48e-6, tau_c0 0.249,
linear prediction 4.6 MPa; iteration 3 runs from it.

### 2.3 Iteration 3 (2026-10-16 00:00-01:50)

Center rho_0 7.64e5, k_m 0.0108, y_c 1.48e-6, tau_c0 0.249:

| gauge strain | experiment (MPa) | FE (MPa) | error |
|---|---|---|---|
| 2% | 87.7 | 94 | +7% |
| 5% | 98.4 | 114 | +16% |
| 10% | 129.3 | 147 | +14% |
| 20% | 180.1 | 201 | +12% |
| 30% | 211.8 | 241 | +14% |
| 40% | 232.8 | 271 | +17% |
| 50% | 238.4 | 295 | +24% |

RMSE 29.1 MPa (from 141). Sensitivities at 2 / 10 / 30%: rho_0 41 / 42 /
31, k_m 9 / 56 / 156, y_c -3 / -23 / -79, tau_c0 3 / 3 / 3. Next center
rho_0 6.86e5, k_m 0.0111, y_c 2.26e-6, tau_c0 0.246 (linear prediction
4.5 MPa, the markers' own spread); iteration 4 decides.

### 2.4 Iteration 4: converged (2026-10-16 01:50-08:15)

Center rho_0 6.86e5, k_m 0.0111, y_c 2.26e-6, tau_c0 0.246:

| gauge strain | experiment (MPa) | FE (MPa) | error |
|---|---|---|---|
| 0.5% | 85.0 | 78 | -8% |
| 2% | 87.7 | 88 | +1% |
| 5% | 98.4 | 107 | +8% |
| 10% | 129.3 | 134 | +4% |
| 20% | 180.1 | 177 | -2% |
| 30% | 211.8 | 207 | -2% |
| 40% | 232.8 | 229 | -1% |
| 50% | 238.4 | 246 | +3% |

RMSE 4.6 MPa (from 29.1), the markers' own spread; the next linear step
predicts 4.5 MPa, so the calibration is converged and the iteration-4 center
is the calibrated set, `decks/ni/cp/calibrated_L1_h040.i` (the further step the
script wrote, cal_center5/cal_gn5, was not run and is deleted). Sensitivities at 2 / 10 / 30%: rho_0 38 / 40 /
28, k_m 9 / 46 / 126, y_c -4 / -25 / -78, tau_c0 3 / 3 / 3 MPa per unit ln p:
tau_c0 is barely identifiable from this curve (it moves the whole curve by
3 MPa per e-fold), the hardening pair k_m / y_c carries the shape.

Wall time was 6.3-6.5 h per run (1.5-2.6 h in iterations 1-3): the five runs
shared the host with the three Al2219 implicit Gc trials and run to 36% bar
strain, where the late steps are slow (below).

**Late-strain solver failure.** In all five runs the global Newton stops
converging beyond 54-60% gauge strain (first cutback at steps 273-306 of
360); the executioner then takes steps after exhausting its 10 cutbacks and
the engineering stress jumps from 252 to 350-520 MPa. Those steps are not a
solution and are not used: every marker is at or below 50% gauge strain and
the test fractures at 46.3% (Fig. 4b). The CP+PFD runs stop at fracture
before this range; any CP-only curve beyond 54% is reported truncated there.

### 2.5 Calibration figures (for the author's review)

`python3 scripts/ni_calibration_figures.py decks/ni/figures/cal_figures.i`
writes the Al2219 calibration set for the calibrated iteration (`figures/ni/cal/cal_curve`,
`cal_oat_panels`, `cal_oat_overlay`, `cal_error`, `cal_sensitivity`) and `cal_convergence.png`
over all iterations (author, 2026-10-16: earlier iterations appear only in the convergence
figure). Curves end at the first load cutback (dotted line). Report: Part III report §4.

## 3. Mesh study on L.1 with the h040-calibrated set (2026-10-16/17)

`decks/ni/cp/mesh_L1.i` (`ni_study.py`): the h040-calibrated set on h020 (gage cell
0.2 mm, 792 hex8), h010 (0.1 mm, 3,168) and h020 with two elements through the
thickness (1,584); grip strain 30% in 300 steps; h040 is the calibration run. The
host froze at about 08:40 on 2026-10-17 (no memory cap on the Ni pipeline yet: the
h020z2 run fell back to sub-stepping at step 229 after ten cutbacks, the same
runaway as docs/10 §10, now seen in a CP-only run); h020 had finished, h010 reached
step 112 (22.7% gauge strain), h020z2 step 248 (clean to 43.5%). Engineering
stress (MPa) at the marker strains, curves to the first cutback:

| case | 2% | 10% | 20% | 30% | 40% | 50% |
|---|---|---|---|---|---|---|
| experiment | 87.7 | 129.3 | 180.1 | 211.8 | 232.8 | 238.4 |
| h040 | 88.5 | 134.0 | 177.0 | 207.3 | 229.3 | 246.1 |
| h020 | 83.4 | 120.5 | 152.9 | 173.7 | 187.2 | 195.8 |
| h010 (partial) | 81.9 | 116.1 | 146.3 | - | - | - |
| h020, two layers | 83.3 | 119.1 | 151.9 | 172.8 | 186.3 | - |

Refinement lowers the curve by 14% (h040 to h020) and a further 4% (h020 to h010)
at 20% gauge strain: converging (ratio about 3.5), limit about 24% below h040; a
second layer through the thickness changes under 1%. The equivalent plastic
strain fields (`figures/ni/mesh/mesh_eqps.png`, `scripts/ni_fields.py
decks/ni/figures/fields_mesh.i`) are nearly the same on h040 and h020 at equal bar
strain (peaks 0.24/0.25, 0.42/0.41, 0.60/0.61); the gage bends into an S as the
crystal shears with clamped grips, and the coarse mesh carries more force for the
same deformation. The h040 calibration therefore absorbed a mesh error of about
20%. Figures: `scripts/ni_mesh_figures.py decks/ni/figures/mesh_figures.i`.

**Author's decision (2026-10-17):** recalibrate on h020; the report shows the h020
calibration only, and the mesh sensitivity keeps h040 as a reference against the
converging finer meshes. No study of the source of the coarse-mesh stiffness: the
gage carries most of the plastic strain and the rate is high, so the denser mesh is
needed in any case.

## 4. CP calibration on h020 (2026-10-17)

`decks/ni/cp/h020/cal_gn1.i`: center = the h040-calibrated set on h020 (the mesh
study's h020 run, copied in as `base`), four perturbations at 2 threads each,
memory cap 16 GB per run; `ni_gauss_newton.py` now writes the next decks next to the
given one (`decks/ni/cp/h020/`, results `results/ni_study/cal_h020/`). About 4 h per
iteration.

### 4.1 Iteration 1 on h020 and a third PINNACLE bug (2026-10-17)

Center (the h040 set on h020): RMSE 29.7 MPa, the curve 12-18% below the markers
from 10% strain. The k_m x1.2 run cut back twice at 14.6% gauge strain and its
curve then sat 46 MPa above the center at 20% (h040's sensitivity predicts about
+19): `cp_quasi_static` kept the cutback's halved time increment for the full load
increments that followed, i.e. twice the strain rate for the rest of the run per
cutback level (+9% flow stress per halving at n = 8.3; 199/172 = 1.16 against
2^(2/8.33) = 1.18). Fixed in PINNACLE ddc76c3 (the time increment now follows the
applied load increment; regression test; ADR 0031 addendum). Earlier results are
unaffected: Al2219 Part I had no cutback, Part II one after its fracture, and the
Ni h040 calibration and mesh study cut back only beyond the curves' cut-off. The
step computed from the corrupted run (k_m sensitivity 310 MPa at 30%) was discarded;
the k_m case is rerun with the fix (the corrupted run kept as
`results/ni_study/cal_h020/gn1_k1_p_dtbug` on the host).

The k_m case rerun with the fix (15:23-20:13) cut back once at the same point; its
time increment returned to 0.02 s after the two half steps and the curve stayed on
trend (one high point, 111.6 against about 98.7 N of reaction, at the first
accepted half step, then back on trend; to be understood before curves with cutbacks
are relied on). Iteration 1 step (center RMSE 29.7 MPa): sensitivities at 2 / 10 / 30%
rho_0 38 / 37 / 25, k_m 9 / 47 / 126, y_c -3 / -26 / -78, tau_c0 3 / 3 / 3 (k_m as on
h040); next center rho_0 7.26e5, k_m 0.01147, y_c 1.54e-6, tau_c0 0.244 (predicted
3.6 MPa). Iteration 2 launched 20:20 (`decks/ni/cp/h020/cal_gn2.i`, 5 x 2 threads).

### 4.2 Iteration 2 on h020: converged (2026-10-17 20:20 - 2026-10-18 04:35)

Center rho_0 7.256e5, k_m 0.01147, y_c 1.542e-6, tau_c0 0.244: RMSE 4.1 MPa (from
29.7), the next step predicting 3.4 MPa; every marker from 2% to 45% within -1% to
+6%, 0.5% at -7%, 50% at +4% (the experiment's softening before failure). This is the
calibrated set on h020 (`decks/ni/cp/h020/calibrated_L1_h020.i`); the further step is
not taken. Sensitivities at 2 / 10 / 30%: rho_0 39 / 39 / 27, k_m 10 / 59 / 161, y_c
-3 / -26 / -83, tau_c0 3 / 3 / 3. The center and tau_c0 runs cut back 5 and 4 times
from 41% gauge strain and recovered every time (all steps converged, valid since
ddc76c3); the k_m run exhausted its cutbacks at 43% (the step's k_m column at 45 and
50% is therefore clamped, which only affects the step not taken).

The y_c run shows two one-step spikes (+9% at 15.9% and 56.8%, no cutback): at each,
the backward-Euler update fell back to sub-stepping inside a line-search trial, the
16 GB cap stopped the memory runaway within 6 s (traceback: `_integrate_stacks` ->
`compute_stacked`, a (2,861 x 12) allocation), and `cp_quasi_static` skipped the
step with the material state left behind, so the next step integrated two increments
over one time increment. Fixed in PINNACLE 9f31bb0 (a raised step is cut back and
retried; test; ADR 0031 addendum). A harness driving the callback's fallback directly
(6,336 points, 45% far off, with the real backward-Euler cascade) keeps resident and
virtual memory flat, so the runaway is still not located; the cap contains it. The
y_c and k_m runs are rerun with 9f31bb0 for clean figures (the old runs kept as
`gn2_*_pre9f31bb0` on the host).

## 5. Mesh study with the h020-calibrated set (2026-10-18, running)

`decks/ni/cp/h020/mesh_L1.i`: h040 (reference), h010 and h020z2 with the calibrated
set; h020 is the calibrated run (copied as `base`). Threads 2 / 8 / 4, memory cap 16 /
24 / 16 GB.

### 5.1 The address-space cap was the wrong measure (2026-10-18)

h020z2 stopped at step 36 (11% gauge strain): the 16 GB cap failed a 58 MB request in
the ordinary tangent assembly (`assembler._tangent_analytic`), and every retry after
it (275 raises, 10 cutbacks) failed the same way; the run is void. The k_m rerun of
§4.2 raised at every step after 231 (load factor 0.77) for the same reason, at 5.8 GB
resident, and was stopped. h010 held 12.6 GB resident against 24.5 GB of address space.

On the host the address space runs about 11 GB above resident memory, a fixed amount:
the ROCm torch libraries map about 7 GB of shared objects and reserve a 4 GB block,
and touch almost none of it (`pmap`). The 16 GB address-space cap therefore left
about 5 GB for the run itself. The local torch build maps far less (address space
2.9 GB above resident).

No leak was found. Locally the calibrated h020 deck rises from 1.5 to 2.5 GB over the
first 13 steps and holds at about 2.4 GB through step 21; under tracemalloc the h040
deck holds 1.48-1.52 GB resident and 130-150 MB of traced Python/numpy memory over 132
steps (to load factor 0.44; no fallback or cutback in that range),
the only churn being the temporaries of the finite-difference tangent and the
backward-Euler update, freed every iteration. h010's resident memory goes both ways
(12.6 GB at step 96, 11.0 GB at step 100): a large working set (25,344 points), not
an accumulation. The Al2219 runaway (1.4 GB to 52 GB within 12 min, docs/10 §10) was
real resident growth; the Ni cap failures were not evidence of it.

Fix (study scripts): `ni_specimen.py` and `pfd_plate.py` now cap resident memory
through a user cgroup scope (`systemd-run --user --scope -pMemoryMax=<cap>G
-pMemorySwapMax=0`, the kernel kills the run at the cap) and fall back to the
address-space limit only where systemd-run is unavailable. The k_m case and h020z2
were rerun under the new cap (2026-10-18 13:00); the failed runs are kept in
`failed_ascap/` of each case directory. h010 keeps running under its 24 GB
address-space cap and is rerun under the new cap if it hits it.

The k_m rerun under the cgroup cap ended the same way the first k_m run did: past 41%
gauge strain every step needed one cutback, the run-total budget (`max_cutbacks` 10)
was spent at step 231 (load factor 0.770, 43.0% gauge strain), and steps 232-234 did
not converge (the executioner moves on after the budget is spent; results.csv keeps
converged steps only). At step 234 the backward-Euler update of one point fell back to
sub-stepping (17:13:35) and the process reached the 16 GB cap 17 s later; the kernel
killed it (systemd: `oom-kill`, memory peak 16.0 G). So the fallback runaway is real
in CP-only runs too, and only after the run has left equilibrium; the cgroup cap
contained it without touching the other runs. Its curve to the last converged step
(241 steps, 43.0%) is the k_m column of iteration 2; GN summary unchanged (RMSE 4.1
MPa, same sensitivities). The author chose to keep this curve and state the early end
rather than rerun with a larger budget (2026-10-18). The figures (`figures/ni/cal/`)
mark the end on the one-at-a-time panels and overlay and draw the k_m sensitivity
only to the 40% marker.

h020z2 (rerun under the cgroup cap) spent its ten cutbacks by step 226 (load factor
0.753, 44.0% gauge strain; 238 converged steps), then fell back to sub-stepping and
was oom-killed at 16 GB (21:16). h010 spent its tenth cutback at step 199 (37.8%) and
converged step 200 at the full increment; it continues, and its next failure ends its
curve. The two-layer and fine meshes therefore cover the markers to 40% (h010 at
least to 38%), h040 and h020 to 50%.

`ni_specimen.py curve` kept every converged row; after the budget is spent the
executioner records the failed step as unconverged and moves on from a state that is
not a solution, so a step converging after it would have entered the curve. The curve
now stops at the first unconverged row (a recovered cutback leaves no row; the
iteration-2 curves are unchanged, 305 / 241 / 300 / 304 / 302 steps).

### 5.2 Mesh study results with the h020 set (2026-10-19)

h010 spent its budget at step 243 (load factor 0.810); its curve ends at step 242,
47.1% gauge strain (253 converged steps, 23 h). Engineering stress (MPa) and the
change against h020 (`figures/ni/mesh_h020/`):

| gauge strain | h040 | h020 | h010 | h020 two layers | h040/h020 | h010/h020 | two layers/h020 |
|---|---|---|---|---|---|---|---|
| 2% | 92.0 | 86.9 | 85.3 | 86.7 | +5.9% | -1.8% | -0.2% |
| 10% | 147.2 | 134.0 | 129.7 | 132.3 | +9.9% | -3.2% | -1.3% |
| 20% | 203.7 | 179.4 | 172.2 | 177.9 | +13.5% | -4.0% | -0.8% |
| 30% | 244.7 | 210.9 | 199.9 | 209.1 | +16.0% | -5.2% | -0.9% |
| 40% | 275.8 | 232.8 | 217.5 | 230.7 | +18.5% | -6.6% | -0.9% |
| 45% | 288.5 | 240.9 | 223.6 | - | +19.7% | -7.2% | - |

The difference shrinks by a factor 0.30-0.37 per halving of h (observed order 1.5-1.8).
Richardson extrapolation from h040/h020/h010 puts h020 +2.7 / +4.9 / +6.0 / +8.4 /
+11.4% above the h -> 0 limit at 2 / 10 / 20 / 30 / 40% (h010 +0.9 to +4.0%). The
mesh is converging, but h020 is within 5% of the limit only to about 10% strain; the
calibrated set absorbs the rest (it is an h020 set). Two layers through the
thickness change the curve by less than 1.3%. Awaiting the author's reading.
h010 then ran away at the address-space cap (56 raised steps, all after its curve
end) and finished at 05:35 with its snapshots (valid to 4.8 s, step 240).

Fields (`figures/ni/mesh_h020/mesh_eqps.png`, bar strain 8 / 16 / 24%): all three meshes
give the same gauge strain for the same bar strain (15.7-16.2 / 29.8-30.1 / 45.1-46.7%)
and the same S-bend of the gage by 24%. h010 resolves diagonal bands of concentrated
slip at the gage-fillet junctions (x about 9-10 and 19-20 mm) that h040 smears and h020
only hints at; its peak equivalent plastic strain is higher (0.25 / 0.36 / 0.51 against
0.20-0.21 / 0.32 / 0.48). The stiffness of the coarse meshes is consistent with these
unresolved bands (an observation, not a measured cause).

## 6. CP+PFD on L.1 h020 (2026-10-19)

Author (2026-10-19): h020 and the h020 set are final; next the G_c / l0 calibration of CP+PFD on
L.1 h020 as for Al2219, without an h-vs-l0 study; then the 15 specimens with FEM, PINN+FEM, PINN
and the surrogate, compared with the experiments (curve, failure strain, fracture shape, axial
strain field against the DIC panels of Fig. 5).

Pipeline: `ni_specimen.py` renders `decks/ni/templates/cp_pfd_specimen_template.i` when the
specimen deck sets Gc (AT2, zeta 1, W_c 0, eta 0, three alternations, tolerance 1e-3, damage
fixed at 0 on the grip faces); `fracture` applies the paper's criterion as a separation test
on the damage snapshots (d > 0.85 nodes leave no intact path from x- to x+), checked on
synthetic bands (straight cut: location 0.25, 90 deg; 45-deg cut: 41 deg; half-width cut: not
separated). Smoke test (local, h040, 4 steps to 1.1% gauge strain, Gc 125, l0 0.8): two
alternations per step, max d 0.012 at 1.1%.

Check run `decks/ni/pfd/check_L1_h020.i`: the paper's G_c 125 MPa mm, l0 0.4 mm (2h), to see the
cost and where L.1 fails before the calibration design is fixed.

Author (2026-10-19): G_c is calibrated to 0.463 at l0 = 0.4, 0.6 and 0.8 mm (2h to 4h), and l0 is
chosen by L.1's fracture shape in Fig. 5 (inclined about 45 deg at about 0.20 of the gage, near
its thinnest station). The check run takes about 60 s per load step at the start (two
alternations).

### 6.1 Check run at the paper's G_c (2026-10-19 09:40-16:55)

`decks/ni/pfd/check_L1_h020.i` (G_c 125 MPa mm, l0 0.4 mm): no failure. The full 300 steps to
60% gauge strain in 7.2 h (two alternations every step, one cutback, no memory issue); max
damage 0.030 / 0.147 / 0.290 / 0.328 / 0.400 at 5 / 20 / 40 / 46.3 / 60%, the stress within 2%
of the CP-only curve (0.985 of it at 46.3%). With the local AT2 balance d = 2 l0 H / (G_c + 2 l0
H), the peak point carries 2 l0 H = 61 MPa mm at 46.3%, so d = 0.85 there needs G_c of about
11 MPa mm at l0 0.4 (and about l0-proportional, 16 and 22 at 0.6 and 0.8, while the response
is local): an order of magnitude below the paper's value, whose l0 is not stated.

Pipeline consequence: the run could only stop on max damage (`stop_damage`), and PINNACLE
writes snapshots at the end of a run, so a failing run would either stop before the section
was cut or run on into the cracked body (where the BE fallback runaway kills it and its
snapshots are lost). PINNACLE 14fa4c6 / e2f4e94 adds `stop_separation = 'x- x+'` with
`stop_separation_damage` (default 0.85) to the coupled executioner (ADR 0031 addendum, tests):
the run ends at the first step whose d > 0.85 nodes cut the grips apart, with that step's
fields recorded and `separated = 1` in the history row. The Ni template uses it;
`ni_specimen.py curve` reports `failure_gauge_strain` from that row. Smoke (local, h040,
G_c 0.5): stopped at step 2 on separation; the snapshot check of `fracture` gives the same
step.

### 6.2 G_c(l0) calibration (launched 2026-10-19 17:40, author's go-ahead)

`decks/ni/pfd/cal_L1_h020.i` with `scripts/ni_pfd_calibrate.py` (base `decks/ni/pfd/base_L1_h020.i`):
l0 = 0.4 / 0.6 / 0.8 mm, a search per l0 in ln G_c: first round G_c = guess x0.5 and x2 (5.5 /
22, 8 / 32, 11 / 44 MPa mm; six runs side by side), then regula falsi between the closest runs
either side of 0.463 (geometric midpoint when the upper run never failed), accept within
+-0.004, at most six runs per l0. Driver checked locally on the h040 smoke deck (target 0.6%,
hit at the first run).

Rounds 1-2 (2026-10-19 17:40 - 2026-10-20 07:10), failure gauge strain (location along the gage):

| l0 (mm) | G_c (MPa mm) -> failure |
|---|---|
| 0.4 | 5.5 -> 13.9% (0.14); 22 -> none by 47.7% (oom-killed there, counted as too high); 11 -> 25.0% (0.17) |
| 0.6 | 8 -> 14.0% (0.16); 32 -> 47.9% (0.22); 27.86 -> 41.3% (0.20) |
| 0.8 | 11 -> 14.4% (0.19); 44 -> 49.9% (0.27); 38.18 -> 43.2% (0.24) |

At G_c / l0 = 14 MPa all three fail at 14%: the response is local there. Round 3: 15.56 / 30.95 /
40.75. `ni_specimen.py fracture` mapped snapshots to steps wrongly (the .pvd stamps a coupled run's
snapshots with their ordinal, and the graph search overwrote the frame index); fixed, and it now
agrees with PINNACLE's separation step on every run (the driver used summary.json, so the
calibration was not affected). Shape figure tool `scripts/ni_pfd_shape_figures.py` (deck
`decks/ni/figures/pfd_shape_L1_round2.i`, Green-Lagrange E_xx from the recorded H at the failure
snapshot, cracked cells black) with the 15 Fig. 5 panels cropped to
`data/digitized/ni_fig5_panels/`. Round-2 shapes: at every l0 the d > 0.85 zone is a band 1-2.5
mm wide, nearly straight across the gage, at 0.17-0.24 of the gage; the experiment's L.1 crack is
narrow and inclined at roughly 55-60 deg (Fig. 5), the paper's simulation also inclined.

Round 3 and a driver error (2026-10-20). l0 0.4 / G_c 15.56 separated at step 173 in the same
step whose mechanics stopped converging (the section lost its load); `curve` kept rows only up to
the first unconverged step, so the separated row was dropped, summary.json had no failure strain
and the driver read the run as "no failure" (too high a G_c) and launched G_c 13.08. The run
failed at 0.339 (last converged step; stress 205 -> 136 MPa there). Fixed: the failure is the
first separated row whatever its mechanics, at the last converged gauge strain when that row has
none; `fracture` takes the summary's strain then; the driver treats a case with summary.json as
done (a run killed before its fields has no fracture.json). The misdirected gc13.08 was stopped
and kept as `void_l0.4_gc13.08_misjudged`; every finished run was refreshed (PINNACLE's failure
step and the snapshot check agree on all of them). The l0 0.4 search resumed with
`decks/ni/pfd/cal_L1_h020_resume_l04.i` (finished runs reused) at G_c 18.5; l0 0.6 / 0.8 resume
with `cal_L1_h020_resume_l0608.i` once their round-3 runs (30.95, 40.75) end.

Round 3 and the l0 0.4 search (2026-10-20, ended 20:40). l0 0.6 and 0.8 hit in round 3, so their
resume deck was not needed: G_c 30.95 -> 0.460 (0.21 of the gage) and 40.75 -> 0.463 (0.25). Their
fracture.json had been written by the pre-fix `fracture` (the pipeline process loaded the script
before the fix), with the snapshot ordinal as step; `curve fracture` was rerun (steps 230 / 231).
l0 0.4 spent its six runs without a hit: 18.5 -> 0.402 (0.17), then 20.17 lost its load at 0.440
(stress 215 -> 141 MPa in one step, max d 0.999, 168 nodes above 0.85) and its mechanics failed on
the next step before the d > 0.85 zone spanned the section, so it stopped unseparated. The driver
read that as "no failure" (too high a G_c); with the budget spent it launched nothing more. Fixed:
a run that ends before the target without separating decides nothing, whether killed or stopped
on the mechanics. Bracket 20.17 (load lost at 0.440) / 22 (not separated by 0.477, max d 0.87):
one run at 21.1 (`decks/ni/pfd/cal_L1_h020_l04_gc21.i`).

| l0 (mm) | G_c (MPa mm) -> failure (G_c / l0, MPa) |
|---|---|
| 0.4 | 5.5 -> 0.139 (14); 11 -> 0.250 (28); 15.56 -> 0.339 (39); 18.5 -> 0.402 (46); 20.17 -> load lost at 0.440 (50); 22 -> none by 0.477 (55) |
| 0.6 | 8 -> 0.140 (13); 27.86 -> 0.413 (46); **30.95 -> 0.460 (52)**; 32 -> 0.479 (53) |
| 0.8 | 11 -> 0.144 (14); 38.18 -> 0.432 (48); **40.75 -> 0.463 (51)**; 44 -> 0.499 (55) |

With W_c = 0 the failure strain is set almost entirely by G_c / l0 (the AT2 strength scale): at
equal G_c / l0 the three length scales fail within about 0.01-0.03 of each other, the larger l0
slightly later, and the calibrated ratio is 51-52 MPa at every l0. At a fixed G_c, l0 therefore
moves the failure strain as G_c / l0 does (e.g. G_c 31 at l0 0.8 is G_c / l0 39, about 0.34), so
a fixed-G_c sweep over l0/h measures that scaling rather than settling. Shapes of the calibrated
runs (`decks/ni/figures/pfd_shape_L1_cal.i`, `_cal_damage.i`): a d > 0.85 zone 2.8 mm (l0 0.6)
and 3.2 mm (l0 0.8) wide, centered at 12.4 / 12.9 mm deformed x, flanks near the axis normal
(fracture-stage angle 87 / 73 deg); the E_xx peak sits at the lower specimen edge inside the zone.

l0 plateau check (author 2026-10-20, option 1). l0 is numerical: on the fixed h020 it is moved until
curve and failure strain stop changing, and G_c is calibrated there. Because the failure strain
follows G_c / l0 here (above), l0 is moved with G_c / l0 held: l0/h = 2, 3, 4 are the calibrated
pairs (l0/h = 2 completes with G_c 21.1), and l0/h = 6 is one run at l0 1.2, G_c 61.1 (= 50.9 x 1.2,
the l0 0.8 ratio; `decks/ni/pfd/cal_L1_h020_l12.i`, launched 2026-10-20 22:20). Curves:
`scripts/ni_pfd_curve_figures.py` (`decks/ni/figures/pfd_curves_L1_cal.i`, every calibration run in
`pfd_curves_L1_all.i`, with max d). The W_c = 0 curves sit about 5% below CP-only from 20% on (max
d 0.5-0.75 across the gage) and drop in one or two steps at failure, from about 226 MPa; the Fig. 3
experiment holds 239-241 MPa from 45 to 50% and drops at 51%. PINNACLE degrades the PK2 stress
before the resolved shear stress (`crystal_plasticity.degrade_pk2`), so damaged material keeps its
flow stress and the load falls only once g(d) is nearly zero at the section.

Plateau result (2026-10-21). l0 0.4 / G_c 21.1 failed at 0.461 (0.17 of the gage): hit. l0 1.2 /
G_c 61.1 failed at 0.472 at 0.32 of the gage, with max d only 0.88 when the d > 0.85 zone first
spanned the section (a wide, shallow band). Figures `decks/ni/figures/pfd_curves_L1_plateau.i`,
`pfd_shape_L1_plateau.i`, `_damage.i`:

| l0/h | l0 (mm) | G_c (MPa mm) | G_c / l0 | failure | location | max stress vs l0 0.8 (5-43%) |
|---|---|---|---|---|---|---|
| 2 | 0.4 | 21.1 | 52.8 | 0.461 | 0.17 | 0.08% |
| 3 | 0.6 | 30.95 | 51.6 | 0.460 | 0.21 | 0.08% |
| 4 | 0.8 | 40.75 | 50.9 | 0.463 | 0.25 | - |
| 6 | 1.2 | 61.1 | 50.9 | 0.472 | 0.32 | 0.30% |

The stress-strain curve does not depend on l0 at all at fixed G_c / l0 (the curves coincide, and
max d against strain too). The failure strain is flat over l0/h = 2-4 and rises at 6; the fracture
zone widens with l0 (2.2 / 2.8 / 3.2 / 3.4 mm at the top face) and its center drifts along the gage
toward the middle (the thinnest station is at 0.20). By the author's rule (first two consecutive
l0/h that barely change the result; take the larger) l0/h = 2 and 3 agree in curve, failure strain
and location within 0.04 of the gage, which gives l0 = 0.6 mm (l0/h = 3); l0/h = 4 still agrees in
curve and failure strain but the zone has moved to 0.25. Proposed to the author: l0 = 0.6 mm, G_c =
30.95 MPa mm for the W_c = 0 condition; the W_c stage runs at l0 = 0.6.

### 6.3 Plastic-work threshold W_c on the CP route (2026-10-20)

The threshold already exists in PINNACLE's `pfd_driving_force` (ADR 0031: H = max(H, psi+ + zeta
max(W_p - W_c, 0))) and applies to any mechanics that supplies W_p; the CP callback supplies it
(tr(sigma Fe Lp Fe^-1) det Fe over the step, also in the backward-Euler update). Checked on the
CP+PFD cube (2x2x2, finite strain, implicit update): zeta 0 -> max H 85,077 Pa, d 0.680; W_c 0 ->
98,911, 0.785; W_c 1e4 -> 88,911 (exactly 1e4 lower), 0.716; W_c 1e12 -> 85,077, 0.680 (elastic
only). PINNACLE 27a8edc adds this as a slow regression test (ADR 0031 addendum, CHANGELOG); no code
change was needed. Scale for the Ni calibration: the gage-mean plastic work along the calibrated
CP curve is about 10 / 26 / 46 / 68 / 83 MPa at 10 / 20 / 30 / 40 / 46.3% gauge strain.

### 6.4 W_c stage (launched 2026-10-21 20:10, author's go-ahead)

Author 2026-10-21: l0 = 0.6 mm fixed (the plateau, §6.2); W_c varied side by side first to learn
the trend, kept as calibration results; once the trend is known, G_c and W_c are recalibrated
together for a sharp crack initiation and propagation. Bracket: the gage-peak plastic work of
the W_c = 0 check run is 61 MPa at 38%, 73.5 at 44% and 76 at 46.3%, so W_c above about 76 MPa
cannot fail at the target; W_c = 55 / 65 / 75 MPa (damage onset near 34 / 40 / 45.5%). G_c
searched per W_c with `scripts/ni_pfd_calibrate.py` (`W_c` held per deck,
`decks/ni/pfd/wc{55,65,75}_L1_h020.i`, results under `results/ni_study/pfd/wc_L1_h020/wc<W_c>/`),
first round x0.5 / x2 around 8 / 4 / 0.6 MPa mm (from the W_c = 0 scaling: failure when W_p -
W_c is about 1.5 G_c / l0). Presentation per W_c: stress-strain curve (and max d) against the
experiment to 46.3%, failure strain, fracture outline on E_xx and d; the measure of the W_c
criterion is how long the curve follows the experiment before it leaves (the CP-only curve is
within 0.8% to 45%).

First result (2026-10-21 20:34): W_c 75 / G_c 0.3 failed at 3.7% gauge strain in the fillet
(location -0.10): with G_c / l0 = 0.5 MPa the unthresholded elastic energy psi+ alone drives the
damage (d 0.04 at 0.25% strain, where psi+ ~ 0.011 MPa and AT2 gives d ~ 2 psi+ l0 / G_c), long
before any plastic work crosses W_c. W_c gates the plastic work only. In the gage near 45% psi+ is
about 0.17 MPa, so G_c / l0 must stay a few MPa above zero (G_c >~ 2 at l0 0.6) or psi+ sets a
diffuse floor of d; that bounds the latest usable onset for each W_c. The driver read the run as
"too low" and continues; the other five first-round runs (G_c 1.2 at W_c 75; 2 / 8 at 65; 4 / 16
at 55) span the usable range.

Report §6 (W_c = 0 condition) written 2026-10-21 from the plateau result: Figs. 11-15
(`decks/ni/figures/pfd_search_L1_cal.i` is new: failure strain against G_c and against G_c / l0),
Table 6, the carried pair l0 0.6 / G_c 30.95 and what the condition does not reproduce (diffuse d,
curve 2-6% low, band normal to the axis); §7 skeleton for the W_c stage. PINNACLE gained
`elastic_prefactor` (chi) on `pfd_driving_force` (ADR 0031 addendum 2026-10-21, author's request: the
MOOSE CP-PFD carries such prefactors); the Ni pipeline passes it as `chi` (default 1, template key
`elastic_prefactor`), the driver holds it per deck. Not used by the running W_c searches.

First round, four of six runs (2026-10-21 22:40). Every small-G_c run cracked at the fillet/gage
transition long before the plastic work reached W_c: W_c 75 / G_c 1.2 (G_c / l0 2 MPa) died at
9.2% (max d 0.67, oom-killed in the mechanics fallback); W_c 65 / G_c 2 (3.3 MPa) separated at
13.8% at -0.16 of the gage (in the fillet); W_c 55 / G_c 4 (6.7 MPa) died at 12.2% with max d 0.15
(a backward-Euler fallback at one point, then the memory cap). W_c 55 / G_c 16 (27 MPa) is at 30%
and W_c 65 / G_c 8 (13 MPa) at 19%, alive. The driving force backed out of the max-d history by the
homogeneous AT2 relation H = (G_c / l0) d / (2 (1 - d)) is the same function of strain at every
G_c (0.012 MPa at 0.5%, ~0.29 at 8%, ~0.44 at 10%), i.e. independent of the damage: it is the
elastic driving energy, and it is about ten times sigma_1^2 / 2E of the nominal stress. AT2 has a
tensile strength sigma_c = sqrt(27 E G_c / (256 l0)) (103 / 205 / 265 / 375 MPa for the four pairs)
and the fillets cracked at nominal stresses of 94 / 127 / 148 / 145 MPa, a factor 1.1-2.6 below
it, consistent with the local elastic energy being that much above the nominal estimate. A check
of what the coupler passes as psi+ on the CP cube is running (scratchpad). Consequence for the
design: with chi = 1 the usable G_c / l0 is bounded below by the fillet strength (>~ 10 MPa), which
caps how late the onset can be; chi = 0 removes the elastic barrier and lets W_c alone set the
onset. The wc75 search stopped itself (gc1.2 undecided); wc55 stops after gc16 (gc4 undecided);
wc65 continues after gc8.
Check done (2026-10-21 23:30): on the CP cube the psi+ passed to H equals 1/2 E^e : C : E^e and
sigma_1^2 / 2E to within the nodal averaging; the coupler is right, and the fillet driving energy is
the local elastic energy. No PINNACLE change.

First round complete (2026-10-22 06:31). Figures `decks/ni/figures/pfd_curves_L1_wc_round1.i`,
`pfd_shape_L1_wc_round1.i`, `_damage.i` (the shape tool now draws the last snapshot of a run that
did not separate, so marked):

| W_c (MPa) | G_c (MPa mm) | G_c / l0 | end | max d at 10 / 20 / 30 / 40% | stress / experiment at 20-45% | where |
|---|---|---|---|---|---|---|
| 0 (carried) | 30.95 | 51.6 | failed 0.460 | 0.29 / 0.51 / 0.65 / 0.75 | 0.97 / 0.96 / 0.95 / 0.95 / 0.94 | gage, 0.21 |
| 55 | 4 | 6.7 | died 0.122 (max d 0.15; CP fallback, memory cap) | 0.02 ... | 1.00 | - |
| 55 | 16 | 26.7 | mechanics failed at 0.480, not separated | 0.02 / 0.06 / 0.10 / 0.38 | 0.995 / 0.993 / 0.996 / 0.992 / 0.992 | gage, diffuse |
| 65 | 2 | 3.3 | separated 0.138 at -0.16 (fillet) | 0.28 / 1 | - | fillet |
| 65 | 8 | 13.3 | mechanics failed at 0.458, max d 0.998 (28 nodes > 0.85), not separated | 0.05 / 0.13 / 0.25 / 0.54 | 0.993 / 0.988 / 0.989 / 0.983 / 0.970 | fillet (-0.1) |
| 75 | 0.3 | 0.5 | separated 0.037 at -0.10 (fillet) | - | - | fillet |
| 75 | 1.2 | 2.0 | died 0.092 (max d 0.67) | - | - | fillet |

The trend the author asked for is there: raising W_c holds the phase field down (d < 0.1 to 30%
at W_c 55 against 3% at W_c 0), the stress-strain curve returns onto the experiment (within 1% to
45% at W_c 55 / G_c 16, where the peak is 241 MPa against the experiment's 239), and the max-d
history takes the expected form: flat, a knee where the plastic work crosses W_c (37% at W_c 55),
then an accelerating rise (0.10 at 30%, 0.38 at 40%, 0.71 at 48%). But with chi = 1 every G_c small
enough to fail by 0.463 lets the elastic energy crack the fillets first (G_c / l0 <= 7 MPa), and at
the one pair that nearly reached the target (65 / 8) the damage localized at the fillet (x = 8 mm,
-0.1 of the gage), not in the gage: the elastic driving energy picks the fillet, the plastic work
the thinnest station. Two numerical notes: (i) with a sharp crack the mechanics Newton fails at
d -> 1 one step before the d > 0.85 zone spans the section (65 / 8: 0.4577, within 0.0053 of the
target; the W_c = 0 runs dropped their load first and separated), so a run that ends on a mechanics
failure with max d > 0.99 should count as failed at its last converged strain (proposed to the
author); (ii) the W_c 55 / G_c 4 death was a CP fallback at one point then the memory cap, with
max d only 0.15. Next round prepared, not launched: chi = 0 at W_c 55 / 65 / 75
(`decks/ni/pfd/wc{55,65,75}_chi0_L1_h020.i`, G_c from 8 / 4 / 0.6), which removes the elastic
barrier so that W_c alone sets the onset and the plastic work alone picks the site.

Second round launched (2026-10-22 07:20, author's go-ahead): chi = 0 at W_c 55 / 65 / 75, six
first-round runs side by side (`driver_chi0_wc*.log`). `ni_specimen.py curve` now records
failure_mode: `separation`, or `mechanics_failure_cut` when the run ends on a mechanics failure
with max d > 0.99 (`cut_damage`; the failing step's damage solve converged, so its max d counts):
W_c 65 / G_c 8 (chi 1) is thereby a failure at 0.458, step 239; W_c 55 / G_c 16 (max d 0.72) is not.

chi = 0, first results (2026-10-22 20:15; figures `decks/ni/figures/pfd_curves_L1_wc_chi0.i`,
`pfd_shape_L1_wc_chi0.i`, `_damage.i`):

| W_c | G_c | end | max d at 20 / 30 / 40 / 45% | onset (d > 0.1) | stress / experiment 20-45% | peak | zone |
|---|---|---|---|---|---|---|---|
| 55 | 4 | separated 0.476 (max d 0.86, 342 nodes > 0.85) | 0 / 0 / 0.72 / 0.84 | 35.8% | 0.996 / 0.996 / 1.00 / 0.986 / 0.984 | 237 | 0.39 of the gage, flank 48 deg, 3 mm wide |
| 65 | 8 | mechanics failed 0.484 (max d 0.70) | 0 / 0 / 0 / 0.55 | 41.3% | 0.996 / 0.996 / 1.00 / 1.00 / 1.002 | 244 | mid-gage, diffuse |
| 75 | 0.3 | mechanics failed 0.471 (max d 0.88, 70 nodes) | 0 / 0 / 0 / 0 | 45.7% | 0.996 / 0.996 / 1.00 / 1.00 / 1.008 | 242 | inclined band, 0.45-0.5 of the gage |
| 75 | 1.2 | mechanics failed 0.474 (max d 0.76) | 0 / 0 / 0 / 0 | 45.9% | same | 242 | inclined band, mid-gage |

With the purely ductile driving force every curve lies on the experiment to the end (within 1-2%
at every marker from 20 to 45%; peaks 237-244 MPa against 239), the phase field is exactly zero
until the plastic work crosses W_c (onset 36 / 41 / 46% for W_c 55 / 65 / 75) and then rises
steeply, and the damage sits in the gage, not at the fillets: the W_c 75 band is inclined at
about 45 deg across the gage at 0.45-0.5, the W_c 55 zone has a 48 deg flank at 0.39. The
experiment's crack is at 0.20 and 55-60 deg. What ends the runs now is the mechanics: at 47-48%
the staggered steps stop converging (cutback budget) whether or not the section has cut, so the
search has to bring failure in before that; the drivers read these as "G_c too high" and
continue (W_c 75 -> G_c 0.1, W_c 65 -> 2; W_c 55 regula falsi after G_c 16). The AT2 saturation
d = 2H / (G_c / l0 + 2H) means a sharp cut needs H = W_p - W_c well above G_c / l0: at W_c 75 the
plastic work beyond W_c reaches only ~2.5 MPa by 47%, so G_c / l0 must be a few tenths of a MPa
(G_c ~ 0.1-0.2) for d -> 1 in the gage within one or two steps.

chi = 0, second results (2026-10-23 00:40). W_c 65 / G_c 2 separated at 0.475 at 0.41 of the gage
with the flank at 58 deg (the experiment: 0.20, 55-60 deg); W_c 75 / G_c 0.1 cracked at the target:
d 0 -> 0.09 -> 0.54 -> 0.80 -> 0.87 -> 0.88 over steps 236-241 (45.2 -> 46.4%), the load fell from
241.5 to 227 MPa, 32 then 60 nodes above 0.85, and step 242 spent the ten cutbacks (dt halved
to 2e-5 s) before the zone spanned the section: a mechanics failure at 0.4639 with max d 0.909.
The driver read it as "no failure", judged G_c too high and started G_c 0.033; that driver was
stopped (the run continues to its end for the record) because a lower G_c only repeats the
mechanics stop. The cut rule (max d > 0.99) does not fire at 0.909. Proposed to the author: a
run whose mechanics fails while max d >= 0.85 and the load has dropped by more than 3% from its
peak counts as failed at its last converged strain (`cut_damage` 0.85 plus a load-drop test);
then W_c 75 / G_c 0.1 fails at 0.464, the target. Figures redrawn
(`pfd_curves_L1_wc_chi0.i`: 55 / 4, 65 / 8 and 2, 75 / 1.2, 0.3, 0.1; shapes: 55 / 4, 65 / 2, 75 / 0.1
and 0.3; the 75 / 0.1 field is its last snapshot, step 240 at 46.1%, before the last two steps).

Author 2026-10-23: option 1 adopted and final (load drop > 3% with max d >= 0.85 = full
separation; `failure_mode = load_drop_cut`, `fracture.json` mode); W_c 75 / G_c 0.1 fails at 0.464
at 0.38 of the gage (18 nodes above 0.85 in its last snapshot, step 240). Crack line added to the
shape figure (`ridge = true`): the max-d node of each row across the gage width, joined and
colored by d, over the 0.85 outline; on the chi = 0 runs the line is inclined at roughly 45-60
deg across the width at 0.38-0.41 of the gage.

chi = 0 searches, state at 2026-10-23 11:00 (W_c 55 / G_c 3.255 still running; figures redrawn with
every run, `pfd_curves_L1_wc_chi0.i`, `pfd_shape_L1_wc_chi0.i`, `_damage.i`, ridge line on):

| W_c | G_c | failure | mode | onset (d > 0.1) | S / exp at 20 / 30 / 40 / 45% | peak | zone center | flank |
|---|---|---|---|---|---|---|---|---|
| 55 | 1.333 | 0.409 | separation | 35.7% | - | 225.5 | 0.39 | 56 deg |
| 55 | 4 | 0.476 | separation | 35.8% | 0.996 / 0.996 / 0.986 / 0.984 | 237.1 | 0.39 | 48 deg |
| 55 | 16 | none by 0.601 (path end) | - | 36% | - | 252.4 | - | - |
| 65 | 0.6667 | 0.440 | separation | 40.7% | - | 234.8 | 0.40 | 60 deg |
| 65 | **1.384** | **0.459** | separation | 40.7% | 0.996 / 0.996 / 1.00 / 0.978 | 234.9 | 0.41 | 58 deg |
| 65 | 2 | 0.475 | separation | 40.7% | - | 236.4 | 0.41 | 58 deg |
| 65 | 8 | none by 0.484 (mechanics) | - | 41.3% | - | 243.8 | - | - |
| 75 | 1.2 / 0.3 | none by 0.474 / 0.471 (mechanics; max d 0.74 / 0.86) | - | 45.9 / 45.7% | 0.996 / 0.996 / 1.00 / 1.008 | 242 | - | - |
| 75 | **0.1** | **0.464** | load drop (5.9%) | 45.7% | 0.996 / 0.996 / 1.00 / 1.008 | 241.5 | 0.38 | 77 deg (zone); crack line inclined |
| 75 | 0.03333 | 0.460 | load drop (17.2%) | 45.5% | same | 241.4 | 0.38 | 83 deg (zone); crack line inclined |

W_c 65 is calibrated by the driver (G_c 1.384, 0.459, within tolerance). W_c 75: G_c 0.1 (0.464,
the author's preferred result) and 0.033 (0.460) both hit under the final definition; G_c 0.1 is
the carried value (the larger G_c with the smaller load drop at the cut). W_c 55 awaits G_c 3.255
(bracket 1.333 -> 0.409, 4 -> 0.476). The W_c 55 / G_c 16 run went to the path end at 60% with no
mechanics trouble, so the 47-48% mechanics stops of the other runs are damage-related. The
`angle_to_axis_deg` of the zone (PCA of the d > 0.85 nodes) reads 77-83 deg for the W_c 75 cases
because their zone is a short segment at one row; the crack line (max-d ridge) is what shows
the inclination there.

