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
is the centre, `decks/ni/cp/cal_gnN.i` runs it with one multiplicative
perturbation per free parameter (x1.2; tau_c0 x1.5) through
`scripts/ni_study.py`, and `scripts/ni_gauss_newton.py` takes the damped
step in log space against the Fig. 3 experiment markers
(`data/digitized/ni_fig3_experiment_markers.csv`, engineering stress
against GAUGE strain, markers below 0.3% strain excluded) and writes the next
centre and GN deck. Runs on the host, torch CPU, 2 threads, 5 side by side,
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
that centre, and the clipping means two or three more iterations.

### 2.2 Iteration 2 (2026-10-15 21:15 - 2026-10-16 00:00)

Centre rho_0 1.17e6, k_m 0.0146, y_c 0.74e-6, tau_c0 0.255:

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
centre (clipped step) rho_0 7.64e5, k_m 0.0108, y_c 1.48e-6, tau_c0 0.249,
linear prediction 4.6 MPa; iteration 3 runs from it.

### 2.3 Iteration 3 (2026-10-16 00:00-01:50)

Centre rho_0 7.64e5, k_m 0.0108, y_c 1.48e-6, tau_c0 0.249:

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
31, k_m 9 / 56 / 156, y_c -3 / -23 / -79, tau_c0 3 / 3 / 3. Next centre
rho_0 6.86e5, k_m 0.0111, y_c 2.26e-6, tau_c0 0.246 (linear prediction
4.5 MPa, the markers' own spread); iteration 4 decides.

### 2.4 Iteration 4: converged (2026-10-16 01:50-08:15)

Centre rho_0 6.86e5, k_m 0.0111, y_c 2.26e-6, tau_c0 0.246:

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
predicts 4.5 MPa, so the calibration is converged and the iteration-4 centre
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

`decks/ni/cp/h020/cal_gn1.i`: centre = the h040-calibrated set on h020 (the mesh
study's h020 run, copied in as `base`), four perturbations at 2 threads each,
memory cap 16 GB per run; `ni_gauss_newton.py` now writes the next decks next to the
given one (`decks/ni/cp/h020/`, results `results/ni_study/cal_h020/`). About 4 h per
iteration.

### 4.1 Iteration 1 on h020 and a third PINNACLE bug (2026-10-17)

Centre (the h040 set on h020): RMSE 29.7 MPa, the curve 12-18% below the markers
from 10% strain. The k_m x1.2 run cut back twice at 14.6% gauge strain and its
curve then sat 46 MPa above the centre at 20% (h040's sensitivity predicts about
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
are relied on). Iteration 1 step (centre RMSE 29.7 MPa): sensitivities at 2 / 10 / 30%
rho_0 38 / 37 / 25, k_m 9 / 47 / 126, y_c -3 / -26 / -78, tau_c0 3 / 3 / 3 (k_m as on
h040); next centre rho_0 7.26e5, k_m 0.01147, y_c 1.54e-6, tau_c0 0.244 (predicted
3.6 MPa). Iteration 2 launched 20:20 (`decks/ni/cp/h020/cal_gn2.i`, 5 x 2 threads).

### 4.2 Iteration 2 on h020: converged (2026-10-17 20:20 - 2026-10-18 04:35)

Centre rho_0 7.256e5, k_m 0.01147, y_c 1.542e-6, tau_c0 0.244: RMSE 4.1 MPa (from
29.7), the next step predicting 3.4 MPa; every marker from 2% to 45% within -1% to
+6%, 0.5% at -7%, 50% at +4% (the experiment's softening before failure). This is the
calibrated set on h020 (`decks/ni/cp/h020/calibrated_L1_h020.i`); the further step is
not taken. Sensitivities at 2 / 10 / 30%: rho_0 39 / 39 / 27, k_m 10 / 59 / 161, y_c
-3 / -26 / -83, tau_c0 3 / 3 / 3. The centre and tau_c0 runs cut back 5 and 4 times
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
h010 keeps running for its field snapshots (written at the end; valid to step 240).
