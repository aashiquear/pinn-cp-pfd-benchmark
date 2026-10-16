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
