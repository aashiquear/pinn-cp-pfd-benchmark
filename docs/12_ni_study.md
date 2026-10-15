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
| implicit, cap 1e-5 (single element, 16 steps of 2.5e-4) | -- | within 0.5% of a 1,024-step explicit path | -- | 12 s |

The explicit scheme is not monotonic in the step count because its two
errors have opposite signs (the trial stress overshoots, the recorded
relaxed-state reaction undershoots) and cross between 128 and 256 steps;
it needs increments near 1e-4 for 5% and is 4% low at 3.7e-5. The
iterated scheme at a 1e-4 sub-increment cap is not better than the
explicit one at the same step count; at 1e-5 it is +8% at 1.5e-4
increments. The implicit scheme is accurate at any step size on the
single element but costs an integration per residual evaluation (about 25
s per load step on this 1,760-point mesh).

Decision: PINNACLE gets a consistent stress update (ADR 0033): with
`state_update = implicit` the callback's `stress()` integrates the state
from the committed one over the step's time increment (the existing
sub-stepper) at every residual evaluation, and the tangent is the
central-difference derivative of that stress per quadrature point; the
default stays explicit so Parts I-II reproduce. Validation: the implicit
FE at 64 steps must reproduce the explicit limit (87 MPa at 3%) and a
single element under uniaxial stress must match the material point.
