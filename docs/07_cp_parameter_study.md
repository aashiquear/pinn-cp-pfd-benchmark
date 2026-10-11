# CP parameter study on the x-loaded plate

Lab-notebook record. Design: `docs/workflow_polycrystal_plate.md` (pipeline)
and `decks/cp_study/` (decks). Every run is in `docs/run_ledger.md`. Host
rocm-server, numpy, float64, one core per run. Ground truth: the nine
digitized experiment markers of Rezwan (2025) Fig. 3 (Vaughan 2024).

Plate 100 x 50 x 1 um loaded along x, one element through thickness, mesh
level hN = N elements along x. Grain structure: the 30-grain phase-field
structure of `results/graingrowth_pilot/n32_bdf1_h40` (2.5 interface times),
orientation seed 20261009, mapped to every mesh by nearest point. Loading:
20% nominal strain at 1e-3 /s (assumed rate) in 256 steps.

## 1. Base case, reference parameter set, h40 (`decks/cp_study/h40_base.i`)

Run `cp-h40-base`: 38 min wall; **the solve reaches 7.77% strain and then
exhausts the constitutive sub-step budget** (100,000 forward-Euler
sub-steps at |dgamma| <= 1e-3) at load step 100 of 256, through ten dt
cutbacks down to 7.6e-4 s. 125 of 258 history rows are converged, 133 are
failed attempts. Up to the wall: zero cutbacks, 3-6 Newton iterations.

| nominal strain | mean sigma_xx (MPa) |
|---|---|
| 0.08% | 55.5 |
| 1.33% | 773.7 |
| 2.58% | 864.0 |
| 5.08% | 1002.6 |
| 7.58% | 1123.9 |

Everything after step 100 is NOT a result: failed steps advance the load
without advancing the state, so the later "converged" rows climb at 60-75
GPa per unit strain, i.e. elastically, to 3,970 MPa at 12.8%. The curve
stage now truncates at the first failed step. Figure
`figures/cp_study_h40_base.png` shows the full history for the record.

Two readings. (i) With the reference set under the latent-hardened Taylor
sum the plate is at 1,124 MPa at 7.6% and still hardening at about 5 GPa per
unit strain; the Kocks-Mecking saturation of this set is rho_sat =
(k1 / 2 y_c)^2 = 4,100 um^-2, i.e. tau_c,sat = 762 MPa with the latent sum
(macroscopically ~2,300 MPa), far above the experiment's 590 MPa plateau.
The calibration will have to bring k1 / y_c down by an order of magnitude
in saturation strength, whichever density convention is used. (ii) The
wall at 7.8% is the same signature the first campaign hit on holed plates
at 6-10% (`archive/docs/04_seventeen_percent_blocker.md`): the explicit
forward-Euler constitutive update cannot absorb the increment at some
quadrature point however small dt is made. That is a framework defect to
fix before any 20% run, independent of parameters; see §2.

**Base case rerun with the direct solver** (`cp-h40-base (direct)`, 68 min
while sharing the host with six other runs): 256/256 steps to 20%, zero
cutbacks. Engineering stress 404 / 815 / 1,030 / 1,147 / 1,257 MPa at 0.57
/ 2 / 7 / 11 / 17% (experiment 345 / 412 / 502 / 553 / 593). The Taylor
surrogate for the same set gives 407 / 881 / 1,164 / 1,320 / 1,475: the FE
polycrystal is 0.99 of Taylor at yield and 0.85 at 17%, which is the
offset the FE correction of §4 has to absorb. Figure
`figures/cp_study_h40_anchors.png`.

## 2. Diagnosis of the sub-step wall

**h20 twin** (`decks/cp_study/h20_base.i`, 20 x 10 x 1, run `cp-h20-base`,
9 min): walls at step 126, 9.83% strain, 1,227 MPa. Up to the wall the two
meshes agree closely (998.6 vs 1,032.0 MPa at 5%, 3.3%), so the pre-wall
curve is close to mesh-converged and the wall moves EARLIER with refinement
(7.8% at h40 against 9.8% at h20) -- the signature of a local stress
concentration the finer mesh resolves better, at a grain-boundary junction
or the loaded edge, not of a global instability. Figure
`figures/cp_study_base_h20_h40.png`.

**The wall is the linear solver.** With PINNACLE `bb774d0`, which logs the
cause of every cutback, the instrumented h20 rerun (`cp-h20-diag`) shows at
step 126: `Newton did not converge ... CG did not converge (info=1386)` on
every cutback. The constitutive sub-stepper never raised. The load stepper's
old cutback line said "CP sub-step budget exhausted" for both causes, which
is how the first campaign's "sub-step wall" was read; whether that campaign's
holed-plate failures were also CG failures is unknown and worth re-checking.

Conjugate gradients requires a symmetric positive-definite matrix. The
finite-strain CP tangent is not symmetric in general and loses positive
definiteness as slip localizes (geometric softening of the plastic
modulus), so CG stalls at the strain where the first points soften --
earlier on the finer mesh, which resolves the concentration better. The
study decks used `linear_method = cg` because the first campaign measured
the sparse direct solve scaling as n^2.91 on 3-D void meshes; this plate
has 7,749 (h40) / 2,079 (h20) displacement DOFs, where a sparse LU is cheap.
Test: `cp-h20-direct`, the h20 base case with `linear_method = direct`:
**256 of 256 steps to 20% strain, zero cutbacks, 15 min** (1,618 MPa at
20%, still hardening -- see §3 for why). The wall was the solver choice.
`direct` is now the study default (`linear_method` in the pipeline deck's
[Mechanics], written into the template); PINNACLE `47a41e0` adds
`linear_method = gmres` (ILU-preconditioned) as the iterative option for
meshes too large for a sparse LU. The CG runs are kept as
`results/cp_study/h40_base_cg` and `h20_base_cg`. Figure
`figures/cp_study_h20_direct_vs_cg.png`.

Open question carried: whether the first campaign's holed-plate "sub-step
wall" (6-10%) was also CG. Its decks used `linear_method = cg` too.

## 3. Why the reference set does not saturate: the density convention (C.15)

The author asked why the stress keeps rising. The hardening law is working
as parameterized: a Taylor material point (30 grains, same seed, archived
`calibrate_cp_polycrystal.taylor_curve`) with the identical set gives 879 /
1,068 / 1,164 MPa at 2 / 5 / 7% against the plate's 864 / 1,003 / 1,124.
What is far away is saturation. With tau_c = tau_c0 + alpha0 mu b
sqrt(sum_b a_ab rho_b) and the latent-hardened sum (a_ab = q = 1.4 off the
diagonal, 12 FCC systems), every density is multiplied by 1 + 11 q = 16.4
inside the square root, a factor 4.05 on tau_c:

| quantity | self term | latent sum |
|---|---|---|
| tau_c(rho0 = 220 um^-2) | 43.7 MPa | 176.7 MPa |
| rho_sat = (k1 / 2 y_c)^2 | 4,100 um^-2 | 4,100 um^-2 |
| tau_c(rho_sat) | 188 MPa | 762 MPa |
| plateau, Taylor factor 3.07 | 578 MPa | 2,340 MPa |
| shear to saturate from rho0 | ~0.28 | ~0.28 |

The self-term plateau is the experiment's 590 MPa plateau; the paper's
Table 2 is therefore a self-term set (as `03_material_point_verification.md`
concluded for C.15). Under the latent sum, chosen for this study, the same
numbers are 4.05x too strong and saturate only after ~28% shear, which is
exactly the curve in §1.

**Latent-sum equivalent of the paper's set.** Requiring the same tau_c at
every point of the path: rho' = rho / 16.4 and, so that d rho'/d gamma
reproduces tau_c(gamma), k1' = k1 / sqrt(16.4) = k1 / 4.05 with y_c and
tau_c0 unchanged (then rho_sat' = rho_sat / 16.4 and tau_c,sat' = tau_c,sat).

| parameter | paper (self term) | latent-sum equivalent |
|---|---|---|
| rho0 | 220 | 13.41 um^-2 |
| k1 (k_m) | 0.3432 | 0.08475 |
| y_c | 0.00268 | 0.00268 um |
| tau_c0 | 0.112 | 0.112 MPa |

Run `cp-h40-rescaled` (direct solver, 72 min, zero cutbacks): engineering
stress 210 / 232 / 281 / 307 / 332 MPa at 0.57 / 2 / 7 / 11 / 17%, i.e.
0.6x the experiment with the wrong hardening slope. Superseded by §4: the
calibration targets the experiment directly.

Taylor (upper-bound) prediction for the rescaled set at 1e-3 /s, engineering
stress: 216 MPa at 0.56%, 242 at 2%, 306 at 7%, 342 at 11%, 375 at 17%,
against the experiment's 345 / 412 / 502 / 553 / 593. So the paper's
numbers do not reproduce the experiment in this framework under EITHER
reading: 2.3x too strong with the latent sum as written, about 0.6x too
weak once re-expressed. Two framework-specific factors sit between the
table and the curve: the density convention (4.05x on tau_c) and the rate
sensitivity (at gamma_dot0 = 1e-5 /s, n = 8.33 and a shear rate of ~3e-3
/s the flow stress is 1.98x tau_c; the experiment's rate is not recorded).
Both are the kind of model- and framework-sensitivity the author
anticipated. The calibration therefore starts from the rescaled set and
raises rho0 (yield level) and k1 (hardening) in the agreed order; the
self-term check at the material point (q_lat = 0) could not be run
directly because the explicit Taylor integrator diverged at that lower
strength, which is itself the sub-step weakness of §2.

## 4. Calibration (decision 2026-10-09 evening)

Objective: the parameter set that minimizes the error against the
experimental curve in THIS framework at THIS mesh (h40). Neither the paper's
numbers nor their re-expression is a target; both are framework- and
model-dependent. The MOOSE decks are convergence tests of the input and are
not compared against further.

Conventions, stated once: stress is ENGINEERING (force over original
area) on both sides -- the markers are from a tensile test, and the FE mean
Cauchy stress is converted with the two mean transverse stretches
(`polycrystal_plate.py curve`, column `eng_stress`); strain rate 1e-3 /s
(assumed, not reported); latent-hardened sum; error = RMSE over the nine
markers plus the per-marker relative error.

Free, in the author's order: rho0, k1 (k_m), y_c, tau_c0. Bounds from the
experiment's features (`decks/cp_study/taylor_screen.i`): rho0 20-80
um^-2, k1 0.02-0.12, y_c 0.001-0.006 um, tau_c0 0-50 MPa.

Procedure: (a) Taylor material-point surrogate (`scripts/taylor_surrogate.py`;
30 grains, same seed as the FE orientations are NOT shared -- the surrogate
has its own 30 random grains) for OAT sensitivities and a first fit inside
the bounds; (b) FE at h40, one-at-a-time around that point and a finite-
difference Gauss-Newton correction, five concurrent runs per iteration;
(c) error analysis of the final set.

### 4a. Surrogate stage (`results/cp_study/taylor_screen/`)

OAT at the start point (rho0 40, k1 0.05, y_c 0.00268, tau_c0 5), engineering
stress at the marker strains (0.56 / 2 / 7 / 11 / 13 / 15 / 17%), experiment
345 / 412 / 502 / 553 / 574 / 588 / 593:

| case | 0.56% | 2% | 7% | 11% | 17% | RMSE |
|---|---|---|---|---|---|---|
| start | 344 | 403 | 400 | 396 | 386 | 130 |
| rho0 = 20 / 80 | 270 / 399 | 281 / 553 | 304 / 532 | 313 / 514 | 317 / 489 | 194 / 75 |
| k1 = 0.02 / 0.12 | 344 / 344 | 394 / 426 | 360 / 503 | 342 / 545 | 318 / 584 | 173 / 15 |
| y_c = 0.001 / 0.006 | 344 / 344 | 410 / 391 | 437 / 352 | 454 / 332 | 474 / 310 | 79 / 180 |
| tau_c0 = 0 / 50 | 329 / 404 | 377 / 622 | 375 / 616 | 372 / 603 | 364 / 583 | 147 / 85 |

What each parameter does (this framework, this law): rho0 moves the yield
level and, through the Taylor term, the whole curve up; k1 sets the slope
after yield and the plateau height; y_c lowers the plateau (recovery) and
with it the late slope; tau_c0 is a uniform offset that also raises the
early hardening because the slip resistance's rate term scales with it.
With k1 at the first bound (0.12) the fit was already at 15 MPa RMSE, so the
bound was widened to 0.30.

Bounded Gauss-Newton (log parameters, 6 iterations, 4 s per curve set):
**rho0 = 37.19 um^-2, k1 = 0.1256, y_c = 0.002666 um, tau_c0 = 5.43 MPa,
RMSE 13.2 MPa** (`fit.json`). The paper's y_c survives essentially
unchanged; rho0 and k1 are 0.17x and 0.37x the table's values under the
latent sum. This is the FE starting point (`decks/cp_study/cal_center1.i`),
FE iteration 1 is `decks/cp_study/cal_gn1.i`.

### 4b. FE iteration 1 at h40 (`decks/cp_study/cal_gn1.i`, `results/cp_study/cal/gn1/`)

Five runs (center = surrogate optimum, one perturbation each: rho0 x1.2, k1
x1.2, y_c x1.2, tau_c0 x1.5), 256 steps to 20% each, all converged, zero
cutbacks, ~75 min wall for the set (seven runs shared the host).

Center against the experiment (engineering stress, markers from 0.1%):

| strain | experiment | FE center | error |
|---|---|---|---|
| 0.14% | 138 | 99 | -28.0% |
| 0.57% | 345 | 338 | -2.0% |
| 2.02% | 412 | 384 | -6.9% |
| 7.01% | 502 | 454 | -9.7% |
| 11.02% | 553 | 491 | -11.2% |
| 13.05% | 574 | 505 | -12.0% |
| 15.06% | 588 | 516 | -12.3% |
| 16.98% | 593 | 525 | -11.5% |

RMSE 53.8 MPa: the FE polycrystal is softer than the surrogate,
as the anchors showed, and the deficit grows with strain.

Sensitivities dS/d(ln p) in MPa at 2 / 7 / 17% (what each parameter does
in the FE plate): rho0 155 / 132 / 98 (level, fading with strain); k1 36 /
172 / 328 (hardening, growing with strain); y_c -11 / -70 / -177 (recovery,
lowers the late curve); tau_c0 28 / 28 / 27 (uniform offset).

Levenberg-Marquardt step (log parameters, bounds of the surrogate deck):
next center rho0 41.15, k1 0.1677, y_c 0.00320 (at its
x1.2 perturbation), tau_c0 5.35; linear prediction RMSE
14.0 MPa. Iteration 2 is `decks/cp_study/cal_gn2.i`.

### 4c. FE iteration 2 and the calibrated h40 set (`results/cp_study/cal/gn2/`)

Center rho0 41.15, k1 0.1677, y_c 0.0032, tau_c0 5.35; five runs, all
converged, zero cutbacks. The next Levenberg-Marquardt step moves no
parameter by more than 2.5% (rho0 41.20, k1 0.1716, y_c 0.00325, tau_c0
5.14) and predicts no gain (14.0 vs 14.3 MPa), so the iteration is stopped.

**Calibrated set (h40, latent-hardened sum, 1e-3 /s, engineering stress):
rho0 = 41.15 um^-2, k1 = 0.1677, y_c = 0.0032 um, tau_c0 = 5.35 MPa.**
Fixed by decision: C11 108.0, C12 61.3, C44 28.5 GPa; mu 25.7 GPa; b
0.000286 um; q 1.4; alpha0 0.4; gamma0 1e-5 /s; n = 8.333.

Error analysis:

| strain | experiment (MPa) | FE (MPa) | error |
|---|---|---|---|
| 0.14% | 138 | 99 | -28.0% |
| 0.57% | 345 | 350 | +1.6% |
| 2.02% | 412 | 407 | -1.1% |
| 7.01% | 502 | 502 | +0.0% |
| 11.02% | 553 | 550 | -0.4% |
| 13.05% | 574 | 567 | -1.1% |
| 15.06% | 588 | 581 | -1.3% |
| 16.98% | 593 | 591 | -0.2% |

RMSE over the eight markers from 0.1%: 14.3 MPa; over the seven from
0.5%: 4.6 MPa (max |error| 7 MPa, 1.6%).
The 0.14% marker is the one miss (-28%): it lies on the digitized curve's
initial slope of ~99 GPa, while the random-texture plate's elastic modulus
is 71 GPa (Al's polycrystal value is ~70 GPa), so that point reflects the
test record or its digitization rather than a parameter of this law.

Linearized 1-sigma parameter uncertainty from the post-0.5% residual (4 MPa)
and the FE Jacobian (a guide to how tightly the data pins each value):

| parameter | value | relative 1-sigma |
|---|---|---|
| rho0 | 41.152 | ±59% |
| k1 | 0.16771 | ±11% |
| y_c | 0.0031974 | ±20% |
| tau_c0 | 5.3492 | ±343% |

Sensitivities at the calibrated point, dS/d(ln p) in MPa at 2 / 7 / 17%:
rho0 160 / 132 / 96; k1 47 / 222 / 403; y_c -14 / -94 / -228; tau_c0 28 /
28 / 27. rho0 sets the yield level, k1 and y_c the hardening slope and its
saturation, tau_c0 a uniform offset: the author's reading of Fig.
`cp_cal_gn1.png`, confirmed by the Jacobian. Figure of the calibrated curve:
`figures/cp_cal_final_h40.png`.

Caveats carried: the strain rate is assumed; the markers are digitized; the
set is for THIS mesh (h40) and this grain realization. §5 tests the mesh.

## 5. Mesh sensitivity with the calibrated set

Decks `decks/cp_study/mesh_h{20,80,160}.i` (h40 is the iteration-2 center).
Same structure, orientations, loading and set; only nx x ny changes.
Engineering stress (MPa) at 2 / 7 / 17 / 20% and RMSE against the markers
from 0.1%:

| mesh | elements | h (um) | elements across a grain (13.7 um) | 2% | 7% | 17% | 20% | RMSE | wall |
|---|---|---|---|---|---|---|---|---|---|
| h20 | 200 | 5.0 | 2.7 | 419.8 | 523.3 | 621.9 | 636.4 | 23.3 | 16 min |
| h40 | 800 | 2.5 | 5.5 | 406.9 | 502.4 | 591.5 | 603.9 | 14.3 | ~40 min |
| h80 | 3,200 | 1.25 | 11 | 401.4 | 492.4 | 575.3 | 586.4 | 19.9 | 259 min |
| h160 | 12,800 | 0.625 | 22 | (running) | | | | | |

h20 -> h40: -3.1% at 2%, -4.0% at 7%, -5.1% at 20%. h40 -> h80: -1.4%,
-2.0%, -2.9%. The coarse mesh is stiffer, as a displacement-based mesh with
2.7 elements across a grain should be, and each halving of h removes about
0.55 of the previous change, so the remaining change beyond h80 is
expected to be about 3% at 20% (1.7% at 7%); h160 tests that. Figure
`figures/cp_mesh_sensitivity.png`. The RMSE rises away from h40 because the
set was calibrated AT h40; the refit on the chosen mesh (next stage) is a
small, mostly rho0/k1 correction of the same size as the mesh shift.
