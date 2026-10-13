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
has 5,166 (h40) / 1,386 (h20) displacement DOFs (41 x 21 x 2 and 21 x 11 x 2 nodes), where a sparse LU is cheap.
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
| h160 | 12,800 | 0.625 | 22 | 398.3 | 486.0 | 563.5 | 573.5 | 26.0 | 1,063 min |

h20 -> h40: -3.1% at 2%, -4.0% at 7%, -5.1% at 20%. h40 -> h80: -1.3%,
-2.0%, -2.9%. h80 -> h160: -0.8%, -1.3%, -2.2% (h160 run `mesh-h160`: 256/256
steps, zero cutbacks, 17.7 h wall on one core). The coarse mesh is stiffer,
as a displacement-based mesh with 2.7 elements across a grain should be, and
every halving softens the plate further. Figure
`figures/cp_mesh_sensitivity.png`.

The convergence is NOT asymptotically clean. The ratio of successive changes
(h80->h160 over h40->h80) is 0.56 at 2% strain but 0.74 at 20%, i.e. an
apparent order of 0.8 at yield falling to 0.44 at 20% strain. At yield the
response converges like a first-order displacement mesh; on the plastic
branch the finer meshes keep resolving sharper slip bands at grain
boundaries and junctions (the CP law has no internal length, so localized
slip is mesh-dependent), and the change per halving decays slowly. A
geometric extrapolation with the local ratio puts the remaining change beyond
h160 at ~4 MPa (1%) at 2% but ~36 MPa (6%) at 20%, and that number is itself
uncertain because the ratio is still rising. Elapsed-time cost: h40 ~40 min,
h80 4.3 h, h160 17.7 h per run.

Reading for the mesh decision (the author's call): the per-halving change is
at most 2.2% at h80->h160 and 2.9% at h40->h80 over the whole curve, so h40
is within ~5% and h80 within ~3% of h160 everywhere, and no mesh in reach is
converged to better than a few percent at 20% strain. The RMSE rises away
from h40 only because the set was calibrated AT h40; the refit on the chosen
mesh is a small, mostly rho0/k1 correction of the size of the mesh shift
(about +2% on rho0 and +3% on k1 per halving, from the §4c sensitivities).

## 6. Parameter sensitivity at h40 (decision 2026-10-11: h40 is the mesh)

The author chose h40 on the within-5% criterion (§5). Study deck
`decks/cp_study/sens_h40.i`: the calibrated center (`cal_center2.i`) and one
parameter at a time at x0.5, x0.8, x1.25 and x2.0 of its calibrated value --
rho0, k1, y_c, tau_c0, sixteen cases plus the base, ten concurrent. Collected
with `cp_parameter_study.py collect` against the markers (RMSE, per-marker
error, stress at 2 / 10 / 17%). Results (all 17 cases 256/256, 0 cutbacks, ~3 h wall for the set at ten
concurrent): the full table is in the report §6 (Table 6). Yield 272 / 331 /
364 / 401 / 495 MPa for rho0 x0.5 / 0.8 / 1 / 1.25 / 2 (sigma_y ~ rho0^0.43);
20% stress 418 / 527 / 604 / 702 / 1,002 MPa for k1 and 788 / 663 / 604 / 548
/ 445 MPa for y_c over the same multiples; tau_c0 shifts the curve 15-32
MPa per e-fold with no change of shape. RMSE ranking at x2: k1 247, y_c 98,
rho0 85, tau_c0 22 MPa. Figures `figures/report/sens_panels.png`,
`sens_summary.png`; summary `results/cp_study/sens/h40/sens_h40_summary.csv`.

## 7. Grain-instance check

Two further grain structures from new Voronoi seedings (uniform in the
plate, rng seeds 20261011 and 20261012, decks written by
`scripts/make_growth_instance.py` from the pilot growth deck and grown with
the identical settings), each with its own orientation draw (same seeds), run
with the calibrated set at h40 (`decks/cp_study/grain_s{A,B}.i`). The curves
are compared with the calibrated structure's to show the set is not tied to
one realization.

First pair (2026-10-11): instance B (30 grains, seed 20261012) is within 2%
of the calibrated structure over the whole curve (593.2 vs 603.9 MPa at 20%,
-1.8% mean). Instance A (32 grains, seed 20261011) was 17.5% stronger from
yield onward (729.5 MPa at 20%); its area-weighted Taylor estimate along x
is only 4.5% above the calibrated set's, so most of the offset was the grain
arrangement, not the orientation draw. The author's decision: such severe
outlier realizations occur and are dropped; instance A and the
texture-swap diagnostic prepared for it were removed, and a replacement
instance (seed 20261013, `grain_sC.i`) is run instead. Instance B's result
is accepted.

Instance C (seed 20261013; 31 grains, d_med 13.7 um): 256/256, 0 cutbacks,
E 71.5 GPa, yield 372 MPa, 419 / 524 / 622 / 635 MPa at 2 / 7 / 17 / 20%,
+4.4% mean (5.2% max) above the calibrated structure; RMSE 23.6 MPa. With
instance B at -1.8% the accepted spread over three realizations is about
+-3% on the flow stress. Report §7, Table 7, Figs. 12-13.

## 8. Field series for the plastic-strain animation

PINNACLE 6517a31 adds `snapshots_every` to `cp_quasi_static`; the four mesh
runs were repeated with `snapshots_every = 4` (`decks/cp_study/fields_h*.i`,
64 frames to 20%) so the equivalent plastic strain can be animated
(`scripts/plastic_strain_animation.py`, deck
`decks/report/plastic_strain_animation.i`). The stress-strain curves of the
reruns must equal the originals; that equality is checked when they land.

## 9. Framework fix found while batching the CP path (2026-10-11)

Making the CP evaluation batched for the GPU (PINNACLE a311a8b, ADR 0027)
exposed, through a batched-versus-per-point parity test, that the batched
sub-stepper on an orientation-field deck had been subsetting the per-point
directions, normals and Schmid tensors but NOT the per-point rotated
stiffness: the CP state evolved under the material-global elastic law while
the residual (stress and tangent) used the rotated one. Every polycrystal FE
run of this study took that path. Al's cubic anisotropy is mild (Zener ratio
1.22), so the effect on the curve is expected to be small, but it must be
measured, not assumed: `check_cfix_h40.i` repeats the calibrated h40 run with
the fix.

**Result (run `check-cfix-h40`, 256/256, 0 cutbacks):** the curve moves by
far more than the mesh or instance spread. Engineering stress before -> fixed:
0.57% 350.8 -> 332.5 (-5.2%), 2% 406.9 -> 375.1 (-7.8%), 7% 502.4 -> 457.9
(-8.9%), 11% 550.1 -> 499.3 (-9.2%), 17% 591.5 -> 535.0 (-9.6%), 20% 603.9
-> 545.6 (-9.7%). RMSE against the markers 14.3 -> 48.1 MPa. The anisotropy
argument was wrong: the elastic stiffness enters the resolved shear stress
at every grain, and a ~10% error in the driving stress integrates into the
hardening over the whole path. Which curve is right: the fixed one. Before
the fix the residual and tangent used the rotated stiffness while the state
evolved under the material-global one -- the solve was inconsistent by
construction; after it the batched path equals the per-point path to 1e-10
(PINNACLE test `test_cp_batched_constitutive_path.py`).

**Consequence:** every polycrystal FE result of this study (base case,
calibration iterations, calibrated set, mesh sensitivity, parameter
sensitivity, grain instances, field series) was produced on the inconsistent
path and is superseded. The decks are unchanged and the framework is now
4-5x faster (ADR 0027), so the campaign can be repeated as it stands:
calibration (two Gauss-Newton iterations of five h40 runs, ~10 min each),
mesh h20/h40/h80/h160, sensitivity (17 runs), grain instances (2 runs),
field series (4 runs). Decision requested from the author before launching
(standing rule). The h160 field-series run still in progress on the
pre-fix code was stopped (its output would have been superseded too).

## 10. Regeneration with the fixed solver (from 2026-10-12)

Author's decision: every result of the report is regenerated with PINNACLE
at a311a8b or later; nothing computed before it stays in `results/` or in
the report, and every run writes the field series (`snapshots_every = 4`).
Order: calibration (restarted from the previous set as the center, iteration
decks `cal_center{N}.i` / `cal_gn{N}.i`), effect of the parameters and
error analysis, mesh sensitivity h20/h40/h80/h160 with the new set; then
the report and the author's mesh decision; parameter sensitivity, grain
instances and the PINN / PINN+FEM comparisons follow on that mesh. Sections
1-9 above are the record of the superseded campaign and are not reused.

### 10.0 Convergence criterion of every run (decided 2026-10-11)

PINNACLE's Newton stop test is ``||R|| <= max(rel_tol ||R_0||, tol)`` with
``||R_0||`` the first iterate's residual. That test is unfair to a warm
start: a predictor that lands within 0.03% of the increment makes ``R_0``
tiny, so the step is held to a far stricter target than a cold start and
runs to the absolute floor -- the linear predictor's 256-step h40 run took
1,617 Newton iterations against 1,595 cold. A 24-step study on the same
deck with the relative test off (`results/cp_study/tol_study/`):

| absolute tol (MPa um^2) | cold iterations | linear-predictor iterations |
|---|---|---|
| 1e-6 | 105 | 92 |
| 1e-4 | 85 | 70 |
| 1e-2 | 66 | 47 |

Protocol for every run of the regenerated study (FE, hybrid, PINN):
``newton_rel_tol = 0``, ``newton_tol = 1e-3`` MPa um^2 (about 3e-7 of the
nodal force scale at h40, 4e-6 at h160), written into every pipeline deck
through the template. The hybrid and PINN-as-solver first checks launched
before this decision keep the old criterion and are superseded by the runs
on the chosen mesh.

### 10.1 Calibration

Iteration 1 (`cal_gn1.i`, center = the previous study's set, five runs,
256/256 each, 0 cutbacks): center RMSE 48.2 MPa against the markers (mean
relative error -11%); perturbations rho0 x1.2 33.5, k1 x1.2 17.9, y_c x1.2
70.7, tau_c0 x1.5 39.9 MPa. Sensitivities dS/dln p at 2 / 7 / 17%: rho0 149
/ 120 / 87, k1 46 / 206 / 368, y_c -13 / -88 / -208, tau_c0 26 / 25 / 24.
Levenberg-Marquardt step -> center 2: rho0 47.99, k1 0.1988, y_c 0.00344,
tau_c0 5.29 (predicted RMSE 14.2 MPa). Iteration 2 launched (`cal_gn2.i`).

Iteration 2 (`cal_gn2.i`, five runs, 256/256, 0 cutbacks): center RMSE 14.4
MPa (eight markers from 0.1%), perturbations rho0 x1.2 22.6, k1 x1.2 50.9,
y_c x1.2 33.6, tau_c0 x1.5 15.7 MPa. Sensitivities dS/dln p at 2 / 7 / 17%:
rho0 160 / 127 / 92, k1 52 / 238 / 419, y_c -15 / -103 / -241, tau_c0 25 /
25 / 24. The next Levenberg-Marquardt step moves rho0 1.1%, k1 1.4%, y_c
0.6%, tau_c0 -5.4% and predicts no gain (14.2 vs 14.4 MPa): iteration
stopped.

**Calibrated set (h40): rho0 = 47.99 um^-2, k1 = 0.1988, y_c = 0.00344 um,
tau_c0 = 5.29 MPa.** Per-marker errors: 0.14% -28.0%, 0.57% +1.9%, 2.0%
-1.8%, 7.0% +0.1%, 11.0% -0.2%, 13.1% -1.0%, 15.1% -1.2%, 17.0% -0.2%.

### 10.2 Mesh sensitivity

Decks `mesh_h{20,80,160}.i` regenerated from `cal_center2.i` (field series
included); h40 is `cal/gn2/base`. Launched 2026-10-12 00:00 MDT.

| mesh | 2% | 7% | 17% | 20% | change from the coarser mesh (2 / 7 / 17 / 20%) | wall |
|---|---|---|---|---|---|---|
| h20 | 413.1 | 518.7 | 616.0 | 629.8 | | 5 min |
| h40 | 403.9 | 502.8 | 591.5 | 603.3 | -2.2 / -3.1 / -4.0 / -4.2% | 14 min |
| h80 | 400.2 | 495.4 | 578.6 | 589.3 | -0.9 / -1.5 / -2.2 / -2.3% | 105 min |
| h160 | 397.9 | 490.6 | 569.2 | 578.9 | -0.6 / -1.0 / -1.6 / -1.8% | ~4 h |

All 256/256 steps, zero cutbacks, field series written. RMSE against the
markers 20.1 / 14.4 / 18.2 / 22.6 MPa (the set was calibrated at h40).
Ratio of successive changes h80->h160 over h40->h80: 0.62 at 2%, 0.75 at
20% (apparent order 0.7 falling to 0.4): first-order-like at yield, slower
on the plastic branch (slip localization without an internal length).
Distance to h160: h40 +1.5 / +2.5 / +3.9 / +4.2%, h80 +0.6 / +1.0 / +1.7 /
+1.8% at 2 / 7 / 17 / 20%. Report §5; the author chooses the mesh.

### 10.5 Hybrid first check at h40 (old parameter set and criterion; framework check only)

`decks/cp_study/hybrid/h40_{cold,linear,network}.i`, PINNACLE a54715a,
`scripts/hybrid_compare.py` (`results/cp_study/hybrid/hybrid_h40_summary.csv`,
figures `figures/report/hybrid_h40_{cost,prediction_error}.png`). All three
curves coincide (the corrector decides the answer).

| run | Newton iterations | Newton time | predictor time | proposal error (median / max, relative to the increment) |
|---|---|---|---|---|
| cold FE | 1,595 | 1,539 s | - | - |
| linear extrapolation | 1,617 | 1,571 s | 0.005 s | 0.0003 / 0.07 |
| network u(x, lambda), data-trained, GPU | 1,896 | 1,813 s | 17,085 s (training) | 1.64 / 3.4 |

Two findings. (1) The data-trained network is a BAD load-step predictor as
designed: it fits the accepted fields to a 5e-6 loss but extrapolates
poorly in the load factor (its proposal is worse than the previous step),
so Newton pays for it, and refitting on all accepted steps makes the
training cost grow quadratically (3.9 s at step 3, 135 s at step 256).
Linear extrapolation is four orders of magnitude better at no cost; a
neural predictor has to beat THAT. (2) Under the default relative stop test
even the near-exact linear proposal saved no iteration; with the absolute
criterion of §10.0 it saves about one per step. The physics-trained
predictor (`mode = physics`) is the route that can be right by
construction; its first check (`h40_pinn_solver.i`, corrector = none) is
running.

### 10.6 PINN-as-solver first check at h40 (old set and criterion; framework check only)

`decks/cp_study/hybrid/h40_pinn_solver.i` (physics-trained predictor,
`corrector = none`, w_data = 0, soft Dirichlet weight 10), 256 steps, 6.5 h
of training on the GPU (median 98 s per step). The network minimises its
own scaled strong-form residual to 1e-5, yet the accepted field is not in
equilibrium: residual_norm / reaction_norm = 1.0 at every step (the
out-of-balance force equals the load), and the mean stress is meaningless
(-7,300 MPa at 2%, +14,000 MPa at 20% against 375 / 546 MPa cold). Reading:
the soft Dirichlet term lets the network's boundary values drift; the
executioner then pins the prescribed degrees of freedom, and the jump at
the loaded face puts the boundary elements far from the field the residual
was evaluated on. With no corrector and no data term the next step trains
from that field, so the error compounds. Conclusion of the check: the
route needs the Dirichlet values built into the network ansatz (a hard
distance-function transform, as `pinn_fem_hybrid bc_mode = hard` does) and
the residual evaluated on the field that is actually accepted; that is the
next framework item for this route, after the FE regeneration. Results kept
in `results/cp_study/hybrid/h40_pinn_solver/` (curve and history only).

### 10.3 Parameter sensitivity at h40 (v2)

`sens_h40.i` on the new set, 17/17 runs 256/256, 0 cutbacks, ~3 h for the
set at ten concurrent. Full table in the report §6 (Table 6). Yield 271 /
330 / 364 / 402 / 496 MPa for rho0 x0.5 / 0.8 / 1 / 1.25 / 2; 20% stress
411 / 524 / 603 / 705 / 1,015 for k1 and 801 / 666 / 603 / 544 / 437 for
y_c; tau_c0 13-30 MPa per e-fold, shape unchanged. RMSE ranking at x2: k1
259, y_c 103, rho0 83, tau_c0 21 MPa.

### 10.4 Grain instances at h40 (v2)

Instance 2 (seed 20261013, 31 grains): +3.5% mean / 4.4% max above the
calibrated structure beyond 0.5% strain, RMSE 20.4 MPa; instance 3 (seed
20261012, 30 grains): -1.4% / 1.6%, RMSE 17.1 MPa. Both 256/256, 0 cutbacks.
Report §7.

### 10.7 Neural-route checks at h40 (new set, study criterion)

Reference: the calibrated FE run (`cal/gn2/base`): 1,506 Newton iterations,
838 s of Newton time over 256 steps. Linear-extrapolation predictor
(`hybrid/h40_linear`): identical curve, proposal within 0.03% of the
increment (median; 6.5% max), **948 Newton iterations (-37%) and 667 s of
Newton time (-20%)** -- under the absolute criterion the warm start pays,
as the tolerance study predicted. This is the bar for the physics-trained
network. Physics predictor (Newton corrector) and PINN as solver: running.

First physics-trained attempts (2026-10-12, killed after 7.8 h): the
PINN+FEM run produced a non-finite proposal at step 7 that poisoned Newton
(NaN residual, singular factorization) for the rest of the run; the
PINN-as-solver run (corrector = none) went 86 steps without a warning, then
the accumulated error of its unverified fields blew up the state update at
step 89 (slip rates 1e28/s); nothing was on disk because the history was
written at the end. PINNACLE fixes (commit after a2f95d6): proposal guard
(`predictor_reject_ratio`), `corrector = auto` with `corrector_tol`, live
`solver_history.csv`, per-step progress log. Both checks relaunched: the
hybrid with `corrector = newton`, the solver route with `corrector = auto`
(tol 0.05 of the reaction), so a step the network cannot solve is corrected
instead of accumulating.

**Relaunch with the guard (2026-10-12 18:08):** the physics-trained network
never solves a plate step. From the first predicted step its scaled
residual loss stays at 1.2-1.4, its proposal is worse than the previous
field (relative error 1.3-3.3 of the increment), the guard rejects it
(residual/reaction 0.999) and both routes fall back to Newton: cold cost
plus 1.2 s of training per step. Diagnostic
`scripts/tools/strong_form_floor.py` (4 FE steps on the torch GPU, 6 s;
network regressed to the FE nodal field of the last step in 3,000 epochs,
11 s): the regressed field reproduces the FE displacement to 0.55% in L2,
yet its FE residual is 0.73 of the reaction (the FE field: 0.048 after the
state advance) and its strong-form loss is 5.8 (median |r| 0.9, max 11.6 on
the mu/L scale). **A smooth global network within 0.5% of the FE field in
displacement is nowhere near equilibrium: the equilibrium of a 30-grain
plate lives in the strain kinks at the grain boundaries, which the network
cannot represent, so the strong-form collocation target is unreachable by
this ansatz regardless of training.** The strong-form, single-network
PINN is therefore not a viable solver or predictor on this problem class;
what remains viable is a network trained on the DISCRETE weak-form (FE)
residual of its nodal values -- the loss Newton itself minimises -- or a
domain-decomposed (per-grain) ansatz. Both are design decisions for the
author.

Both strong-form runs were stopped after 13-14 predicted steps (every
proposal rejected, residual/reaction 0.999, physics loss 1.2-1.4,
corrector always Newton; a profile put 0.2 s per epoch in the twelve
autograd passes of the strong-form residual, i.e. about 2 min per step for
a prediction that is discarded). Their partial live histories are kept
(`results/cp_study/hybrid/h40_{physics,pinn_solver}/solver/solver_history.csv`).

**Field comparison, FE vs the linear-predictor hybrid** (`decks/report/neural_compare.i`,
`scripts/field_compare.py`): displacement, equivalent plastic strain and
sigma_xx identical element by element at 2 / 8 / 15 / 20% (relative L2
difference 0 to the printed precision), as the corrector guarantees;
figures `figures/report/neural_h40_curves.png`,
`neural_h40_fields_{2,8,15,20}pct.png`.

### 10.8 Network trained on the assembled (FE) residual at h40 (author's option 1)

Decision of 2026-10-12: train the network on the discrete weak-form
residual of its own nodal values (docs/08, option a). PINNACLE edecbd0,
ADR 0029: `fem/residual_bridge.py` (forward: the FE assembler at the frozen
state with the Dirichlet values overwritten; backward: the transpose of the
consistent tangent, so the gradient is exact for the discretisation at one
residual-and-tangent assembly per loss evaluation) and `load_step_predictor
mode = fem_residual` (the network corrects the linear extrapolation, output
layer started at zero so epoch 0 is the extrapolation itself; loss = squared
residual-to-reaction ratio, the test `corrector = auto` applies; early stop
at `residual_target`; `optimizer = adam | lbfgs`; the best iterate is
returned). Bridge gradient verified against central differences on every
free degree of freedom (`tests/unit/fem/test_residual_bridge.py`); on the
two-element 3-D example the verified solver route accepts every step at a
ratio of 3e-4 to 6e-4 with the stress within 1e-5 of the cold run.

**Timing and optimizer checks, 8 steps at h40 (to 0.6% strain, yield
onset), `decks/cp_study/hybrid/check_fem_h40*.i`, FE arm numpy with
`OMP_NUM_THREADS=8`, network torch on the GPU.** One residual-and-tangent
assembly costs 0.28 s (a Newton iteration with its direct solve 0.45 s);
the warm-up steps take 2 Newton iterations, the plastic steps 3. The linear
extrapolation's residual ratio is 0.004-0.005 while the plate is elastic
(steps 3-4) and 0.03-0.36 once grains start to yield (steps 5-8).

| optimizer, epochs per step | ratio before -> after training, steps 3-8 | assemblies per step | predictor time per step |
|---|---|---|---|
| Adam 2e-3, 40 then 20 | 0.0048 -> 0.0048; 0.19 -> 0.057; 0.047 -> 0.029; 0.12 -> 0.12; 0.36 -> 0.31; 0.27 -> 0.22 | 22 | 6.2 s |
| Adam 1e-2, 40 then 20 | diverges from step 4: the carried-over correction wrecks the baseline (ratio 0.7-1.3 before training, 0.79-0.99 after) | 22 | 6.2 s |
| L-BFGS, 3 then 2 (<= 40 inner iterations) | 0.0048 -> 0.0008; 0.0042 -> 0.0031; 0.034 -> 0.023; 0.14 -> 0.11; 0.34 -> 0.29; 0.26 -> 0.20 | 15-54 | 6-16 s |

With the exact weak-form gradient the network improves the extrapolation
by 20-30% in residual at 10-50 assemblies per step, against 3 Newton
iterations (1.3 s) that solve the step. Whether this is the optimizer
budget or the ansatz is the capacity check below (L-BFGS with up to 400
inner iterations per step, with and without Fourier features).

**Capacity check (`check_fem_h40_lbfgs20{,_fourier}.i`, L-BFGS, 20 epochs
= up to 400 inner iterations per step):**

| step (strain) | ratio before | plain network, after (assemblies) | Fourier features (freq 4), after (assemblies) |
|---|---|---|---|
| 3 (0.23%) | 0.0048 | 0.0008 (32) | 0.0006 (84) |
| 4 (0.31%) | 0.0042 | 0.0029 (75) | 0.0022 (136) |
| 5 (0.39%) | 0.034 | 0.034 (122, no improvement) | 0.019 (360) |
| 6 (0.47%) | 0.16 | 0.11 (235) | 0.081 (453) |
| 7 (0.55%) | 0.32 / 0.44 | 0.29 (187) | 0.44 (182, no improvement) |
| 8 (0.63%) | 0.32 / 0.29 | 0.20 (168) | 0.14 (453) |

Up to 453 assemblies (2 min on the FE arm) per step do not bring a plastic
step below the 0.05 tolerance; the best plastic ratios are 0.08-0.29 and
two steps do not improve at all. **The limit is the ansatz, not the
optimizer budget or the loss:** with the exact weak-form gradient, the
global smooth network of (x, lambda) cannot represent the correction to
the linear extrapolation once grains start to yield, for the same reason
the strong-form check failed (the correction is as kinky at the grain
boundaries as the field itself). Full runs to 20% strain
(`h40_fem_residual.i`: L-BFGS 2 epochs per step, corrector = newton;
`h40_fem_solver.i`: L-BFGS 4 epochs per step, corrector = auto,
residual_target 0.05) launched 2026-10-12 to measure the two routes
end to end: results below.

**Full runs to 20% (2026-10-12, 1.0 h and 1.5 h on the host; report §8,
Tables 8-9, Figs. 14-16).** PINN+FEM with the FE-residual-trained network:
945 Newton iterations against 948 with the linear extrapolation alone
(251 steps identical), fields identical to FE, 1,574 training assemblies
(median 4 per step: L-BFGS finds nothing to take from a proposal already
at r = 1e-3), 349 s of predictor time added. Verified PINN as solver
(r_tol = 0.05, L-BFGS 4 epochs per step): 132 of 256 steps solved without
Newton, evenly from 1% to 20% strain, 650 Newton iterations in all; curve
within 5.1 MPa (0.85%) of FE, fields within 0.1% (u), 0.3% (plastic
strain), 0.6% (sigma_xx) in relative L2 at 2/8/15/20%; cost 16,264
assemblies, 4,964 s of training, 5,228 s total against 838 s cold FE.
Finding: the FE-residual loss makes the network's output verifiable and
lets it solve half the steps, but on a polycrystal a global smooth network
of (x, lambda) adds nothing to the linear extrapolation as a predictor and
costs 6-8x the FE run as a solver; a better predictor needs the grain
structure as input (docs/08).

### 10.9 Grain-aware network inputs (author's request, 2026-10-12)

PINNACLE 9ae5ff3 (ADR 0029 amendment): `pinn/node_features.py` and
`load_step_predictor node_features = 'orientation grain_boundary'` (per
node: the nine crystal-to-sample rotation entries averaged over the node's
elements, and one minus the largest fraction of the node's elements in
one grain). Decks `check_fem_h40_lbfgs20_grain.i`, `h40_fem_residual_grain.i`,
`h40_fem_solver_grain.i`: identical to the coordinate-only decks but for
the one key.

Capacity check (8 steps, L-BFGS 20 epochs): ratio after training at steps
3-8 is 0.0004 / 0.0013 / 0.009 / 0.043 / 0.099 / 0.076 against 0.0008 /
0.0029 / 0.034 / 0.108 / 0.287 / 0.198 for the coordinate-only network at
the same budget (440-460 assemblies per plastic step): 3-4x lower, below
the 0.05 tolerance at step 6. The representational limit of §10.8 is
lifted by the inputs, as expected.

Full runs (1.3 h each, two concurrent plus the check for the first 15
min): hybrid 923 Newton iterations (945 coordinate-only, 948 linear),
fields identical, 716 s of predictor time (2,696 assemblies); verified
solver 131/256 accepted, 635 iterations, 4,279 s of training (13,865
assemblies, 16 s per step median), curve within 0.33 MPa (0.05%) of FE
against 5.1 MPa coordinate-only, fields within 0.3% / 0.5% (plastic
strain / sigma_xx) at the worst strain and 0.01% / 0.02% at 20%. Report
§8 (Tables 8-9, Figs. 14-16). Finding: grain-aware inputs make the
network a far better corrector at the same tolerance, but the training
cost per step (10-20 s) still exceeds the Newton solve it replaces (1-3 s);
the solver route is correct and accurate at 5-6x the FE cost, the hybrid
route saves 2.5% of the iterations at a net wall-time loss.

### 10.10 Parameter sensitivity of the neural routes carried forward

Author's decision 2026-10-12: carry the linear-extrapolation hybrid and
the verified grain-aware solver. Decks `decks/cp_study/sens_h40_linear.i`
and `sens_h40_grain_solver.i` run the sixteen one-at-a-time cases of
`sens_h40.i` (x0.5, x0.8, x1.25, x2 of each calibrated parameter) through
`hybrid/h40_linear.i` and `hybrid/h40_fem_solver_grain.i`; analysis
`scripts/neural_sensitivity.py decks/report/neural_sensitivity.i`
(`results/cp_study/sens_neural/sens_neural_h40.csv`,
`figures/report/sens_neural_h40_cases.png`).

**PINN+FEM, linear extrapolation (16/16 done 2026-10-13 00:30, four
concurrent at two threads):** every case reproduces its FE curve to the
solver tolerance (largest stress difference 3e-6 MPa over the sixteen
cases, whose FE effects at 20% range from -32% to +68%); the Newton
iteration count is 941-974 per case against 1,468-1,536 for the FE cases,
a saving of 37% (median) that does not depend on the parameter. The
predictor's benefit is a property of the load path, not of the material
constants.

**PINN as solver, grain-aware, verified (16/16 done 13:21 2026-10-13,
two concurrent at four threads):** every case within 0.3-1.3 MPa
(<= 0.2%) of its FE case, 129-138 steps of 256 solved without Newton,
611-647 Newton iterations, 4,370-4,780 s per run; no trend with the
parameter (report §9.2, Table 11, Fig. 17).


### 10.11 Hardware and time benchmark (report §10)

`decks/benchmark/hardware_h40.i`, `scripts/hardware_benchmark.py`, 16 steps
per case, one case at a time (2026-10-13 13:30-13:55; the first pass had
to be repeated for a parsing bug in the script). Findings: numpy FE arm
0.45 s per Newton iteration at 1, 2, 4 and 8 threads (no threading in the
CP assembly); torch CPU 0.067 s (6.7x faster), torch GPU 0.079 s (5.7x),
at 7-9x the memory; linear predictor 37 vs 46 iterations at no cost;
grain-aware verified solver 8/16 steps without Newton, 12.6-12.9 s of
training per step whether the network is on the GPU or the CPU (the cost
is the 42 numpy assemblies per step). GPU memory not captured on this
host. Implication: all absolute route times of §10.7-10.10 (numpy arm)
are about 6x what the torch CPU arm would give; ratios unchanged.

### 10.12 Plate surrogate (PINNACLE ADR 0030) on h40

Design `decks/cp_study/surrogate/design_h40.i`: 24 LHS rows over rho0, k1,
y_c, tau_c0 (x0.5..x2, log), run by `scripts/surrogate_design.py` through
the study pipeline (4 concurrent, 2 threads: 623 s per sample, 14,961 s in
all), collected into the executioner's layout (curve = eng_stress vs
nominal_strain; final displacement field). Surrogate decks
`plate_h40_surrogate_{curve,curve_cpu,field}.i` (data source, 18/6
split, 3x64 tanh, Adam 4000 + L-BFGS 50, `evaluate` = 1,000 random rows).
Results (report §11, Table 10, Figs. 18-19): curve surrogate held-out
3.3% mean / 6.6% max on the GPU (2.0 / 5.2% on the CPU, same seed), field
surrogate 0.97% / 1.04%; training 28-121 s; 0.18-0.38 ms per curve query and 0.60 ms per
field query (the network is too small for the GPU to pay); break-even 14 FE runs (3
verified-solver runs). Accuracy is data-limited (sample 21, the lowest
hardening, is the worst case). Comparison `decks/report/surrogate_compare.i`,
`results/cp_study/surrogate/surrogate_h40_summary.csv`.
