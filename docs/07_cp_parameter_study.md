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

This is a derived starting point for the calibration in the author's order
(rho0, then k_m), not a fit: run `cp-h40-rescaled`
(`decks/cp_study/h40_rescaled.i`).

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
