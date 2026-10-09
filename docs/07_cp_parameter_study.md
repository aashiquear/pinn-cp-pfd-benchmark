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

(in progress: h20 twin running; instrumentation of the failing point)
