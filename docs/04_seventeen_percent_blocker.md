# Reaching the experimental 17% strain

The Gc calibration is anchored at the experimental fracture strain (17.0%,
592.6 MPa, the last Vaughan 2024 marker, where the test was stopped by
catastrophic fracture). A model that cannot be loaded that far cannot be
calibrated against it, so this was resolved before any Gc work.

## Two independent defects, separated by bisection

**Rate mismatch.** The slip law is rate dependent (reference rate
gamma0 = 1e-5 /s, n = 8.33). The Taylor calibration integrated at 1e-3 /s;
the plate deck hard-coded `t_end = 10 s`, which on a 17% ramp is 1.7e-2 /s
-- 17x faster, a factor 1.4 in stress. `t_end` is now derived as
strain / strain_rate, with the rate an explicit argument defaulting to the
calibration value.

**Ramp decoupled from `t_end`.** The displacement function hard-coded
`t/10.0`, the ORIGINAL end time. Deriving `t_end` without changing the
function made the ramp t_end/10 times too steep, so the first load step
demanded the whole deformation. This was self-inflicted, introduced by the
rate fix above. The function now reads `t/@TEND@`, verified at two strain
levels to reach exactly the target strain at `t_end`.

## What the bisection showed

Identical mesh, rate, and parameters; the only difference is the
orientation field.

| configuration | cutbacks | strain reached | outcome |
|---|---|---|---|
| single crystal, plain plate | 0 | 17.00% | complete |
| 12-grain polycrystal, same plate | 180 | 10.20% | failed |

**Attribution, corrected.** An earlier write-up said the original 180
cutbacks were "mostly the broken ramp, and only partly grains". That is
wrong and the data in the same run contradicts it: the polycrystal produced
**180 cutbacks both before and after** the ramp fix -- an unchanged count.
The ramp bug was decisive for the single crystal (180 -> 0) and accounts for
none of the polycrystal's cutbacks. Grain interaction accounts for
essentially all of them.

## Consequence for the two calibration tracks

They are no longer blocked by the same thing:

* **Gc (single crystal [-1 4 9])** -- unblocked. The single crystal reaches
  17% with zero cutbacks in about 5 minutes at d=1.
* **CP parameters (polycrystal)** -- still limited to ~10.2%. The
  experimental curve is a polycrystal measurement, so this track needs the
  grain-interaction stepping resolved before it can be calibrated over the
  full range.

## A withdrawn number

An earlier report that the finite-element polycrystal was "3-6% softer than
experiment, the expected direction" is withdrawn. It came from a run at the
mismatched rate, so it is the near-cancellation of two large opposing errors
-- grain interaction softening the response and the 17x rate stiffening it
-- not a measurement of the Taylor-to-finite-element offset.

## A trap worth recording

`CrystalPlasticity` reads `config["material"]`. A FLAT config dict silently
takes SI-unit defaults (C11 = 108.2e9 Pa) in a um/MPa deck, producing
stresses around 1e10 MPa and rho around 1e18 with no error raised. Two
separate diagnostics in this study were invalidated this way before it was
recognized.

## The limit is stress concentration, not the polycrystal

Later runs broadened the picture. Every configuration containing a stress
concentrator stalls before 17%; the homogeneous one does not.

| configuration | initial cutbacks | reached 17%? |
|---|---|---|
| plain plate, single crystal | 0 | yes, verified to 0.22% against the material point |
| plain plate, 12 grains, d=1 (3 elem/grain) | 180 | no, 10.2% |
| plain plate, 12 grains, d=2 (~14 elem/grain) | 6 | no, stalled later |
| plain plate, 12 grains, d=3 (33 elem/grain) | 1 | no, stalled later |
| holed plate, single crystal [-1 4 9], d=4/5/6 | 779 / 163 / 202 | no |

Two readings follow, and only the first is established.

**Elements per grain matters, and the original specimen was pathological.**
At d=1, 69% of elements straddled a misorientation discontinuity, the
smallest grain was 3 elements, and only 2 of 12 grains were interior.
Raising resolution cut the initial cutback count from 180 to 1. An earlier
suggestion to add MORE grains is wrong: at fixed mesh it worsens the same
quantity (108 grains gives 2 elements/grain and 75% of elements on a
boundary).

**Resolution delays the wall but does not remove it.** The d=2 and d=3 runs
stalled at higher strain with cutbacks climbing again, and the holed SINGLE
crystal stalls with no grains at all. So the controlling factor is the
presence of a stress concentrator, not grain structure specifically.

## Mechanism, quantified, with what is NOT established

The slip law raises a local overstress to the power n = 8.33, so a
concentration is amplified enormously in the slip rate:

| local Kt | tau/tau_c | gamma_dot (1/s) | sub-steps needed |
|---|---|---|---|
| 1.0 (homogeneous) | 1.46 | 2.3e-4 | ~0 |
| 2.0 | 2.92 | 7.6e-2 | 57 |
| 3.0 (circular hole) | 4.38 | 2.2 | 1,674 |
| 3.5 | 5.11 | 8.0 | 6,050 |

Exhausting the 100,000 sub-step budget requires tau/tau_c = 7.15, i.e.
**4.9x** the homogeneous value. The elastic Kt = 3 of a circular hole gives
1,674 sub-steps, comfortably inside budget. **So the elastic concentration
alone does not explain the failures**; additional localization is involved,
either plastic strain concentrating beyond the elastic field or
misorientation across a grain boundary. Which one is NOT established, and
should not be asserted without measuring the local field at the failing
point.

## Integration tolerance is now explicit

`max_shear_increment` (1e-3) and `max_substeps` (100,000) were hard-coded in
`CPFieldCallback.update_state`. They set both the accuracy/cost trade-off and
whether a solve completes at all, so they are now constructor parameters with
unchanged defaults (pinnacle 06e6c17). **Whether the cap can be loosened is
unresolved**: a hand-rolled emulation of the sub-stepping produced 269 MPa
where the verified material point gives 313 MPa under identical conditions,
so it did not reproduce a validated result and cannot be used to judge the
cap. A faithful test must drive the real callback.
