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
