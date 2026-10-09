# Results: phase-field polycrystal plate and the first full CPFE run

Lab-notebook record. Workflow in `workflow_polycrystal_plate.md`, commands in
`reproduce_polycrystal_plate.md`. All runs on the ROCm host, numpy backend,
float64, CPU. PINNACLE `feat/pinn_benchmark` at `141ebd8` unless stated.

## 1. The phase-field grain-growth route works; the earlier stall was the time scheme

Earlier work (commits 23b1f69, 9f4e014) found the Allen-Cahn mobility
"inert" in the 8-order-parameter polycrystal configuration and parked the
route. That configuration ran `scheme = convex_splitting` at dt = 1.0 t_int,
i.e. dt = 8/(L mu). The convex-splitting stabiliser PINNACLE derives for the
Moelans energy over its declared band [-1.5, 1.5] with N = 8 components is

    S = mu (3 B^2 - 1 + 6 gamma (N-1) B^2) = 147.5 mu     (B = 1.5, gamma = 1.5),

so the semi-implicit update `(eta - eta_n)/dt = -L f'(eta_n) - L S (eta - eta_n) + ...`
has lag factor dt L S = 1180 at the old step and 37 at dt = 0.25/(L mu).
In that regime each step moves by about -f'/S whatever dt or L is -- which is
exactly "mobility does nothing".

Verification set, 16-grain seeding of `results/graingrowth/n16/deck.i`
(grid 158 x 79, w_GB = 1.904 um, M_GB sigma = 1, 1/(L mu) = 0.453), horizon
one interface time T = 3.624 in every case (`decks/graingrowth/verify/`,
outputs in `results/graingrowth_verify/`):

| case | scheme | M_GB | dt | steps | energy at T | drop from 28063.0 | gb_area |
|---|---|---|---|---|---|---|---|
| D | convex_splitting | 1 | 3.624 | 1 | 28062.96 | 0.04 | 404.8 |
| A | convex_splitting | 1 | 0.1133 | 32 | 28033.50 | 29.5 | 399.1 |
| B | convex_splitting | 1 | 0.2265 | 16 | 28046.31 | 16.7 | 401.1 |
| C | convex_splitting | 2 | 0.1133 | 16 | 28046.31 | 16.7 | 401.1 |
| E | bdf1 | 1 | 0.1133 | 32 | 27955.99 | 107.0 | 376.9 |
| F | bdf1 | 1 | 0.2265 | 16 | 27955.73 | 107.3 | 376.7 |
| G | bdf1 | 2 | 0.1133 | 16 | 27955.73 | 107.3 | 376.7 |

Readings. (i) Under convex splitting the per-step drop is 0.92 (A) and 1.04
(B) -- nearly independent of dt -- and B equals C to eight digits: the scheme
sees L dt but is saturated. (ii) Under bdf1, E and F agree to 0.2% (dt
convergence) and F equals G to eight digits (mobility connected, result a
function of M_GB t only). (iii) bdf1 delivers 3.6x the relaxation of
convex splitting over the same physical time, and 2,700x that of the old
one-step-per-interface-time setting. Cost: 18 s/step at 12.7k nodes x 8
components (direct sparse LU), 4-5 Newton iterations, no cutbacks.

Consequence for PINNACLE (to be filed): the `convex_splitting` scheme is
energy-stable but, with a Gershgorin bound of order 100 mu for multi-OP
energies, it is a damped fixed-point iteration whose reported time is not the
physical time. The executioner should report dt L S and warn when it exceeds
order one. Not changed in this session.

Also noted, not pursued: the "1.8x kinetics deficit" of commit edb46a3 was
measured on a deck whose `ac_bulk` was bound to a `double_well` material
(mobility 1.0, barrier 0.1) while `ac_interface` was bound to `gb_evolution`
(L = 26.7, kappa = 0.0019): two different mobilities and two different bulk
scales in one equation. The bundled PINN example
`grain_shrink_law_2d.i` has the same construction. That measurement does
not bear on the polycrystal route, where both kernels bind `gb`.

## 2. Pilot: 32-grain seeding, bdf1, three horizons

Seeds of `meshes/plate/poly_d3_n32.npz` on the 100 x 50 um plate; grid
203 x 102 (h = 0.493 um), w_GB = 1.485 um, 8 order parameters, M_GB = sigma
= 1, bdf1 at dt = 0.1377 = 0.5/(L mu). Decks `decks/graingrowth/pilot/`.

| horizon | steps | t / t_int | wall | free energy start -> end | gb_area (|grad eta| integral) |
|---|---|---|---|---|---|
| h40 | 40 | 2.5 | 28 min | 36070.8 -> 35958.2 | 592.4 |
| h80 | 80 | 5 | 56 min | 36070.8 -> 35936.1 | 570.2 |
| h160 | 160 | 10 | 1 h 52 min | 36070.8 -> 35893.5 | 527.6 |

Morphology of the labelled structure (`scripts/grain_morphology.py`,
interior grains only for shape metrics; junction statistics over all
junctions), against the Voronoi seeding sampled on the same grid
(`results/graingrowth_verify/voronoi_n32_baseline_morphology.json`):

| structure | grains | junctions | rms deviation from 120 deg | within 15 deg | circularity | boundary length (um) |
|---|---|---|---|---|---|---|
| Voronoi seeding | 32 | 45 | 30.3 | 35% | 0.715 | 648.6 |
| h40 (2.5 t_int) | 30 | 42 | 14.3 | 69% | 0.835 | 600.7 |
| h80 (5 t_int) | 30 | 42 | 14.1 | 71% | 0.862 | 579.8 |
| h160 (10 t_int) | 27 | 36 | 15.3 | 61% | 0.879 | 536.7 |

Figures `figures/polycrystal_plate_grains_h40.png`,
`figures/polycrystal_plate_grains_h40_h80.png` and
`figures/polycrystal_plate_grains_horizons.png` (colours are per panel, not
matched by identity). Doubling the horizon to 5 t_int keeps 30 grains,
shortens the boundary network a further 3.5% and rounds the grains
(circularity 0.835 -> 0.862) while the junction spread is unchanged at 14
degrees: by 2.5 t_int the junctions have equilibrated and what continues is
slow curvature-driven migration, including one grain shrinking toward
extinction (smallest area 21 -> 5.6 um^2). At 10 t_int three more grains
have gone (27 left, median equivalent diameter 14.3 um), the boundary
network is 17% shorter than the seeding's and circularity is 0.879; the
junction spread (15.3 deg) does not improve further because the grains now
meeting the walls and the newly formed junctions around the vanished grains
are still relaxing. Horizon is therefore a choice of grain count, not of
junction quality: 2.5 t_int already delivers equilibrated junctions. Two small grains were consumed; the boundaries
are curved and the junctions have moved toward 120 degrees. The remaining
spread of 14.3 degrees includes the grid's own angle-measurement floor
(about 4-6 degrees at this resolution, measured earlier on a hexagonal
tiling) and the plate-edge junctions, which the no-flux boundary pins at
90 degrees to the wall.

The h40 structure (30 grains, median equivalent diameter 13.7 um, i.e. about
7 mechanics elements across a grain) is the one the first CPFE run uses.

## 3. The MOOSE reference run

`moose/CPFE_out.e`: 1,250 HEX8, 2,652 nodes, single crystal Euler (0,0,0),
pulled in y at 0.15 um/s; 40 steps of 0.25 s written (the file preallocates
600), reaching 3.0% nominal strain. Mean sigma_yy: 47.8 MPa at 0.075%, 611.7
at 0.975%, 744.9 at 1.575%, 808.8 at 3.0% (`data/reference/moose_cpfe_plain_plate_stress_strain.csv`).
The deck header mentions phase-field damage but the deck has no damage
variable or kernel: it is a CP-only run. Its y_c = 0.00286 um differs from
the paper's table value 0.00268 (audit §G3); the deck's value is used here as
the reference.

## 4. Pipeline smoke test (16-grain one-step field, 2 load steps)

`scripts/polycrystal_plate.py` label -> orient -> deck on the one-step n16
field, `pinnacle validate` and a 2-step run to 0.15% strain: 10 s wall.
Before the translator fix (`141ebd8`) the run reproduced the MOOSE
single-crystal stresses 47.79 / 95.73 MPa to five digits, which exposed that
`[OrientationField]` was being dropped by the deck translator. After the fix
the 16-grain plate gives 54.13 / 108.44 MPa, sigma_xx = 0 and sigma_xy =
-0.06 / -0.11 MPa: a stiffer elastic response than the (0,0,0) crystal, as
random orientations of a cubic crystal with E_<001> at the soft end should.

## 5. Full CPFE run on the 30-grain plate (reference CP set)

`results/polycrystal_plate/n32_h40/cpfe.i`, generated from
`decks/polycrystal_plate/plate_n32_h40.i`: 50 x 25 x 1 HEX8 on 100 x 50 x 1
um, orientations from `orientations_points.csv` (30 grains, seed 20261009,
uniform on SO(3)), the MOOSE deck's parameters, 3.0% at 3e-3 /s in 40 steps.
PINNACLE `141ebd8`, numpy, 8 threads: **10 min 22 s wall, 40 of 40 steps
converged, zero cutbacks, at most 4 Newton iterations per step.** One
warning: the per-point strain increment reached 8.3e-4 against the G20
limit of 8.0e-4 (the mean increment is 7.5e-4), so the next run should use
48 or more steps. Curve in `results.csv`; figure
`figures/polycrystal_plate_curve_n32_h40.png`.

| nominal strain | PINNACLE 30-grain plate, mean sigma_yy (MPa) | MOOSE single crystal (0,0,0) (MPa) | experiment (digitized) |
|---|---|---|---|
| 0.075% | 52.5 | 47.8 | -- |
| 0.14% | -- | -- | 137.7 |
| 0.56% | -- | -- | 344.9 |
| 0.975% | 686.0 | 611.7 | -- |
| 1.575% | 865.2 | 744.9 | -- |
| 2.02% | -- | -- | 411.9 |
| 3.0% | 961.5 | 808.8 | -- |

sigma_xx and sigma_zz stay at zero to the printed precision, as the free
lateral faces require.

Readings, all tier 3 or weaker (no experiment is matched here):

1. The pipeline delivers what the plan asked for: a phase-field-grown
   polycrystal plate that runs a full crystal-plasticity solve from a deck
   plus two CSV companions, with the stress-strain curve written by the
   executioner.
2. With the published/reference parameter set, both codes sit far above the
   experiment (about 2.2x at 2%), and PINNACLE's polycrystal is 19% above
   MOOSE's single crystal at 3%. The direction is expected for a random
   polycrystal against a (0,0,0) crystal pulled along a cube axis (Taylor
   factor ~3.07 against 2.45 for <001> multi-slip), so the two numbers are
   not in contradiction, but they are also not a code-to-code comparison:
   different microstructure, and a different Taylor-term convention (C.15:
   PINNACLE forms tau_c from the latent-hardened forest sum, giving
   tau_c(0) = 176 MPa for this set against 43.7 MPa for the self term). A
   like-for-like check needs the single-crystal (0,0,0) plate in PINNACLE,
   which the same template runs by pointing the point file at one grain.
3. The elastic slope of the plate is 70 GPa (polycrystal) against 64 GPa for
   the MOOSE <001> crystal, consistent with Al's E_<001> being the soft
   direction. The digitized experiment's first two markers imply a steeper
   initial slope (about 99 GPa) than either, which is a digitization or
   test-compliance question for the calibration stage, not a model one.
4. Calibration is the next step, as the plan says: the CP parameters are
   model- and framework-sensitive and the reference set is a starting point.

## 6. Backend measurement on the host GPU (2026-10-09)

Host: Radeon 8060S (gfx1151, 20 CUs, unified memory), ROCm 7.2.4. The torch
wheel must match the runtime: `torch==2.13.0+rocm7.2` works; the rocm7.1
wheel segfaults on the first device tensor. The user's August probe records
the APU's float64 dense throughput at 0.37 TFLOPs on the GPU against 0.375
on the CPU (fp32: 2.36 against 1.05).

Same 4-step CP plate deck (h40 mesh, 800 elements, 30 grains), identical
stresses to eight digits on all three:

| backend | wall | max RSS |
|---|---|---|
| numpy, CPU (8 threads) | 32 s | 0.2 GB |
| torch, CPU | 145 s | -- |
| torch, GPU | 549 s | 2.2 GB |

Reading: the CP field solve is dispatch- and transfer-bound. `CPFieldCallback`
keeps the quadrature-point state as a Python list of per-point dicts and
crosses to numpy (`backend.to_numpy`) at more than fifteen places per
evaluation, so a device backend pays a host-device round trip per point per
sub-step. ADR 0015 recorded the same structure on an 8-element bicrystal
(17 s numpy, unfinished after 50 min on torch). Making the CP path
GPU-resident -- state as stacked device tensors, the sub-stepper and tangent
batched on device, assembly without per-point marshalling -- is a framework
project with this deck as its benchmark; on this APU the ceiling for float64
is parity with the CPU, so the gain would be in batching, not in the device.
Until that lands, the parameter study runs on numpy, per the measured-policy
rule of the guidelines.

Profile of one load step on torch-GPU (`cProfile`, 130 s total): 108 s in
`CPFieldCallback.tangent`, which called `CrystalPlasticity._evaluate`
224,000 times -- 35 single-point evaluations per quadrature point (6,400
points) through `_stress_numpy`, each a handful of 3x3 `einsum` / `inv` /
`det` launches (1.1 million `einsum` calls, 20 s of launch time alone;
`torch.as_tensor` host crossings 11 s). The batched finite-difference
tangent that commit 7156159 added is not the path a torch backend takes.
Scoped fix for the GPU track: route `tangent` and `stress` through the
batched constitutive path on every backend, with the point state held as
stacked device tensors rather than a list of per-point dicts; expected to
remove the per-point launch cost, which is the whole 17x. Not started.

Thread scaling (same 4-step numpy deck, `OMP_NUM_THREADS` = 1 / 4 / 16):
32.2 / 32.5 / 32.4 s at 99-108% CPU. A run uses one physical core whatever
the thread count, because the per-point constitutive loop is Python-bound
and the sparse CG solve is a small fraction. The host has 16 physical cores
(32 threads), so the practical policy is one core per run and up to ~14
runs concurrently; more cores per run buy nothing.
