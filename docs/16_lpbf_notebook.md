# 16. Part V notebook: laser powder bed fusion process benchmark

Lab notebook of Part V. Plan: `docs/15_lpbf_process_modeling_plan.md`. Final results go to
`docs/reports/lpbf_process_pinnacle_benchmark_report.md`; every run is a row of
`docs/run_ledger.md`. PINNACLE branch `feat/multi-gpu-and-am-process`.

## 1. Author's decisions (2026-11-25)

D1 NIST AM Bench alloy 718 is the benchmark; the 316L single track of
(grilli 2022 elimination) is shipped too, with coarse results and a full run, as its own
report section. D2 flat-surface flow scope. D3 conduction-driven chain first, then flow,
reported separately, one standalone deck per case. D4 orientations: plate grains random or
from the measured plate slab, nuclei random. D5 heat source fitted on two tracks, five held
out. D6 multi-GPU for the whole framework, NVIDIA CUDA and AMD ROCm (room for Apple Metal).
D7 crystal-plasticity parameters: literature set (prediction), then fit (calibration).
D8 316L kept. D9 NIST data are shown with citation.

## 2. Alloy 718 thermal data and their provenance

Read on 2026-11-25 from the PDFs in `/mnt/extended/projects/ambench-data/literature`.

| quantity | value | provenance | source |
|---|---|---|---|
| solidus, liquidus | 1255 C, 1337 C (1528.15 K, 1610.15 K) | `reported_text` (DTA) | (hosaeus 2001 liquid), section 3 |
| density at room temperature | 8190 kg/m3 | `reported_table` (their Table 2 caption, from the supplier) | (hosaeus 2001 liquid); (agazhanov 2019 thermophysical) measured 8170 +/- 80 |
| specific heat, solid | 425, 447, 468, 489, 510, 531, 635, 635, 635, 634 J/(kg K) at 298, 400, 500, 600, 700, 800, 1100, 1200, 1300, 1400 K | `reported_table` | (agazhanov 2019 thermophysical), Table 2 |
| specific heat at the solidus | 659 J/(kg K) | `reported_table` (DSC, extrapolated value in parentheses in the source) | (hosaeus 2001 liquid), Table 2 |
| specific heat, liquid | 775 J/(kg K) | `reported_text` (slope of the liquid enthalpy) | (hosaeus 2001 liquid) |
| conductivity, solid | 9.94, 11.59, 13.24, 14.91, 16.61, 18.34, 22.72, 23.61, 24.47, 25.32 W/(m K) at the same temperatures | `reported_table` | (agazhanov 2019 thermophysical), Table 2 |
| conductivity, liquid | 26.80 W/(m K) at 1337 C rising to 33.13 at 1800 C | `derived`: Wiedemann-Franz, L = 2.45e-8 V2/K2, from the volume-corrected resistivity of Table 2 (1.472 to 1.533 uOhm m) | (hosaeus 2001 liquid), Table 2 and section 4 |
| heat of fusion | 196.85 kJ/kg | `reported_table` (DSC) | (agazhanov 2019 thermophysical), Table 1 |
| heat of fusion, alternative | 290 kJ/kg between solidus and liquidus (enthalpy 679 to 969 kJ/kg, pulse heating), which contains the sensible heat of the 82 K range | `derived` | (hosaeus 2001 liquid), Table 2 |

Model choices: the mesh is fixed, so the volumetric capacity is `rho0 c_p(T)` with the
room-temperature density (mass per reference volume is constant); the solid conductivity
is held at its last tabulated value (1400 K) up to the solidus; the smaller, DSC heat of
fusion is used and the larger one is a sensitivity case.

Not yet sourced (`assumed`, values of the author's toolkit): boiling temperature 3186 K,
latent heat of vaporization 6.4e6 J/kg, molar mass 0.0588 kg/mol, recondensed fraction 0.18,
emissivity 0.4, film coefficient 10 W/(m2 K), plate temperature 296.65 K for the bare-plate
tracks. Their effect is checked before the calibration is reported.

## 3. Verification of the thermal components (PINNACLE 355b8d6)

| case | result |
|---|---|
| two-phase Stefan front, implicit | front +0.7 %, temperature within 0.004 of the exact profile (range -0.5 to 1) |
| the same, explicit | front +0.9 % |
| one track against Rosenthal, h = beam radius | 5.6 % rms, 17 % at most, between 0.3 and 0.8 from the beam |
| energy supplied by the beam | recovered to 0.02 % (numpy), 0.02 % (Torch, 64,000 elements) |
| Torch against numpy, same mesh | equal to 1e-9 |

Accuracy tier 2 (analytical solutions).

## 4. Case K1, first check (one case, source not calibrated)

Deck `decks/lpbf/in718/k1/k1_case0_h10.i`: baseline track (285 W, 960 mm/s, D4sigma 67 um),
10 um elements, conduction only, absorptivity 0.70 and penetration depth 52 um (the values
the author's toolkit fitted with its own code; here they are a starting point).

| isotherm of the peak temperature | depth (um) | width (um) |
|---|---|---|
| liquidus | 120.0 | 121.7 |
| solidus | 124.5 | 128.6 |
| measured (weaver 2022 trackdata), mean of six sections | 139.7 | 136.3 |

Depth 11 to 14 % low, width 6 to 11 % low, before any calibration (tier 4 comparison of an
uncalibrated model). Peak temperature 5,666 K, far above boiling: the evaporative sink
spread over a 10 um band does not hold the surface temperature down at this resolution.
141 s on the workstation CPU (480 steps, 0.29 s per step, numpy).

## 5. Case K1 at 10 um: seven tracks, heat-source fit (conduction only)

Targets: `data/lpbf/in718/k1_single_tracks.csv`, derived by `scripts/lpbf_k1.py targets`
from the NIST workbook (mean and sample standard deviation of the six sections of the
three single tracks of a case). Study decks `k1_study_h10_start.i` (starting source) and
`k1_study_h10.i` (after the fit). Host GPU, 23 s per case.

| case | P (W) | v (mm/s) | D4sigma (um) | depth measured | depth start | depth fitted | width measured | width start | width fitted |
|---|---|---|---|---|---|---|---|---|---|
| 0 (fit) | 285 | 960 | 67 | 139.7 | 122.1 | 118.2 | 136.3 | 125.4 | 133.7 |
| 1.1 | 285 | 960 | 49 | 227.2 | 133.0 | 126.6 | 106.2 | 120.4 | 125.9 |
| 1.2 (fit) | 285 | 960 | 82 | 102.4 | 114.9 | 111.8 | 141.7 | 131.1 | 140.0 |
| 2.1 | 285 | 1200 | 67 | 109.7 | 113.0 | 109.9 | 112.9 | 113.8 | 121.5 |
| 2.2 | 285 | 800 | 67 | 176.5 | 129.8 | 124.3 | 156.1 | 137.3 | 143.7 |
| 3.1 | 325 | 960 | 67 | 166.1 | 129.4 | 123.7 | 134.3 | 133.6 | 141.4 |
| 3.2 | 245 | 960 | 67 | 116.9 | 114.4 | 111.2 | 129.4 | 118.6 | 124.7 |

Lengths in um; computed values are the mean of the liquidus and solidus isotherms of the
peak temperature. Start: absorptivity 0.70, penetration 52 um. Fitted: 0.7724, 45.01 um
(one Gauss-Newton step on the relative depth and width errors of cases 0 and 1.2).

- The fit moves the widths of the two fitted cases to within 2 % and leaves their depths
  at -15 % and +9 %. A conduction model with one absorptivity and one penetration depth
  cannot follow the measured change of depth with spot size (140 um at 67 um, 102 um at
  82 um, 227 um at 49 um): the computed depth changes by 15 um over that range, the
  measured one by 125 um. The fit was stopped after one step for that reason (predicted
  rms 9.2 % against 10.3 %).
- Over all seven cases the depth rms error is 21.5 % (start) and 23.4 % (fitted); the
  width rms error is 8.6 % in both. The three cases with the deepest measured pools
  (1.1, 2.2, 3.1: small spot, low speed, high power) are 22 to 44 % too shallow. Their
  depth-to-width ratios are 2.1, 1.1 and 1.2: keyhole pools, outside a model with a flat
  surface.
- Peak temperatures of 5,000 to 7,600 K: the evaporative sink, spread over a 10 um band,
  does not hold the surface near the boiling point. In the experiment the recoil pressure
  opens the keyhole instead.
- Tier 4 comparison. This is the result of the conduction-driven stage; the flow stage is
  reported separately (decision D3).

## 6. Case K2, build column: deck, cost and a void first check

Deck `decks/lpbf/in718/k2/k2_column_h10.i` (header lists every choice). Temperature and
30 order parameters are advanced together on the explicit path; 129,600 elements;
0.54 to 0.66 s per step on the host GPU in float64, i.e. about 20 min per layer.

Numerical parameters of the grain stage at 10 um: boundary width 40 um (4 elements);
`L(T_l) kappa = 3e-6 m2/s` (45 % of the stability limit `h^2 / (6 dt)`); melting rate
2e5 /s (`P dt = 0.5`); regrowth rate 1e5 /s, front speed `2 sqrt(D L kappa) = 1.1 m/s`,
above the scan speed, so the front follows the isotherm except at the tail of the pool;
melt switch 2 K below the liquidus. Assumed physical parameters: site density
3.7e13 /m3, critical undercooling 5 +/- 2 K, powder conductivity 0.3 W/(m K), plate at
353.15 K, far-field film coefficient 2e4 W/(m2 K). Between layers the block is reset to the
plate temperature (`layer_reset_temperature`).

**Void check (kept as a record of the defect).** The first one-layer run filled the
space of layers 2 to 4, which holds no material during layer 1, with grains: the regrowth
term acted wherever the temperature was below the liquidus. The grain summary reported
75 % of the "built" volume within 15 degrees of <100> along the build direction; the
section image showed the cause at once (the whole region above the surface in one color).
PINNACLE 0d7befe gates the regrowth by the material state (`process_state` on the
kernel); the column decks name it. No number of that run is used.

## 7. 316L case: what is reported, derived, assumed

From (grilli 2022 elimination), published version, Table 2 and sections 2 to 4
(`reported_table` / `reported_text`): domain 200 x 160 x 108 um; meshes of 4 um (54,000
elements) and 8 um (7,000); 150 W, 700 mm/s; flow rule `gamma_dot0 = 1e-7 /s`, `m = 0.1`;
hardening `h0 = 3,839 MPa`, exponent 2.5, resistance 213 MPa at 303 K, saturation 302 MPa;
`C11, C12, C44 = 204.6, 137.7, 126.2 GPa` with slopes -90.33, -45.10, -51.78 MPa/K;
resistance factor `k_A + k_B exp(-k_C (T - T0))`, 0.53, 0.47, 0.008; melting 1648.15 K;
ambient 298 K; volumetric expansion 44.73e-6 /K, slope 0.01011e-6 /K2; residual stiffness
0.01 and 0.001; temperature sampled every 45 us; time step 1 us (0.1 us around a
temperature update); no normal displacement on x = 0, x = 200 um, y = 0, y = 160 um and
the bottom; cost 9.5 h on 12 CPUs (8 um), 99 h on 16 CPUs (4 um).

`derived` for the decks: linear expansion 14.91e-6 /K and slope 3.37e-9 /K2 (a third of
the volumetric values); one stiffness factor with slope -3.93e-4 /K (mean relative slope
of the three constants); dislocation-density parameters that give the paper's initial
resistance, saturation resistance and initial hardening rate (deck header): the paper's
phenomenological hardening law exists in PINNACLE (`kalidindi`) only as a point model
without a batched field path.

`assumed`: beam radius 35 um, absorptivity 0.5, penetration 30 um (the paper refers to its
thermal-fluid papers); 316L thermal properties (liquidus 1723.15 K, heat capacity
3.95e6 to 6.32e6 J/(m3 K), conductivity 14 to 30 W/(m K), latent heat 2.13e9 J/m3);
no powder layer; shear modulus 75 GPa, Taylor coefficient 0.3, latent hardening 1.4;
844 Voronoi grains of 16 um with uniform random orientations (seed 1) in place of the
paper's phase-field grains with EBSD orientations; load steps of 45 us.

## 8. Device pool (PINNACLE ADR 0050)

`[Global] devices = 'cuda:0 cuda:1 ...'` (or `num_gpus`) splits the element loop of
`pf_fem_explicit` over the devices. On the host's one GPU, pools of two slots of the GPU
and of GPU + GPU + CPU reproduce the single-device field to 4e-15 (13.1 s single, 14.7 s
and 18.6 s pooled: the overhead of the transfers, no speed-up on one physical device).
Scaling on several GPUs is not measured. The mechanics solve and network training are not
on the pool.

## 9. Case K2: one layer and four layers at 10 um (PINNACLE 996d9d3)

Two further one-layer runs after the void check of section 6 are also void as grain
results: the grains of the plate coarsened (a mobility sized for the solidification front
moved the boundaries of the cold solid). `mobility_reference_temperature` with
`activation_energy = 1.5e6 J/mol` makes the mobility fall below the liquidus; the decks
give the mobility at the liquidus.

Valid runs (`results/lpbf/in718/k2/layer1`, `column_h10`; one seed each):

| quantity | one layer | four layers |
|---|---|---|
| wall time, steps | 1,226 s, 2,240 | 4,820.9 s, 9,040 (0.533 s per step) |
| nucleation events | 621 | 3,704 |
| grains of at least 8 elements (fragments below) | 42 (36) | 378 (340) |
| grains not touching the plate | 14 | 258 |
| equivalent diameter, median / mean (um) | 71.3 / 89.1 | 36.8 / 47.2 |
| height, median (um); aspect, median | 85; 0.58 | 40; 0.76 |
| consolidated fraction above the plate | 0.215 | 0.912 |
| build direction within 15 degrees of <100>, <110>, <111> | 0.076, 0.304, 0.329 | 0.192, 0.157, 0.195 |
| surface time above melting, p25 / p50 / p75 (ms) | p50 2.31 | 1.70 / 2.39 / 3.02 |
| cooling rate at freezing, surface, p25 / p50 / p75 (K/s) | p50 7.6e4 | 5.3e4 / 6.3e4 / 8.0e4 |
| gradient at freezing, all melted elements, p50 (K/m) | | 3.56e6 |

Measured: median diameter 11.5 um; texture 0.10, 0.21, 0.10 (random 0.10, 0.20, 0.14);
thermography of leg L9 (`data/lpbf/in718/k2_thermography_leg9_summary.csv`): time above
melting 1.28 ms, solid cooling rate 1.9e5 K/s. The cooling rates are different quantities
(NIST: solid, over a temperature window; here: at the moment of freezing).

Reading of the figure `figures/lpbf/k2_grains_column_h10.png`: plate grains continue
across the plate surface; the grains are blocky (a boundary is 4 elements wide); the
corners of the block at the height of layer 2 are not consolidated. Size and texture are
not results at this mesh.

## 10. 316L case: the relaxation stage was missing in the first mechanics run

First run (`results/lpbf/ss316l/mechanics_h8_walls_held`, 178 steps to 8.01 ms, 5,670 s,
439 Newton iterations, no cutback): mean stress of the volume xx 1,797 MPa, yy 1,986 MPa,
zz 0.8 MPa; von Mises mean 2,111, 95th percentile 3,152, maximum 3,686 MPa; equivalent
plastic strain mean 0.016, maximum 0.053. The paper's Figure 15 shows slice means between
-400 and +130 MPa and local values between about -900 and +1,500 MPa.

Cause, found by reading the paper's section 4.1 again: its result is taken after a
relaxation simulation in which the walls x = 200 um and y = 160 um are removed. The first
run ended with all four walls fixed, so its mean stress is the reaction of the walls on a
volume that was heated to at least 800 K everywhere and cooled. A second contribution is
expected from the flow rule (reference rate 1e-7 /s, exponent 0.1: at a strain rate of
1 /s the resolved stress is five times the slip resistance); it is not separated here.

PINNACLE: `release_time` on a Dirichlet condition (applied by `cp_quasi_static`, refused
elsewhere). Decks: walls xmax and ymax released at 8.01 ms, 22 further steps of 45 us at
the final temperature (the paper does not print the length of its relaxation stage;
assumed).

**Two-stage run** (`results/lpbf/ss316l/mechanics_h8`, PINNACLE 4827247, 200 steps,
2,957 s, 463 Newton iterations, no cutback). Step 178 (walls fixed) reproduces the first
run: mean xx 1,796.8 MPa, yy 1,986.1 MPa. Step 179 (walls released): mean xx 7.4, yy 3.9,
zz 0.6 MPa, three Newton iterations; steps 180 to 200: one iteration each, no change of
stress or plastic strain (elastic unloading; the assumed length of the relaxation stage
has no influence). After relaxation: von Mises mean 691, 95th percentile 1,177, maximum
1,887 MPa; element sigma_xx from -1,290 to 2,392 MPa; equivalent plastic strain mean
0.016, maximum 0.053.

Layer means against the height (z in um: sigma_xx, sigma_yy in MPa): 3.9: -429, -472;
11.6: -400, -461; 19.3: -332, -423; 27.0: -240, -361; 34.7: -132, -240; 42.4: -19, -22;
50.1: 116, 327; 57.9: 182, 422; 65.6: 235, 417; 73.3: 248, 368; 81.0: 241, 284;
88.7: 223, 188; 96.4: 228, 83; 104.1: 183, -53. Paper, Figure 15d, read by eye: about
-420 / -400 at the bottom, sign change near 30 / 22 um, maxima about +130 / +95, zero at
the top.

## 11. Case K1 at 5 um, case 0 (cost and mesh check, PINNACLE 4827247)

706 s on the host GPU for 588,000 elements and 1,920 steps (23 s at 10 um: the cost is
linear in elements times steps, 6.3e-7 s per element step). Depth 117.6 um and width
133.3 um (mean of the two isotherms) against 118.2 and 133.7 um at 10 um: the 10 um
result is converged in the mesh within 0.5 %. The depth error against the measurement
(-15 %) is therefore not a mesh error.

## 12. Flow stage: inputs, verification, first track (PINNACLE ADR 0052)

**Decisions of the author, 2026-11-25:** multi-GPU and flow next; coarse meshes only
(fine decks are shipped, the author runs them on 8 x H200 after the push); an animation
stack is kept for every run; three seeds for stochastic cases; the process-parameter
study is later; K3 starts from the toolkit's calibrated tensile CP set and gets its own
calibration study with PINNACLE.

**Flow inputs of alloy 718** (deck `k1_flow_template.i`, study `k1_flow_study_h10.i`):

| quantity | value | provenance |
|---|---|---|
| liquid density | 7,578 kg/m3 | `reported_table` of (amin 2024 heatsource), read through a literature extraction tool, not from the PDF; (nishi 2024 density) measured 7,430 |
| viscosity | 5.3 mPa s | the same table of (amin 2024 heatsource) |
| d gamma / dT | -1.1e-4 N/(m K) | `assumed` (value in common use; not sourced) |
| volumetric expansion of the liquid | 5e-5 /K | `assumed` (value of the author's toolkit) |
| drag constant A | 1e8 /s, eps 1e-3 | `assumed` (toolkit) |
| recoil coefficient | 0.54 with the evaporation constants of the heat loss | `assumed` |
| artificial sound speed | 25 m/s | numerical |
| upwind factor | 1 (first order) | numerical |

**Verification** (PINNACLE examples, tier 2 / 3): lid-driven cavity Re 100, 40 x 40:
-0.2148 / +0.1802 / -0.2540 against -0.2109 / +0.1753 / -0.2453 (ghia 1982 cavity);
side-heated cavity Ra 1e4: Nu 2.2433, u 16.21, v 19.68 against 2.243, 16.178, 19.617
(devahldavis 1983 natural); thermocapillary layer: interior within 0.007 of the closed
form, surface speed -0.212 against -0.25 (band).

**Two defects found on the first track and corrected** (kept as a record):

1. Drag as an explicit point term: non-finite fields after 7 steps (a liquid node next to
   solid quadrature points sees `dt D = 1e4`). Now a nodal backward-Euler step.
2. Heat advection as `C u . grad T`: the stored energy grew to 1.68 times the supplied
   energy by 1.2 ms and went on growing after the beam had left; the pool was 220 um wide
   and never froze. The discrete velocity is divergence free only weakly. Now the
   conservative form: stored / supplied = 0.917 from 450 us to the end (the rest is the
   surface loss). The run with the first form is void.

**Case 0 with flow, no recoil, conduction-calibrated source (one run, 507 s):** depth
116.1 um, width 162.7 um (conduction 118.2 / 133.7; measured 139.7 / 136.3); largest
speed 1.35 m/s; peak temperature 6,100 to 6,500 K. With recoil (trial of 150 us): largest
speed 3.3 m/s, peak temperature 4,870 K.

**Seven tracks with flow and recoil** (`results/lpbf/in718/k1/flow_h10`, 505 to 512 s
each): depth / width in um and error against the measurement: case 0 124.4 (-11.0 %) /
134.3 (-1.4 %); 1.1 131.1 (-42.3) / 116.0 (+9.2); 1.2 119.2 (+16.3) / 147.1 (+3.8);
2.1 111.5 (+1.6) / 122.7 (+8.7); 2.2 135.9 (-23.0) / 142.9 (-8.5); 3.1 134.0 (-19.3) /
140.6 (+4.7); 3.2 114.2 (-2.3) / 126.4 (-2.3); rms 21.0 % / 6.3 %. Peak temperatures
4,312 to 5,293 K. Case 0: stored / supplied energy 1.004 with the losses on; with the
losses off 1.009 with flow and 1.008 without (300 us runs): the flow adds 0.1 %.

## 13. Case K3: targets, starting set, calibration (2026-11-25)

Targets from the raw NIST records (`scripts/lpbf_k3_targets.py`,
`data/lpbf/in718/k3_tensile_mean_curve.csv`, `k3_tensile_summary.json`): 8 specimens,
gauge 3 mm; E 128.1 +/- 5.3 GPa (100 to 400 MPa fit; NIST 123.6), yield 687.3 +/- 3.3 MPa
(NIST 690.7), ultimate 951.6 +/- 9.2 MPa, strain at ultimate 0.246 +/- 0.020; mean curve
to 21.9 %.

Starting set from the toolkit's `calibrations/in718/tensile_cp.json` (Voce: g0 227.217,
gs 666.333, h0 1407.207 MPa, a 3.1017, q 1.4, gdot0 1e-3, m 0.05, C 250/160/120 GPa;
Taylor, 400 grains, rms 18.5 MPa there): tau_c0 224.72 MPa, k1 0.12866, y_c 3.6411 nm
(rho0 1e10, alpha0 0.3, mu 81 GPa, b 0.254 nm assumed). Only this tensile set of the
toolkit is used.

PINNACLE's 0-D inverse route (`cp_constitutive`, `[Calibration] mode = inverse`) was
tried with a Taylor polycrystal and a strain table: 1.8 GPa at 0.1 % strain, negative
stresses; not used, not investigated.

Gauss-Newton (log parameters, 5 % finite differences, damping 1e-3, 88 points from 1 %):
rms 66.58, 31.26, 3.59, 2.75, 2.662, 2.661 MPa. Calibrated: tau_c0 263.365 MPa,
k1 0.0526896, y_c 1.48908 nm. Jacobian condition 39 to 95.

Three seeds with the calibrated set: rms 2.66 / 7.16 / 8.46 MPa; yield 696.5 / 690.3 /
698.5; 5 %: 832.4 / 824.9 / 837.4; 10 %: 891.1 / 883.9 / 898.4; 20 %: 941.4 / 936.4 /
953.3 MPa; initial slope 219.6 / 217.9 / 219.7 GPa.

**316L hardening correction.** While deriving equation 18 the same derivation of the
316L deck was found wrong: its k1 = 4.0e-3 (y_c 0.52 nm) gives the saturation of the
paper's law but a hardening rate lower by S = 1 + 11 q = 16.4. Corrected to k1 = 0.065155,
y_c = 8.2581 nm; the earlier result is kept in `results/lpbf/ss316l/mechanics_h8_k1_low`
and the case is rerun.

## 14. Case K2 on three seeds; PINNACLE 037245e and f720f35 (2026-11-25)

Decks `decks/lpbf/in718/k2/k2_column_h10.i` (seed 1), `_seed2.i`, `_seed3.i`; results
`results/lpbf/in718/k2/column_h10`, `column_h10_seed2`, `column_h10_seed3` (not shipped).

| quantity | seed 1 | seed 2 | seed 3 |
|---|---|---|---|
| grains / fragments | 378 / 340 | 375 / 344 | 361 / 383 |
| not touching the plate | 258 | 249 | 244 |
| nucleation events | 3,704 | 3,516 | 3,640 |
| median / mean diameter (um) | 36.76 / 47.16 | 37.67 / 47.78 | 38.55 / 47.90 |
| consolidated fraction | 0.9125 | 0.9126 | 0.9121 |
| texture 100 / 110 / 111 | 0.192 / 0.157 / 0.195 | 0.134 / 0.110 / 0.082 | 0.288 / 0.267 / 0.115 |
| wall (s) | 4,821 | 4,855 | 4,831 |

The thermal summary is identical in the three runs. The three runs predate the
nucleation correction of PINNACLE f720f35 (`undercooling_nucleation process_state`); the
standalone column decks use it.

PINNACLE after the flow stage: 037245e (ADR 0053, crystal-plasticity points on the
device pool; identical on CPU slots, 2e-13 relative on repeated GPU slots; 16.3 / 17.9 /
26.1 s for one, two CPU and two GPU slots: no gain on one GPU) and f720f35 (example chain
`examples/applications/lpbf_process_chain`, four decks: 64, 197, 9, 441 s).

## 15. Third stage (2026-11-25/26): surface correction, defects, K4, EBSD cube, four GPUs

PINNACLE f4eefa4 (ADR 0054, 0055, documentation). Decks, scripts and results of every
run: `docs/run_ledger.md`.

- **Surface under unbuilt layers (ADR 0054).** Second run of the column
  (`k2_column_h10_r2.i`): 391 grains, median 37.7 um, time above melting 2.375 ms
  (before: 378, 36.8 um, 2.39 ms). The correction changes the 10 um result little.
- **Column with flow.** Check of 1.4 ms: finite, 2.65 m/s, 1.016 s per step; the pool
  reaches the sides of the 600 um block, so the deck has velocity walls on the cut sides.
- **Keyhole number.** Material values at the liquidus: rho c_p = 6.3473e6 J/(m3 K),
  k = 26.8 W/(m K); T_0 = 296.65 K (tracks), 353.15 K (build). The depth law
  e / r0 = 0.4 (Ke - 1.4) was recomputed on the 87 rows of Supplementary Data 1 of
  (gan 2021 keyhole): R2 = 0.977. The absorptivity there is measured (0.12 to 0.67); the
  values 0.50 and 0.35 used for alloy 718 are `assumed`.
- **Keyhole-number source depth.** Decay length = e / 2, smallest value 10 um, source
  absorptivity 0.7724 (`derived`, the conduction fit). Depth rms 41.8 % (0.7724 in Ke)
  and 21.1 % (0.50); fitted depth: 23.4 %.
- **Lack of fusion.** Hatch 110 / 150 / 190 / 230 um: unmelted share 0 / 0.0001 / 0.106 /
  0.294; overlap index crosses 1 between 150 and 190 um.
- **316L, three seeds.** Layer means within 26 MPa.
- **K3 on the measured structure.** Cube of NIST mds2-2767 (2 um volume, every second
  voxel, x, y from 400, 424 um of the volume corner, z 200 to 320 um above the plate):
  yield 638.9 MPa (measured 690.7).
- **K4.** `assumed`: expansion 13e-6 /K + 5e-9 /K2, resistance factor 1.0 / 0.85 / 0.60
  at 800 / 923 / 1,073 K and constant above, cut-off 1,273 K, stiffness slope -3.5e-4 /K.
  A resistance factor of 0.05 with the slip exponent 20 did not converge (factor 1e26 on
  the slip rate).
- **Four GPUs (Modal, 4 x A10G).** Identical fields; no useful speed-up (Table 17a, b of
  the report); 4.66 USD metered, covered by credits.
- **Datasets** for the defect indicators: `data/lpbf/defects/raw/*/SOURCE.md` (the data
  files are not tracked).
