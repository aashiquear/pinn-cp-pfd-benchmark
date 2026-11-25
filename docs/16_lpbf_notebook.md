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
