# 15. Part V plan: laser powder bed fusion process modeling

Started 2026-11-25, revised the same day (benchmark changed to the NIST AM Bench alloy 718 data,
section 2). Status: plan for the author's review. No PINNACLE code has been written and
no simulation has been run for this part.

Aim (author, 2026-11-25): add laser powder bed fusion (LPBF) process modeling to PINNACLE and
benchmark it end to end against a published process chain: thermal field with melt-pool flow,
grain structure by phase field, crystal plasticity for the stress left by the process and for
the properties of the built material. Order of
work: this plan; then the capabilities in the local PINNACLE checkout; then a coarse-mesh
capability check on `rocm-server`; then a standalone full-accuracy deck that a user can run on
larger hardware.

Where things go:

| content | location |
|---|---|
| new capability (components, tests, examples, ADRs, docs) | `/mnt/extended/projects/pinnacle`, branch `feat/multi-gpu-and-am-process` (local first) |
| decks, scripts, figures, benchmark numbers, report | this repository, local first |
| meshes, if a generated box is not enough | `/mnt/extended/projects/hexaform` |

`PROJECT_GUIDELINES.md` section 3 names `feat/pinn_benchmark` as the framework branch. For this
part the author named `feat/multi-gpu-and-am-process`. On 2026-11-25 that branch is identical to
`master` at `ddf1e9e` (the merge of `feat/pinn_benchmark`); it holds no process-modeling or
multi-GPU code yet.

Sources of this plan: the literature note `docs/literature/ScholarGPT_LPBF_process_modeling.md`
(a navigation aid; it states that its references were not retrieved live), a literature search
of 2026-11-25 (45 references and data records added to `bib/references.bib`, each DOI checked
against Crossref or DataCite), a capability survey of PINNACLE at `ddf1e9e`, the local NIST
data in `/mnt/extended/projects/ambench-data` with its README and the NIST record READMEs, and
a reading of the thermal, flow, phase-field and crystal-plasticity stages of the author's LPBF
toolkit (github.com/aashiquear/AM-Orientation-Shift-TOOLKIT, JAX route and MOOSE route),
including its benchmark document for the same NIST data.

## 1. The process chain

The published chains have the same four links, coupled one way (review of the candidates in
section 2):

1. **Thermal field** of the scanning laser: conduction with latent heat, or conduction with
   melt-pool flow (Marangoni shear, recoil pressure, buoyancy, a mushy-zone drag), with or
   without a deforming free surface (khairallah 2016 meltflow).
2. **Grain structure**: a phase-field model (yang 2021 grain; chadwick 2021 grain;
   liu 2018 sitespecific) or a cellular automaton (gandin 1994 cellular; rolchigo 2022 exaca)
   driven by the transient temperature field: melting erases grains, the pool boundary regrows
   epitaxially, nuclei appear in the undercooled liquid.
3. **Mechanics during the process**: crystal plasticity with thermal eigenstrain and
   temperature-dependent slip resistance on the grain structure of link 2, loaded by the
   temperature history of link 1; molten material carries no stress
   (grilli 2022 elimination; hu 2022 residual).
4. **Properties**: a mechanical test of the as-built volume (liu 2020 integration;
   rolchigo 2025 credibility).

Related methods that this plan does not adopt but uses as reference points: kinetic Monte Carlo
grain structure (rodgers 2017 montecarlo); the comparison of phase-field and cellular-automaton
grain competition in a melt pool (elahi 2023 competition); a graph-network form of the
phase-field grain model (xue 2022 graph); GPU finite elements for the thermal problem
(xue 2023 jaxfem; dugast 2021 partscale); the MOOSE Navier-Stokes module as an open
finite-element flow solver (lindsay 2023 navierstokes).

Orientation rule for this part (author, 2026-11-25): grain orientations are drawn at random
where a grain is created (substrate, powder, nucleus), are inherited by epitaxial growth, and
do not change after solidification. All three have published precedent: random orientations for
substrate, powder and nuclei in (yang 2021 grain) and as a stated source of scatter in
(rolchigo 2022 uncertainty); one fixed orientation per grain handed to the mechanics stage in
(grilli 2022 elimination), (hu 2022 residual) and (liu 2020 integration); orientation frozen
below the melting point in (secer 2024 polycrystalline). Measured grains do carry orientation
gradients (polonsky 2020 gradients); the chain does not model them, and the report will say so.

## 2. Benchmark

Revised 2026-11-25 after the author pointed to the local NIST AM Bench data
(`/mnt/extended/projects/ambench-data`, 79 GB). The first version of this plan proposed the
316L chain of (hu 2022 residual); it is now the secondary case (section 2.5).

### 2.1 Proposal: NIST AM Bench 2022 alloy 718, process to microstructure to tensile property

The AM Bench 2022 builds of alloy 718 (AMB2022-01, bridge parts with legs on the NIST
metrology testbed), the bare-plate tracks of the same alloy (AMB2022-03) and the tensile tests
of specimens cut from the 2.5 mm legs of one of those builds (AMB2025-02) form one measured
chain: commanded scan path, melt pool, surface thermal history, grain structure in three
dimensions, stress-strain curve (levine 2024 outcomes). Every link is an experiment (tier 4),
the raw data are public, and most of them are already local. AMB2025-02 was posed by NIST as
exactly this problem: predict the tensile properties from the AMB2022-01 process and
microstructure data, with no mechanical property supplied.

| link | measured anchor | record and paper | local | values checked on 2026-11-25 |
|---|---|---|---|---|
| process input | commanded position, power and trigger at 100 kHz for every layer; plate thermocouples; powder certificate | (lane 2022 builddata) | yes | not yet read |
| melt pool | depth and width of seven single-track cases (power, speed and spot varied), six sections per case | (weaver 2022 trackdata), (weaver 2024 meltpool) | yes | baseline case, 285 W, 960 mm/s, 67 um spot: depth 139.7 um, width 136.3 um (mean of the six sections of the NIST workbook, `derived`) |
| thermal history | time above melting and solid cooling rate per pixel and layer, builds B6 to B8 | (lane 2023 thermographydata) | maps yes, videos no (243 GB) | not yet read |
| grain structure | 3-D EBSD of a 2.5 mm leg in build coordinates, 938 x 969 x 529 voxels of 1 x 1 x 2.12 um, reaching into the plate; large-area 2-D maps; the line-intercept grain-size challenge solution | (schwalbach 2024 sectioningdata), (levine 2022 microstructuredata), (levine 2024 microstructure) | yes | challenge sheet opened (region 1, lines along X: 15,037 intercepts below 100 um, 45 between 100 and 200 um) |
| track-scale grains | EBSD of the single tracks and pads at 0.25 um | (levine 2023 trackebsddata) | yes | not yet read |
| tensile property | eight specimens, 1e-3 /s: E 123.6 +/- 5.9 GPa, 0.2 % yield 690.7 +/- 3.1 MPa, ultimate 951.6 +/- 9.2 MPa, strain at ultimate 0.247 +/- 0.020 (engineering) | (benzing 2025 tensiledata) | yes, raw curves and answer sheet | answer sheet and readme read (`reported_table`) |
| part-scale residual strain | energy-dispersive diffraction on the bridge mid-plane, contour method, deflection after cutting | (levine 2022 residualdata), (phan 2025 residual) | yes | not a target, see 2.3 |
| porosity | X-ray CT of the same leg, 2 um | (schwalbach 2024 sectioningdata) | yes | not a target, see 2.3 |
| thermophysical properties | conductivity, heat capacity, expansion, melting enthalpy of the solid; enthalpy and density of solid and liquid | (agazhanov 2019 thermophysical), (hosaeus 2001 liquid) | PDFs yes | to be extracted by hand before a deck |

Recorded by an earlier analysis of the 3-D EBSD volume (`ambench-data/agent/`, not re-derived
here): about 38,000 grains above the plate with a median equivalent diameter of 11.5 um, and a
texture close to random (fractions of the build direction within 15 degrees of <100>, <110>,
<111>: 0.10, 0.21, 0.10 against 0.10, 0.20, 0.14 for random orientations). Build parameters as
recorded in the author's toolkit from the NIST challenge description, to be read again from the
NIST page and the scan file: 285 W, 960 mm/s, hatch 110 um, layer 40 um, beam diameter 77 um
(the single tracks used 67 um), 90 degree rotation between layers, 312 layers.

### 2.2 Why it is suitable

- **Every stage has a tier 4 anchor with open raw data**, including the last one. The 316L
  chain validates the mechanics stage against lattice strains that depend on an unpublished
  EBSD texture; here the measured grain structure itself is public.
- **Each stage can be tested alone.** The mechanics stage can run on the measured 3-D grain
  structure and on the predicted one, which separates the error of the constitutive stage from
  the error handed down by the thermal and grain stages.
- **Seven melt-pool cases** over power, speed and spot size test the thermal model away from
  the point it was fitted on.
- **The author's toolkit already ran this benchmark** with its own thermal, grain and
  crystal-plasticity codes on the same GPU, so a code-to-code comparison (tier 2) exists for
  the conduction-driven chain, with measured costs.
- **Orientations need no unpublished input.** The plate grains can be drawn from the measured
  plate texture or at random; nuclei are random (section 1).

### 2.3 What the data do not give, and limits to state

- **No grain-scale stress during the process.** The residual strains of
  (levine 2022 residualdata) belong to the whole 75 mm bridge and its plate; a column of
  0.65 mm cannot predict them. The stress left by the process is computed and reported, and its
  check is the secondary case of section 2.5.
- **The baseline track is as deep as it is wide** (139.7 against 136.3 um). That is a
  keyhole-like pool. The toolkit's conduction model, fitted on two cases, gave pools 11 to
  30 % too narrow and, for the three deepest cases, 22 to 41 % too shallow. A flat-surface flow
  model corrects the width (Marangoni) but cannot dig a keyhole; depth will rest on a
  calibrated volumetric source until a deforming surface exists. This is the largest expected
  disagreement of the thermal stage and it is known in advance.
- **Grain size needs a fine grid.** A median grain of 11.5 um is two cells at 5 um; the
  toolkit obtained 16.8 um at that resolution. The coarse check can show the mechanisms
  (remelting, epitaxy, nucleation, a near-random texture); the grain-size comparison is
  meaningful at about 2 um or below, which is the standalone deck.
- **The bridge is not modeled.** The model domain is the column the 3-D EBSD covers
  (650 x 650 um, 0 to 682 um above the plate, 17 layers), with the commanded scan path of that
  location.
- **Elastic modulus.** The reported 123.6 GPa is far below the polycrystal modulus that
  single-crystal constants of alloy 718 give (the toolkit's model gave 201 GPa). Yield,
  ultimate strength and strain at ultimate are the targets; the modulus is reported and not
  fitted.
- **Crystal-plasticity parameters of as-built alloy 718 are not part of the data.** They come
  from the literature (not yet searched) or from a fit; a fit on the tensile curve turns the
  yield stress from a prediction into a calibration, and the report keeps the two apart.
- **Porosity** is measured but no stage of this plan forms pores (no powder, no free surface).

### 2.4 Cases

| case | content | domain | coarse check | standalone deck |
|---|---|---|---|---|
| K1 | single tracks, seven cases, bare plate | 700 x 300 x 345 um | 10 um (72,450 cells), baseline case first | 5 um (579,600 cells) or finer, all seven |
| K2 | build column: thermal history and grain structure | 600 x 600 x 800 um, 15 layers on a 200 um plate slab | 10 um (288,000 cells), 4 layers | 650 x 650 x 880 um, 17 layers, 5 um (2.97 million cells) for the thermal stage, 2 um for the grains |
| K3 | tensile test of a cube cut from the K2 grain structure, and of the same cube cut from the measured 3-D EBSD | 100 to 200 um cube | 10 um | 2 to 5 um, several cubes |
| K4 | stress left by the process in the K2 column | as K2 | 10 um, 4 layers | as K2 |

The domains of K1 and K2 are the ones the toolkit used, so that its results can serve as the
code-to-code reference. Its reduced box of 240 um proved too small for the physics (the pool
fills it and heat piles up at insulated sides); K2 therefore keeps the 600 um width and
reduces resolution and layer count instead.

### 2.5 Secondary case: process stress against a published simulation

The single track of (grilli 2022 elimination), from the 316L chain of (hu 2022 residual)
(thermal-fluid flow, phase-field grains (yang 2021 grain), crystal plasticity in MOOSE, compared
with the lattice strains of (chen 2019 microscale)). It tests the one stage the AM Bench data
cannot: thermal eigenstrain, temperature-dependent slip resistance and the treatment of molten
material in a crystal-plasticity solve. The paper prints what a run needs (`reported_text` and
`reported_table`, read from the arXiv version; to be extracted again by hand before a deck):

| quantity | value |
|---|---|
| laser power, scan speed | 150 W, 700 mm/s, one track |
| mechanics mesh | first-order hexahedra, 4 um (54,000 elements) and 8 um (7,000 elements) |
| cost reported by the authors (MOOSE) | 99 h on 16 CPUs at 4 um; 9.5 h on 12 CPUs at 8 um |
| temperature hand-off | thermal-fluid field sampled every 45 us |
| time step | 1 us, 0.1 us around each temperature update during the scan |
| grains | 4,522 grains from the phase-field model, voxels mapped to elements |
| slip law | power law, reference rate 1e-7 /s, exponent m = 0.1, hardening 3,839 MPa |
| slip resistance | 213 MPa at the reference temperature, saturation 302 MPa, temperature dependent |
| elastic constants (C11, C12, C44) | 204.6, 137.7, 126.2 GPa; slopes -90.33, -45.10, -51.78 MPa/K |
| melting | stiffness degradation, or element elimination and reactivation |

Agreement with it is tier 3. Its thermal-fluid field and its EBSD orientations are not public,
so the comparison is of stress patterns and magnitudes under a substitute temperature field
and random orientations, not of numbers to a tolerance. The stage errors that chain reports
(melt pool within 9 %, grain measures about 10 %, lattice strain within 0.5e-3, as restated in
(hu 2025 process)) are the published yardstick for what such a chain achieves. Still to obtain
for this case: the full texts of (hu 2022 residual) and (yang 2021 grain), the supplement of
(chen 2019 microscale), and the beam size and absorptivity of the single track.

### 2.6 Other candidates considered

| candidate | chain | why not |
|---|---|---|
| (rolchigo 2025 credibility), (rolchigo 2022 uncertainty) | OpenFOAM heat transfer, ExaCA, ExaConstit on the AM Bench 2018 IN625 artifact | cellular automaton grains; the published ensemble is exascale; the alloy 718 data are richer and local |
| (briffod 2025 mapping) | conduction, cellular automaton, crystal-plasticity FE, Hastelloy-X | data not open; abstract only read |
| (liu 2020 integration) | FE thermal, phase field, EVP-FFT, Ti-6Al-4V | needs the alpha/beta transformation stage |
| (sun 2026 coupled) | two-way phase-field and crystal-plasticity coupling | beyond a first capability; abstract only read |
| AM Bench 2018 IN625 tracks (lane 2020 meltpool) | melt pool and cooling rate only | superseded by the alloy 718 tracks, which belong to the same chain |

## 3. Modules required and how much PINNACLE has

Survey of PINNACLE `ddf1e9e`. "Share present" is this plan's estimate of how much of the module
exists as registered, tested components; it is a judgment, not a measurement. Component names
are PINNACLE registry names.

### 3.1 Thermal

| id | module | share present | what exists | what is missing |
|---|---|---|---|---|
| T1 | transient conduction, FE | 50 % | `time_derivative` + `diffusion` + `body_source` on `fem_transient` (implicit Euler, fixed step) | the solve is linear and dense (n x n arrays, `numpy.linalg.solve`): no 3-D process mesh fits |
| T2 | temperature-dependent conductivity | 40 % | `thermal_conductivity` (constant, linear, table, parsed) under the steady Newton solver | refused by `fem_transient`; no nonlinear transient path |
| T3 | temperature-dependent heat capacity and density | 0 % | `time_derivative` takes one constant | material hook for rho c_p(T) |
| T4 | latent heat (enthalpy or apparent capacity, liquid fraction) | 0 % | none | enthalpy form with a smooth liquid fraction, consistent tangent |
| T5 | moving surface Gaussian flux | 70 % | `gaussian_flux` boundary condition, constant velocity | scan-path input; its example deck is translated but not run by the test suite |
| T6 | volumetric source | 40 % | `ellipsoid_path` (one symmetric ellipsoid, constant velocity) | Gaussian beam with exponential depth decay (the toolkit's source), double ellipsoid |
| T7 | scan path (tracks, layers, hatch rotation, laser on/off, from a file) | 5 % | a `parsed` function can encode one raster by hand | a scan-path function component reading way-points `time, x, y, z, power, radius, on` |
| T8 | convection boundary | 100 % | `convective` | - |
| T9 | radiation boundary | 50 % | `radiative` under the steady Newton solver | transient use |
| T10 | evaporative cooling (Hertz-Knudsen) | 0 % | none | boundary condition |
| T11 | powder, liquid, solid state with an irreversible melted flag | 0 % | none | state variable that switches conductivity |
| T12 | layer and track activation | 0 % | none | deposited flag from the scan path (quiet-element form) |
| T13 | adaptive time step | 30 % | cutback on Newton failure in `pf_fem_transient` | step ladder tied to laser on/off; growth during cooling |
| T14 | graded 3-D hexahedral mesh | 40 % | `box_mesh` (uniform), Exodus and Gmsh import, hexaform | grading on `box_mesh` |
| T15 | melt-pool and solidification measures per step | 5 % | extreme values, a volume fraction above a threshold | pool width, depth, length from an isotherm; cooling rate; gradient G and front speed R at freezing |
| T16 | verification against the moving point source solution | 0 % | listed as an acceptance criterion of PINNACLE phase 42, no test | example and test (rosenthal 1946 moving) |

### 3.2 Melt-pool flow

| id | module | share present | what exists | what is missing |
|---|---|---|---|---|
| F1 | incompressible Navier-Stokes | 0 % | none | velocity and pressure solve |
| F2 | advection of heat by the flow | 30 % | `advection` (prescribed constant velocity, no stabilization) | velocity from F1, stabilization |
| F3 | buoyancy (Boussinesq) | 0 % | none | body force |
| F4 | mushy-zone drag (voller 1987 mushy) | 0 % | none | Carman-Kozeny term on the liquid fraction of T4 |
| F5 | Marangoni shear on the pool surface | 0 % | none | boundary traction from the surface temperature gradient |
| F6 | recoil pressure | 0 % | none | boundary term from the vapor pressure of T10 |
| F7 | deforming free surface (level set or volume of fluid) | 0 % | none | not planned in the first pass (section 5, decision D2) |
| F8 | verification (lid-driven cavity, melting with natural convection) | 0 % | none | two example decks with tests |

### 3.3 Grain structure

| id | module | share present | what exists | what is missing |
|---|---|---|---|---|
| M1 | Allen-Cahn and multi-order-parameter grain growth | 90 % | `grain_growth_energy`, `gb_evolution`, `pf_fem_transient` (sparse Newton), `pf_fem_explicit` (explicit, GPU path) | - |
| M2 | temperature field as an input of the phase-field kernels | 10 % | a driving force that depends on time only; a fixed linear temperature gradient | every kernel below needs T(x, t) from the thermal stage |
| M3 | melting: liquid erases the grain parameters | 0 % | none | liquid penalty on the liquid fraction |
| M4 | solidification drive at the pool boundary (epitaxial regrowth) | 0 % | none | kernel; form to follow (yang 2021 grain) after extraction |
| M5 | grain-boundary mobility as a function of the local temperature | 30 % | Arrhenius mobility at one uniform temperature | field-valued temperature |
| M6 | nucleation in the undercooled liquid | 60 % | `nucleation_inserter`, `nucleation_penalty`, `langevin_noise` | rate law in the local undercooling |
| M7 | growth selection along the temperature gradient for cubic crystals | 40 % | `cubic_kinetic_anisotropy`, `ac_interface_anisotropic` | selection factor of the crystal axes against grad T |
| M8 | many grains on few order parameters | 30 % | `grain_tracker` reports remaps and does not apply them | applied remap, or a fixed orientation set (36 in (yang 2021 grain)) |
| M9 | orientation per grain: random, Bunge convention | 90 % | `per_grain_orientation` (uniform), `[OrientationField]` from a point or grain file (ADR 0026, 0036) | hand-off file written by the grain stage |
| M10 | substrate and powder initial grains | 60 % | Voronoi and grown initial conditions | powder layer as a field |
| M11 | grain measures: size, aspect ratio, texture along the build direction | 30 % | feature counts and volumes | aspect ratio, pole-density measures |

### 3.4 Mechanics

| id | module | share present | what exists | what is missing |
|---|---|---|---|---|
| S1 | crystal plasticity, finite strain, implicit update, orientation field | 90 % | `stress_divergence_cp` with `kalidindi` and `dislocation_density_cp`, total-Lagrangian balance, consistent tangent (ADR 0033 to 0036) | - |
| S2 | thermal eigenstrain in the multiplicative split | 0 % | none | F = F^e F^p F^th |
| S3 | temperature-dependent slip resistance and elastic constants | 0 % | none | the forms of (grilli 2022 elimination) Table 2 as material options |
| S4 | temperature field as a load of the mechanics solve | 0 % | none | one-way thermal to mechanics executioner mode |
| S5 | molten material carries no stress; state reset on solidification | 0 % | none | stiffness degradation first, element elimination and reactivation second |
| S6 | J2 plasticity with temperature-dependent yield (cheap check of S2, S4, S5) | 60 % | `j2_isotropic_hardening`, small and finite strain | temperature input, eigenstrain |
| S7 | restart of a mechanics run from a stored state (stress after the process, then a tensile test) | 50 % | checkpoint and restart on the coupled executioners | initial state from a file for a new deck |
| S9 | tensile test of a cube of the grain structure, predicted or measured (case K3) | 80 % | `cp_quasi_static`, `[OrientationField]` from a point file, the route of Parts I to III | a voxel grain map as input (conversion of the 3-D EBSD volume is a benchmark script) |
| S8 | stress and strain measures: lattice strain per grain family, von Mises, residual stress maps | 20 % | stress extremes, `cp_field` output | grain-family average of the elastic strain along a direction |

### 3.5 Coupling, input and output, hardware

| id | module | share present | what exists | what is missing |
|---|---|---|---|---|
| C1 | field hand-off between solves and meshes | 5 % | initial condition from a file at identical node coordinates | a field function that reads a stored time series and interpolates in space and time |
| C2 | field time series output | 60 % | VTK series from the phase-field and quasi-static executioners | series from the thermal executioner; a compact array format for large grids |
| C3 | chain driver | 30 % | staggered and monolithic executioners for mechanics with damage | sequential thermal, grain, mechanics stages from one deck or one script that reads decks |
| C4 | sparse linear solvers | 50 % | SciPy sparse LU, CG with Jacobi, GMRES with ILU | sparse assembly of the thermal operators; a preconditioner that scales; device solver |
| C5 | single-GPU execution | 50 % | Torch backend, batched crystal plasticity, explicit phase-field path on the device | device path for the thermal and flow stages |
| C6 | multi-GPU execution | 0 % | `num_gpus` is validated and then unused | domain split of the explicit structured paths; cannot be tested on `rocm-server` (one GPU) |
| C7 | network routes for the thermal problem | 30 % | `pf_pinn`, `vpinn_transient`, `parametric_surrogate` (data-driven for transient decks) | a moving-source thermal network in the laser frame (zhu 2021 meltpool; liao 2023 hybrid; hosseini 2023 singletrack) |

### 3.6 Summary

| link | share present (estimate) | largest gap |
|---|---|---|
| thermal | about 30 % | a nonlinear, sparse, transient solve with latent heat |
| melt-pool flow | below 5 % | all of it |
| grain structure | about 50 % | coupling to the temperature field |
| mechanics | about 45 % | temperature in the constitutive law and in the split |
| coupling and hardware | about 30 % | field hand-off between stages; multi-GPU |

Three defects found by the survey, to be fixed with the first stage: the moving-source example
decks are not run by any test; a `diffusion` kernel bound to a field-dependent conductivity
appears to pass through `pf_fem_transient` with the dependence dropped (read from the code, not
yet run); two documentation pages say that `num_gpus` above 1 logs a warning, and no code does.

## 4. Development stages in PINNACLE

Each stage is one or more PINNACLE phases under the framework's rules: registered components
with documented `params()`, weak form and tangent, unit tests with an analytical or limiting
case, an example deck that the test suite runs, a translation golden, regenerated syntax pages
and schema, a physics note, a CHANGELOG entry and an ADR (next free number 0047). Nothing
specific to the benchmark enters the framework. Each stage ends with a verification case whose
number is recorded in this repository.

| stage | modules | content | verification case (accuracy tier) |
|---|---|---|---|
| A | T1-T4, T8-T10, T13, C4 | nonlinear transient thermal path: sparse assembly, enthalpy with latent heat, k(T) and rho c_p(T), radiation and evaporation, step control | 1-D solidification front against the Stefan solution (2); nonlinear slab against the steady Newton solver (1) |
| B | T5-T7, T11, T12, T14-T16, C2 | scan path, sources, material state and activation, graded box, melt-pool measures, series output | moving point source (rosenthal 1946 moving) (2); case K1, baseline track against (weaver 2024 meltpool) (4) |
| C | C1, M2-M11 | field hand-off; grain growth driven by the temperature field | static pool: epitaxial regrowth and selection against the gradient (2); case K2 at four layers against the toolkit's result (2) and the 3-D EBSD measures (4, limited by resolution) |
| D | S2-S8, C3 | thermal eigenstrain, temperature-dependent crystal plasticity and J2, molten-state treatment, thermal load mode, state hand-off | constrained bar, stress = -E alpha dT (2); case K3 against (benzing 2025 tensiledata) (4); case K4; single track of (grilli 2022 elimination) at 8 um (3) |
| E | F1-F6, F8 | flat-surface melt-pool flow coupled to stage A | lid-driven cavity (3); melting with natural convection (3, 4); width and depth of the seven K1 tracks (4) |
| F | C7 | network routes on the thermal stage (PINN, PINN with FE correction, parametric surrogate over power, speed and spot of the K1 cases), then on the chain where they apply | against the FE result of the same deck (1), never as validation |
| G | C5, C6 | device paths for stages A and E; multi-GPU split | identical results on one and several devices (1); scaling table |

Order: A, B, C, D give a complete conduction-driven chain and the first end-to-end number
(cases K1 to K4 at the coarse size). E replaces the effective liquid conductivity by a computed
flow and is judged on the seven K1 tracks. F and G follow. The flow stage is
placed after the first chain result on purpose: the chain can be checked while the largest new
solver is built.

Design choices proposed (each becomes an ADR; the author may overrule):

- **Thermal kernels on the existing block Newton executioner.** Temperature becomes a variable
  of `pf_fem_transient` with point terms and exact tangents, so latent heat, k(T), radiation
  and later the flow coupling use one solver, and a temperature variable can be co-solved with
  phase-field variables. `fem_transient` stays as the linear path.
- **An explicit lumped-mass device path for the large structured cases.** At a 5 um mesh the
  explicit conduction limit of alloy 718 is of order 1 us (h^2 / 6a with a about 4.5e-6 m^2/s),
  against a laser step of 5 um / 0.96 m/s = 5.2 us: a handful of sub-steps per laser step, so
  the explicit path of ADR 0039 and 0040 is the cheap route on the GPU for K1 and K2.
- **Flow by a projection method** on the same hexahedral mesh, flat surface, Marangoni shear as
  a boundary traction, drag from the liquid fraction; the toolkit's JAX flow solver is the
  reference for scope (it is a staggered-grid finite-volume code and is rewritten, not ported).
- **Molten material:** stiffness degradation first (no change of topology), element elimination
  and reactivation as the second option, both from (grilli 2022 elimination).
- **Grain model form:** the toolkit's form (multi-order-parameter grain-growth energy with a
  liquid penalty, a solidification drive, a selection factor against the temperature gradient
  and a nucleation rate in the undercooling), because it has been run on this benchmark and
  gives the code-to-code reference; the form of (yang 2021 grain) is the alternative once its
  full text is read.

## 5. Decisions asked of the author

| # | question | proposal |
|---|---|---|
| D1 | benchmark | NIST AM Bench alloy 718: AMB2022-03 tracks, AMB2022-01 build column, AMB2025-02 tensile tests (cases K1 to K4); the 316L single track of (grilli 2022 elimination) as the secondary case for the process stress |
| D2 | scope of the flow stage | flat surface with Marangoni, recoil, buoyancy and drag; no deforming free surface, no powder particles. The baseline track is keyhole-like, so depth stays tied to a calibrated source (section 2.3); a deforming surface would be a later part |
| D3 | stage order | A, B, C, D first (conduction-driven chain on K1 to K4), then E (flow) |
| D4 | orientations | plate grains: uniform random with recorded seeds, or drawn from the measured plate slab of the 3-D EBSD (public data); nuclei: uniform random; several realizations |
| D5 | heat-source calibration | absorptivity and penetration depth fitted on two of the seven tracks, the other five held out (the toolkit's design); or all seven reported without a fit first |
| D6 | multi-GPU | designed and unit-tested here on one device and on CPU partitions; scaling measured by the author on other hardware |
| D7 | crystal-plasticity parameters of as-built alloy 718 | literature set first and the tensile result reported as a prediction, then a fit on the measured curve reported as a calibration; elastic constants from the literature, modulus not fitted |
| D8 | secondary 316L case | keep it as the tier 3 check of stage D, or drop it and report the process stress of K4 without an external comparison |
| D9 | part number and names | this is Part V; report `docs/reports/lpbf_process_pinnacle_benchmark_report.md`, decks under `decks/lpbf/`; NIST data stay in `ambench-data` and are read by path, not copied (their licence statement is cited in the report) |

## 6. Benchmark work in this repository

- `decks/lpbf/` holds one deck per verification case and per chain stage, and
  `decks/lpbf/standalone/` the full-accuracy decks. Scripts read decks (`scripts/lpbf_*.py`).
- Order: one checked case per stage before any series (CLAUDE.md rule 5); every run is a row
  of `docs/run_ledger.md`.
- Results in this document's successor notebook and in the report; reproduction steps in
  `docs/reproduce_lpbf.md`; workflow in `docs/workflow_lpbf.md`.
- Every chain result ships images of the temperature, grain and stress fields and an animation
  of the scan; network results over several seeds.
- Hardware: `rocm-server` (16 cores, 32 threads, 62 GB, one Radeon 8060S, Torch 2.13 ROCm 7.2),
  idle on 2026-11-25, PINNACLE there on `feat/pinn_benchmark` at `d41c76f`; the new branch is
  synced before the first run.

Sizes are in section 2.4. No PINNACLE timing exists for these stages. Measured indications from
the toolkit's own codes on this benchmark (its manual and benchmark document, not remeasured
here): 0.017 s per thermal step and 0.65 s per grain-growth step at 129,024 cells on the same
GPU in single precision; the 2.30 million cell, 15 layer column without crystal plasticity in
about four hours on a desktop GPU; grain-growth cost growing as the cell count to the power
1.7; crystal plasticity four to five times the grain-growth cost per layer. The full column at
5 um was estimated there at about 13 h per layer on the Radeon 8060S and was not run. Cost is
measured on the first case of each stage before a series is proposed.

## 7. Risks

- The limits of section 2.3: keyhole-like baseline track, grain size below the coarse grid,
  no measured grain-scale process stress, crystal-plasticity parameters not in the data.
- Values of the secondary case and of the other candidates come partly from abstracts and from
  one extraction tool; the alloy 718 values not marked as checked in section 2.1 come from the
  toolkit's documents and an earlier analysis. No value enters a production deck before it is
  read in the primary source or the NIST file.
- The thermography videos (243 GB) are not local; the time-above-melting and cooling-rate maps
  are. A pool length from the videos would need a download.
- The implicit crystal-plasticity update has not been run with a temperature that changes by
  more than 1,000 K within a few steps; step control near melting is the expected difficulty,
  as the published run needed a tenfold smaller step at each temperature update.

## 8. References used

Papers: (agazhanov 2019 thermophysical), (briffod 2025 mapping), (chadwick 2021 grain),
(chen 2019 microscale), (dugast 2021 partscale), (elahi 2023 competition),
(gandin 1994 cellular), (grilli 2022 elimination), (hosaeus 2001 liquid), (hu 2022 residual),
(hu 2025 process), (khairallah 2016 meltflow), (lane 2020 meltpool), (levine 2024 microstructure),
(levine 2024 outcomes), (lindsay 2023 navierstokes), (liu 2018 sitespecific),
(liu 2020 integration), (phan 2025 residual), (polonsky 2020 gradients),
(rodgers 2017 montecarlo), (rolchigo 2022 exaca), (rolchigo 2022 uncertainty),
(rolchigo 2025 credibility), (rosenthal 1946 moving), (secer 2024 polycrystalline),
(sun 2026 coupled), (voller 1987 mushy), (weaver 2024 meltpool), (xue 2022 graph),
(xue 2023 jaxfem), (yang 2021 grain), (hosseini 2023 singletrack), (liao 2023 hybrid), (zhu 2021 meltpool).

NIST data records: (benzing 2025 tensiledata), (lane 2022 builddata),
(lane 2023 thermographydata), (levine 2022 microstructuredata), (levine 2022 residualdata),
(levine 2023 trackebsddata), (schwalbach 2024 sectioningdata), (weaver 2022 trackdata).

Keys and reading status of each (full text, abstract, metadata only, data record) are in
`bib/citation_keys.md`.
