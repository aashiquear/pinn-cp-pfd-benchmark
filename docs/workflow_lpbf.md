# Workflow of Part V (laser powder bed fusion process chain)

How the cases of `docs/reports/lpbf_process_pinnacle_benchmark_report.md` are organized.
Commands are in `docs/reproduce_lpbf.md`; results and their reading are in the report.

## 1. Pieces

| piece | where | role |
|---|---|---|
| PINNACLE components | `/mnt/extended/projects/pinnacle`, branch `feat/multi-gpu-and-am-process` | `pinnacle.components.additive` (ADR 0047 to 0050): material, kernels, services |
| NIST records | `/mnt/extended/projects/ambench-data` (not copied here) | measurements; read by path |
| derived targets | `data/lpbf/in718/k1_single_tracks.csv` | written by `scripts/lpbf_k1.py targets` from the NIST workbook |
| decks | `decks/lpbf/in718/k1`, `decks/lpbf/in718/k2`, `decks/lpbf/ss316l` | one PINNACLE deck per run; study decks and a template for the generated ones |
| scripts | `scripts/lpbf_k1.py`, `lpbf_column.py`, `lpbf_mechanics.py`, `lpbf_orientations.py`, `lpbf_figures.py` | generate decks, run, reduce, draw |
| results | `results/lpbf/...` (not shipped) | PINNACLE output and reduced tables |
| figures | `figures/lpbf/` | shipped |

## 2. The chain

1. **Thermal deck.** Material tables, scan path, state service, the four heat kernels.
   Alone for the single tracks (K1) and the 316L track; together with the grain variable
   for the build column (K2). Output: per-element record (`process_state.npz`: state,
   peak temperature, melt count, freezing gradient, cooling rate and direction) and, when
   asked for, the nodal temperature series (`temperature_series.npz`).
2. **Grain stage.** Not a separate run: the order parameters are a second variable of the
   thermal deck. Output: final order parameters and the orientation of each
   (`fields.npz`), reduced to grains by `lpbf_column.py grains`.
3. **Mechanics deck.** Crystal plasticity with `thermal_load`; reads the temperature
   series through a `field_series` function and the grain orientations through an
   `[OrientationField]` point file. Output: stress, plastic strain, history
   (`lpbf_mechanics.py`).

## 3. Study decks and generated decks

`k1_study_<tag>.i` lists the seven NIST cases, the mesh size, the back end and the heat
source; `lpbf_k1.py <study> decks` fills `k1_template.i` and writes
`cases/<tag>/case<line>.i`. A generated deck is complete: `pinnacle run` runs it. The
fit (`lpbf_k1.py <study> fit`) writes its runs under `results/.../fit/` and a
`k1_fit_step.json` with the residuals, the Jacobian and the step; the step is then typed
into the next study deck (the calibrated one names where it came from).

Column and mechanics decks are written by hand. Check runs change a deck key on the
command line (`lpbf_column.py <deck> run utilities.scan.layers=1 ...`); the override is
recorded in `wall.json`.

## 4. Rules followed

- One checked case before a family of runs; every run is a row of `docs/run_ledger.md`.
- A value enters a deck with its provenance (`docs/16_lpbf_notebook.md`, sections 2 and 7).
- A grain statistic is quoted only after its section image has been looked at
  (notebook, section 6: a first run filled empty space with grains and the summary
  numbers did not show it).
- Runs on the GPU host use `backend = torch`, `use_gpu = true`, `precision = float64`,
  with `TORCH_BLAS_PREFER_HIPBLASLT=0` and `TORCH_ROCM_AOTRITON_ENABLE_EXPERIMENTAL=1`.

## Additions of the second stage (2026-11-25)

- **Flow.** A track deck with flow is the conduction deck plus two variables (velocity,
  pressure) and four kernels; `k1_flow_template.i` is filled by the same script from a
  study deck that also gives the artificial sound speed, the speed bound for the thermal
  step and the surface-tension coefficient. Runs go through `scripts/lpbf_run.py`.
- **Animations.** Every run made through `lpbf_run.py` or `lpbf_column.py` keeps a field
  series (`snapshots_every` in the deck). `scripts/lpbf_animate.py` turns it into a GIF and
  a PNG frame stack under `figures/lpbf/animations/`; a report figure of a transient is a
  frame of that stack.
- **Tensile case.** `scripts/lpbf_k3_targets.py` derives the target curve from the raw
  records; `scripts/lpbf_k3.py` runs a cube deck, calibrates (one deck per run of every
  iteration, written next to the start deck), draws the curves, and cuts a cube of grains
  out of a column run.
- **Seeds.** A grain realization is a deck of its own (`*_seed<n>.i`) with its own output
  directory; three seeds per stochastic case.

## Additions of the third stage (2026-11-26)

- **Defect indicators.** A report is a `[Utilities]` block of a thermal deck, listed
  after the process state: `defect_indicators` (lack of fusion, keyhole mode, balling)
  and `cracking_indicators` (one block per solid-fraction curve). `scripts/lpbf_defects.py`
  turns the reports into the tables and figures. A run that is to carry cracking
  indicators has to cool until nothing is liquid.
- **Dwell.** The scan path has a dwell between tracks and between layers. In the decks
  for the multi-GPU machine the next layer is spread at the end of the layer dwell
  (`recoat = dwell_end`), the scan is held at the start of every track until the
  material is solid (`hold_scan_until_solid`), and the run ends when the path is over
  and all is solid. `scripts/lpbf_column.py <deck> dwell` checks it after a run.
- **Several GPUs.** `global.rank_devices="cuda:i cuda:j ..."` cuts the mesh of any deck
  of the chain into one part per card (PINNACLE ADR 0057, library `gpumpi`); the result
  equals the single run. The cards are named at run time (`scripts/gpu_free.py`), not in
  the deck. `scripts/lpbf_parts_check.py` compares a single run with its cut run.
- **Memory.** Every run on a shared or small machine goes through
  `scripts/memory_guard.py`. Fields of two runs are compared only after their shapes
  have been checked (a scalar field is `(n,)`).
- **Generated decks of the wide column**: `scripts/lpbf_column_wide.py <element size>`;
  the header of each deck states elements, memory per card and steps.

