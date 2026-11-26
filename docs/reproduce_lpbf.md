# Reproducing Part V (laser powder bed fusion process chain)

Requirements: PINNACLE on branch `feat/multi-gpu-and-am-process` (0d7befe or later)
importable (`pip install -e`, or `PYTHONPATH=<pinnacle>/src`), Python 3.10+, NumPy, SciPy,
PyTorch for the GPU and mechanics runs, matplotlib for figures. Commands run from the
root of this repository. On an AMD GPU export `TORCH_BLAS_PREFER_HIPBLASLT=0` and
`TORCH_ROCM_AOTRITON_ENABLE_EXPERIMENTAL=1`.

## 1. Verification (PINNACLE)

```bash
cd <pinnacle>
pytest -q tests/unit/components/test_additive_thermal.py tests/unit/components/test_additive_grains.py \
          tests/unit/components/test_additive_thermomechanics.py tests/unit/framework/test_device_pool.py \
          tests/integration/fem/test_additive_process.py
pinnacle run examples/fem/additive_process/melting_front_stefan_1d.i      # and the other decks there
```

## 2. Single tracks (case K1)

```bash
# targets from the NIST workbook (path in the study deck; record mds2-2718)
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10.i targets
# decks, runs (GPU: 23 s per case at 10 um), table, figure
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10.i decks
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10.i run
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10.i collect
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10.i figure
python3 scripts/lpbf_figures.py track results/lpbf/in718/k1/h10/case0/process_state.npz 139.7 136.3 figures/lpbf/k1_case0_section_h10.png
# the starting source and the fit step that led to the calibrated deck
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10_start.i decks
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10_start.i run
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_study_h10_start.i fit
```

A single generated deck runs alone: `pinnacle run decks/lpbf/in718/k1/cases/h10/case0.i`.
For a CPU run set `backend = numpy` and `use_gpu = false` in the study deck (141 s per
case at 10 um on eight threads). Standalone decks at 5 um: `k1_study_h5.i`.

## 3. Build column (case K2)

```bash
python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10.i run       # four layers, 10 um
python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10.i grains
python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10.i figure
# one-layer checks with another selection exponent
python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10.i run utilities.scan.layers=1 \
        executioner.end_time=5.6e-3 kernels.grains.selection_exponent=1.0 \
        global.output_dir=../../../../results/lpbf/in718/k2/layer1_sel1.0
```

Standalone decks of the full measured column: `k2_column_full_h5.i` (2.97 million
elements, float64, about 80 GB, 162,000 steps) and `k2_column_full_h2p5.i` (23.8 million
elements, float32). Add `devices = 'cuda:0 cuda:1 ...'` to `[Global]` to split the element
loop over several GPUs. They have not been run for this report.

## 4. 316L single track

```bash
cd decks/lpbf/ss316l && pinnacle run track_thermal_h8.i && cd ../../..      # 138 s on a CPU
python3 scripts/lpbf_orientations.py decks/lpbf/ss316l/grains_seed1.csv 200e-6 160e-6 108e-6 16e-6 1
python3 scripts/lpbf_mechanics.py decks/lpbf/ss316l/track_mechanics_h8.i run   # scan, cooling, wall release, relaxation
python3 scripts/lpbf_figures.py stress results/lpbf/ss316l/mechanics_h8 \
        results/lpbf/ss316l/thermal_h8/process_state.npz figures/lpbf/ss316l_stress_h8.png
```

The mechanics deck has two stages in one run: 178 steps with the side walls fixed, then
the walls x = 200 um and y = 160 um are released (`release_time` on the two boundary
conditions) and 22 relaxation steps follow. `history.csv` holds the mean stress of every
step; the stress field written is the one after the relaxation.

The 4 um decks (`track_thermal_h4.i`, `track_mechanics_h4.i`) are the paper's fine mesh
and have not been run for this report.

## 5. Single tracks with melt-pool flow (case K1, section 6.5 of the report)

```bash
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_flow_study_h10.i decks     # 7 decks from k1_flow_template.i
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_flow_study_h10.i run       # 510 s each on the GPU host
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_flow_study_h10.i collect
python3 scripts/lpbf_k1.py decks/lpbf/in718/k1/k1_flow_study_h10.i figure
# animation (GIF and PNG frame stack) of one track
python3 scripts/lpbf_animate.py track results/lpbf/in718/k1/flow_h10/case0 \
        figures/lpbf/animations/k1_flow_h10_case0.gif
# one case without the recoil pressure
python3 scripts/lpbf_run.py decks/lpbf/in718/k1/cases/flow_h10/case0.i \
        kernels.marangoni.recoil_coefficient=0.0 \
        global.output_dir=../../../../../../results/lpbf/in718/k1/flow_h10_norecoil/case0
```

Every case is run by `scripts/lpbf_run.py`, which keeps the field series of the deck's
`snapshots_every` in `fields.npz`. The verification decks of the flow kernels are PINNACLE
examples (`examples/fem/additive_process/lid_driven_cavity_2d.i`,
`heated_cavity_convection_2d.i`, `thermocapillary_layer_2d.i`).

## 6. Tensile case and calibration (case K3, section 6.6)

```bash
# targets from the raw NIST records
python3 scripts/lpbf_k3_targets.py "<ambench-data>/nist/DATA_AMB2025_02/raw tensile data" \
        data/lpbf/in718/k3_tensile_mean_curve.csv data/lpbf/in718/k3_tensile_summary.json
# grains of the cube (three seeds)
python3 scripts/lpbf_orientations.py decks/lpbf/in718/k3/grains_seed1.csv 120e-6 120e-6 120e-6 20e-6 1
# starting set, then the Gauss-Newton calibration (4 runs per iteration, in parallel)
python3 scripts/lpbf_k3.py run decks/lpbf/in718/k3/k3_cube_n6_start.i
python3 scripts/lpbf_k3.py calibrate decks/lpbf/in718/k3/k3_cube_n6_start.i \
        data/lpbf/in718/k3_tensile_mean_curve.csv 6
# calibrated set on three grain seeds, and the figure
for s in 1 2 3; do python3 scripts/lpbf_k3.py run decks/lpbf/in718/k3/k3_cube_n6_seed$s.i; done
python3 scripts/lpbf_k3.py figure data/lpbf/in718/k3_tensile_mean_curve.csv \
        figures/lpbf/k3_tensile_calibration.png "start=results/lpbf/in718/k3/start_n6" \
        "seed 1=results/lpbf/in718/k3/calibrated_n6_seed1"
# the chain: a cube of the simulated column
python3 scripts/lpbf_k3.py cube results/lpbf/in718/k2/column_h10 \
        decks/lpbf/in718/k3/grains_k2_column_seed1.csv 240e-6 240e-6 220e-6 120e-6
python3 scripts/lpbf_k3.py run decks/lpbf/in718/k3/k3_cube_k2structure_seed1.i
```

The decks `k3_cube_n6_seed<n>.i` hold the calibrated values; `cal_<i>_<name>.i` are the
decks of the calibration iterations.

## 7. Grain realizations and animations of the column (case K2)

```bash
for s in 2 3; do
  python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10_seed$s.i run
  python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10_seed$s.i grains
  python3 scripts/lpbf_column.py decks/lpbf/in718/k2/k2_column_h10_seed$s.i thermal
done
python3 scripts/lpbf_animate.py grains results/lpbf/in718/k2/column_h10 figures/lpbf/animations/k2_column_h10.gif
```

## 8. Hosts and revisions

Every run with its deck, script, PINNACLE revision, host, precision and wall time is a
row of `docs/run_ledger.md` (rows `lpbf-*`).
