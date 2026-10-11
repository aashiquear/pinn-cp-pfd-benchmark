# Reproduce: phase-field polycrystal plate and its CPFE run

Everything below runs on the ROCm host (`ssh rocm-server`) with the study
repository at `~/projects/pinn-cp-pfd-benchmark` and PINNACLE at
`~/projects/pinnacle` on branch `feat/pinn_benchmark`. Commands are given
from the study root. Local file paths in decks are relative to the deck.

## Environment

```
python3 -m venv --without-pip ~/venvs/pinnacle
curl -sS https://bootstrap.pypa.io/get-pip.py | ~/venvs/pinnacle/bin/python -
cd ~/projects/pinnacle && git checkout feat/pinn_benchmark && ~/venvs/pinnacle/bin/pip install -e ".[dev,meshio]"
export PATH=~/venvs/pinnacle/bin:$PATH
```

PINNACLE revision: at or after `14dd45c` (needs ADR 0026's
`[OrientationField] source` and the `cp_quasi_static` history columns).
Backend numpy, float64, CPU. Set `OMP_NUM_THREADS=8` per run.

## 1. Verification of the growth scheme (optional, 20 min)

```
for c in D_M1_dt8.00_n1 A_M1_dt0.25_n32 B_M1_dt0.50_n16 C_M2_dt0.25_n16 \
         E_bdf1_M1_dt0.25_n32 F_bdf1_M1_dt0.50_n16 G_bdf1_M2_dt0.25_n16; do
  pinnacle run decks/graingrowth/verify/$c.i
done
```
Compare the `energy` column of `results/graingrowth_verify/*/solver_history.csv`.

## 2. Grain growth pilot (three horizons, ~1 h wall, concurrent)

```
for n in 40 80 160; do pinnacle run decks/graingrowth/pilot/n32_bdf1_h$n.i & done; wait
```
Outputs in `results/graingrowth_pilot/n32_bdf1_h{40,80,160}/` (`fields.vtk`,
`solver_history.csv`, `postprocessors.csv`).

## 3. Label, orient, deck

```
python scripts/polycrystal_plate.py decks/polycrystal_plate/plate_n32_h40.i label orient deck
```
Writes `results/polycrystal_plate/n32_h40/{grains.npz, label.json,
grains.csv, orientations_points.csv, orient.json, cpfe.i}`.

## 4. CPFE run and curve

```
cd results/polycrystal_plate/n32_h40 && pinnacle validate cpfe.i && pinnacle run cpfe.i && cd -
python scripts/polycrystal_plate.py decks/polycrystal_plate/plate_n32_h40.i curve
```
`results.csv` holds step, load factor, nominal strain, cutbacks, time and
the mean stress/strain components.

## 5. MOOSE reference curve

```
python scripts/tools/exodus_summary.py moose/CPFE_out.e --pull-axis y \
    --csv data/reference/moose_cpfe_plain_plate_stress_strain.csv
```

## 6. CP parameter study (calibration, mesh, sensitivity, grain instances)

All study decks are under `decks/cp_study/` and ship with the repository.

```
# calibration iterations (center + one perturbation per parameter)
python scripts/cp_parameter_study.py decks/cp_study/cal_gn1.i prepare run collect
python scripts/cp_gauss_newton.py decks/cp_study/cal_gn1.i          # writes cal_center2.i, cal_gn2.i
python scripts/cp_parameter_study.py decks/cp_study/cal_gn2.i prepare run collect

# mesh sensitivity with the calibrated set (h40 = results/cp_study/cal/gn2/base)
for m in h20 h80 h160; do
  python scripts/polycrystal_plate.py decks/cp_study/mesh_$m.i label orient deck
  (cd results/cp_study/mesh/$m && pinnacle run cpfe.i)
  python scripts/polycrystal_plate.py decks/cp_study/mesh_$m.i curve
done

# parameter sensitivity at h40 (base + 16 one-at-a-time cases)
python scripts/cp_parameter_study.py decks/cp_study/sens_h40.i prepare run collect

# grain instances: new seedings grown, then the calibrated set on each
pinnacle run decks/graingrowth/instances/n32_bdf1_h40_sA.i
pinnacle run decks/graingrowth/instances/n32_bdf1_h40_sB.i
for t in A B; do
  python scripts/polycrystal_plate.py decks/cp_study/grain_s$t.i label orient deck
  (cd results/cp_study/grain_instance/s$t && pinnacle run cpfe.i)
  python scripts/polycrystal_plate.py decks/cp_study/grain_s$t.i curve
done
```

## 7. Field series and the plastic-strain animation

```
for m in h20 h40 h80 h160; do
  python scripts/polycrystal_plate.py decks/cp_study/fields_$m.i label orient deck
  (cd results/cp_study/fields/$m && pinnacle run cpfe.i)      # writes fields_step*.vtk + fields.pvd
done
python scripts/plastic_strain_animation.py decks/report/plastic_strain_animation.i
```

## 8. Report figures

```
python scripts/report_figures.py decks/report/cp_calibration_figures.i   # -> figures/report/
```

## Provenance record

| Item | Where |
|---|---|
| Growth seeds, grid, material, scheme, dt, horizon | the growth deck header and body |
| Labelling threshold and speck size | `[Labeling]` in the pipeline deck |
| Orientation seed and the drawn triples | `[Orientations]`, `grains.csv` |
| CP parameters, loading, mesh, tolerances | `[Mechanics]`, copied verbatim into `cpfe.i` |
| PINNACLE revision | `git -C ~/projects/pinnacle rev-parse HEAD` at run time; recorded in the results document |
