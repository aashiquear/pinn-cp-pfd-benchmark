# Reproduce: Al2219 polycrystal plate, Parts I and II

Commands that rebuild every result of
`docs/reports/cp_pfd_pinnacle_benchmark_report.md` from the shipped decks. Run
from the study root with PINNACLE `feat/pinn_benchmark` at 9548123 or later
(Bunge Euler angles, ADR 0036; total-Lagrangian balance, ADR 0035; implicit
update, ADR 0033) and hexaform `feat/plate-hole`. Units are micrometer, MPa,
second. Solver output is written under `results/cp_study/bunge/` and
`results/pfd_study/` and is not shipped; the figures under `figures/report/`
and `figures/pfd/` are. Each step has a row `al-bunge-*` in
`docs/run_ledger.md` with its host, cores and run time. The J2 runs of Part II
(report §13) are those of the first study (rows `pfd-j2-*`); their decks are
under `decks/pfd_study/j2/`.

The compute host has PINNACLE and torch but no matplotlib and no hexaform:
meshes, surrogate designs, query files and figures are made on the
workstation and the results are copied between the two with rsync.

## Environment

    export PATH=$HOME/venvs/pinnacle/bin:$PATH
    export PYTHONPATH=$HOME/projects/pinnacle_bunge/src     # PINNACLE 9548123
    export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2              # every single `pinnacle run`

A torch run without a thread cap opens one thread per core; several such runs
side by side are slower than one by a factor of fifty. The study drivers set
the cap from the `threads` key of their deck. `pinnacle run` resolves relative
paths against the working directory, so a generated deck is run from its own
folder.

## Part I: crystal plasticity on the plain plate

### 1. Grain structures (report §2.6, §3.6)

    pinnacle run decks/graingrowth/pilot/n32_bdf1_h40.i
    pinnacle run decks/graingrowth/instances/n32_bdf1_h40_sB.i
    pinnacle run decks/graingrowth/instances/n32_bdf1_h40_sC.i

### 2. Calibration (§4)

    python scripts/cp_parameter_study.py decks/cp_study/bunge/cal_gn1.i all
    python scripts/cp_gauss_newton.py    decks/cp_study/bunge/cal_gn1.i     # writes cal_center2.i, cal_gn2.i
    python scripts/cp_parameter_study.py decks/cp_study/bunge/cal_gn2.i all
    python scripts/cp_gauss_newton.py    decks/cp_study/bunge/cal_gn2.i     # predicted gain 0.1 MPa: converged

The calibrated set is `decks/cp_study/bunge/cal_center2.i`; its run is
`results/cp_study/bunge/cal/gn2/base`.

### 3. Mesh sensitivity (§5) and grain instances (§7)

    for d in mesh_h20 mesh_h80 mesh_h160 grain_sB grain_sC; do
      python scripts/polycrystal_plate.py decks/cp_study/bunge/$d.i label orient deck
      out=$(grep out_dir decks/cp_study/bunge/$d.i | awk '{print $3}')
      (cd decks/cp_study/bunge/$out && pinnacle run cpfe.i)
      python scripts/polycrystal_plate.py decks/cp_study/bunge/$d.i curve
    done

Use 4 threads for h80 and 8 for h160 (26 min and 1.9 h).

### 4. Parameter sensitivity (§6)

    python scripts/cp_parameter_study.py decks/cp_study/bunge/sens_h40.i all

### 5. Neural routes (§8) and their sensitivity (§9)

    for r in linear fem_residual fem_solver fem_residual_grain fem_solver_grain; do
      python scripts/polycrystal_plate.py decks/cp_study/bunge/hybrid/h40_$r.i label orient deck
      (cd results/cp_study/bunge/hybrid/h40_$r && OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 pinnacle run cpfe.i)
      python scripts/polycrystal_plate.py decks/cp_study/bunge/hybrid/h40_$r.i curve
    done
    python scripts/hybrid_compare.py decks/report/bunge/hybrid_compare.i
    python scripts/field_compare.py  decks/report/bunge/neural_compare.i
    python scripts/cp_parameter_study.py decks/cp_study/bunge/sens_h40_linear.i all
    python scripts/cp_parameter_study.py decks/cp_study/bunge/sens_h40_grain_solver.i all
    python scripts/neural_sensitivity.py decks/report/bunge/neural_sensitivity.i

### 6. Hardware benchmark (§10)

    python scripts/hardware_benchmark.py decks/benchmark/bunge/hardware_h40.i          # alone on the host
    python scripts/hardware_benchmark.py decks/benchmark/bunge/hardware_h40.i --plot   # figure, workstation

### 7. Surrogate (§11)

    python scripts/surrogate_design.py  decks/cp_study/bunge/surrogate/design_h40.i design
    python scripts/cp_parameter_study.py results/cp_study/bunge/surrogate/h40_data/study.i all
    python scripts/surrogate_design.py  decks/cp_study/bunge/surrogate/design_h40.i collect
    python scripts/surrogate_queries.py decks/cp_study/bunge/surrogate/design_h40.i 1000 results/cp_study/bunge/surrogate/queries_1000.csv
    for t in curve curve_cpu field; do pinnacle run decks/cp_study/bunge/surrogate/plate_h40_surrogate_$t.i; done
    python scripts/surrogate_compare.py decks/report/bunge/surrogate_compare.i

### 8. Figures of Part I

    python scripts/report_figures.py          decks/report/bunge/cp_calibration_figures.i
    python scripts/plastic_strain_animation.py decks/report/bunge/plastic_strain_animation.i

## Part II: phase-field ductile fracture on the perforated plate

### 9. Meshes (§12.3, §16)

    for m in h20_hole5 h40_hole5 h80_hole5 h40_ell_x2.5_y5 h40_ell_x1.25_y5 h40_ell_x3.54_y7.07 \
             h40_ell_x2.5_y10 h40_ell_x2.5_y5_a45 h40_ell_x2.5_y5_a90; do
      python scripts/plate_hole_mesh.py decks/pfd_study/mesh/$m.i
    done

### 10. Calibration of the fracture energy (§14.1)

    python scripts/pfd_calibrate.py decks/pfd_study/cp/bunge/cal_h40.i

The value at 5 um is entered in `decks/pfd_study/cp/bunge/base_h40.i`
(`scripts/pfd_fill_base.py` does it from `calibrated.json`).

### 11. Studies (§14 to §17)

    for s in neural_h40 sweep_l_h40 grains_h40 shape_h40 sens_h40 sweep_h_l5; do
      python scripts/pfd_study.py decks/pfd_study/cp/bunge/$s.i all
    done

This is `scripts/tools/al_bunge_pfd_chain.sh` (5.5 h on the host).

### 12. Hardware benchmark (§18) and surrogate (§19)

    python scripts/surrogate_design.py  decks/pfd_study/cp/bunge/surrogate/design_h40.i design      # workstation
    python scripts/surrogate_queries.py decks/pfd_study/cp/bunge/surrogate/design_h40.i 1000 results/pfd_study/cp/bunge/surrogate/queries_1000.csv
    bash scripts/tools/al_bunge_pfd_chain2.sh                                                        # host
    python scripts/hardware_benchmark.py decks/benchmark/bunge/hardware_pfd_h40.i --plot
    python scripts/surrogate_compare.py  decks/pfd_study/cp/bunge/surrogate/compare_h40.i

The chain runs the hardware cases alone on the host, then the 24 samples,
their collection and the two trainings.

### 13. Figures of Part II

    python scripts/pfd_calibration_figures.py decks/pfd_study/cp/bunge/cal_h40_figures.i
    for s in neural_h40 sweep_l_h40 grains_h40 shape_h40 sens_h40 sweep_h_l5; do
      python scripts/pfd_study_figures.py decks/pfd_study/cp/bunge/${s}_figures.i
    done
    python scripts/pfd_fields.py decks/pfd_study/report/bunge/fields_cal_h40.i
    python scripts/pfd_fields.py decks/pfd_study/report/bunge/fields_shape_h40.i

The field decks need the VTK series of `neural_h40/base` and of the seven
shape cases, and the J2 series of `results/pfd_study/j2/sens_h40_wc200/base`.

## Provenance record

| Item | Where |
|---|---|
| Growth seeds, grid, material, scheme, step, horizon | the growth deck |
| Labeling threshold, orientation seed, drawn triples | `[Labeling]`, `[Orientations]` of the pipeline deck; `grains.csv` beside the run |
| CP parameters, loading, mesh, tolerances, balance, update | `[Mechanics]` of the pipeline deck, copied into `cpfe.i` |
| Fracture parameters, staggered settings | `[Plate]` of `decks/pfd_study/cp/bunge/base_h40.i`, copied into `pfd.i` |
| Framework revision, host, cores, time | `docs/run_ledger.md` |
| First study (exchanged angles, earlier balance) | `archive/decks/al2219_exchanged_angles/`, `archive/docs/al2219_report_exchanged_angles_2026-11-16.md` |
