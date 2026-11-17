# Reproduce: Ni single-crystal dog-bone specimens

Commands that rebuild every result of
`docs/reports/ni_single_crystal_pinnacle_benchmark_report.md` from the shipped
decks. Run from the study root with PINNACLE `feat/pinn_benchmark` at 9548123
or later (Euler angles in Bunge's convention, ADR 0036) and hexaform
`feat/plate-hole` (1cb1671 or later). Units are mm, MPa, s. Solver output is
written under `results/ni_study/bunge/` and is not shipped; the figures under
`figures/ni/` are. Each step has a row `ni-bunge-*` in `docs/run_ledger.md`
with its host, cores and run time.

## 1. Meshes

    python3 scripts/ni_specimen_mesh.py decks/ni/mesh/h020.i
    python3 scripts/ni_specimen_mesh.py decks/ni/mesh/h010.i

Writes `results/ni_study/mesh/<level>/<specimen>/mesh.e` for the 15 specimens
and `figures/ni/mesh_<level>.png`.

## 2. Crystal plasticity on L.1 (report §4)

    python3 scripts/ni_specimen.py decks/ni/cp/bunge/L1_h020_rate1e-3.i deck run curve
    python3 scripts/ni_study.py    decks/ni/cp/bunge/cal_gn1.i all
    python3 scripts/ni_specimen.py decks/ni/cp/bunge/cal_center2a.i deck run curve
    python3 scripts/ni_pfd_curve_figures.py decks/ni/figures/bunge/cal_curve.i

The rate decks (`rate1e-4`, `rate3e-4`, `rate1e-3`) carry the paper's set.
`cal_gn1.i` perturbs the three free values around it; `cal_center2a.i` is the
carried set.

## 3. Fracture energy on L.1 (report §6)

    python3 scripts/ni_pfd_calibrate.py decks/ni/pfd/bunge/cal_L1_h020.i
    python3 scripts/ni_pfd_curve_figures.py decks/ni/figures/bunge/pfd_curves_L1_cal.i
    python3 scripts/ni_pfd_shape_figures.py decks/ni/figures/bunge/pfd_shape_L1_cal_damage.i

## 4. Mesh check (report §5)

    python3 scripts/ni_specimen.py decks/ni/pfd/bunge/mesh_L1_h010.i deck run curve fracture
    python3 scripts/ni_pfd_curve_figures.py decks/ni/figures/bunge/pfd_curves_L1_mesh.i

## 5. The fifteen specimens: FEM, PINN+FEM, PINN (report §7, §8.1 to §8.3)

    python3 scripts/ni_study.py    decks/ni/final/bunge/fem_15_h020.i all
    python3 scripts/ni_study.py    decks/ni/final/bunge/hybrid_15_h020.i all
    python3 scripts/ni_specimen.py decks/ni/final/bunge/pinn_base_h020.i deck run curve fracture

After a study, each case directory is post-processed with
`python3 scripts/ni_specimen.py specimen.i curve fracture`
(`scripts/tools/ni_bunge_chain.sh` runs the three stages and this step in
order on the compute host). Figures:

    python3 scripts/ni_final_figures.py          decks/ni/figures/bunge_final_fem.i
    python3 scripts/ni_final_specimen_figures.py decks/ni/figures/bunge_final_fem_specimens.i
    python3 scripts/ni_final_figures.py          decks/ni/figures/bunge_final_routes.i
    python3 scripts/ni_final_specimen_figures.py decks/ni/figures/bunge_final_hybrid_specimens.i
    python3 scripts/ni_pfd_curve_figures.py      decks/ni/figures/bunge/pfd_curves_L1_routes.i
    python3 scripts/ni_pfd_shape_figures.py      decks/ni/figures/bunge/pfd_shape_L1_routes.i

## 6. Surrogate (report §8.4)

    python3 scripts/ni_surrogate_design.py  decks/ni/final/bunge/surrogate_design_h020.i
    python3 scripts/ni_study.py results/ni_study/bunge/final/surrogate/data/study.i all
    python3 scripts/ni_surrogate_targets.py decks/ni/final/bunge/surrogate/targets_h020.i
    (cd decks/ni/final/bunge/surrogate && pinnacle run surrogate_h020.i)
    python3 scripts/ni_surrogate_figures.py decks/ni/figures/bunge_final_surrogate.i

The design step writes the thirty synthetic specimens, their meshes and the
study deck. `pinnacle run` resolves the deck's relative paths against the
working directory, so it is started from the deck's directory.

## Superseded material

Decks of the runs made before the Euler-convention fix are in
`archive/decks/ni_exchanged_angles/`, their figures in
`archive/figures/ni_exchanged_angles/` and `archive/figures/ni_first_attempt/`.
