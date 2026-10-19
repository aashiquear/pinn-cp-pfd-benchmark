# Reproduce: Ni single-crystal dog-bone specimens (Part III)

Commands that rebuild every Part III result from the versioned decks. Run from
the study root with PINNACLE `feat/pinn_benchmark` (f83c378 or later) and
hexaform `feat/plate-hole` (1cb1671 or later) installed; on the compute host
`source ~/venvs/pinnacle/bin/activate`. Units are mm, MPa, s.

## 1. Meshes (hexaform `dogbone`)

    python3 scripts/ni_specimen_mesh.py decks/ni/mesh/h040.i   # also h020.i, h010.i, h020z2.i

Writes `results/ni_study/mesh/<level>/<specimen>/mesh.e` (+ `mesh.json` with
counts and quality) for the 15 specimens of Table 2 and `figures/ni/mesh_<level>.png`.

## 2. One specimen (CP only)

    python3 scripts/ni_specimen.py decks/ni/cp/calibrated_L1_h040.i all

Stages `deck` (renders `decks/ni/templates/cp_specimen_template.i` into
`<out_dir>/cp.i`), `run` (`pinnacle run cp.i`, threads from the deck) and
`curve` (`results.csv`: gauge strain, grip strain, engineering stress, Newton
count, cutbacks; `summary.json`).

## 3. Calibration (Gauss-Newton, docs/12 §2)

    python3 scripts/ni_study.py decks/ni/cp/cal_gn1.i prepare
    python3 scripts/ni_study.py decks/ni/cp/cal_gn1.i run       # 5 runs side by side
    python3 scripts/ni_study.py decks/ni/cp/cal_gn1.i collect
    python3 scripts/ni_gauss_newton.py decks/ni/cp/cal_gn1.i    # writes cal_center2.i, cal_gn2.i

repeated for `cal_gn2.i` .. `cal_gn4.i`. The centre of iteration 4 is the
calibrated set, `decks/ni/cp/calibrated_L1_h040.i`.

## 4. Figures

    python3 scripts/ni_calibration_figures.py decks/ni/figures/cal_figures.i

Writes `figures/ni/cal/cal_*.png` (report Part III §4).

## CP+PFD (Part III stage 6)

A specimen deck that sets `Gc` (and `l0`) is a CP+PFD case: `scripts/ni_specimen.py DECK
deck run curve fracture` renders `decks/ni/templates/cp_pfd_specimen_template.i` into
`<out_dir>/pfd.i`, runs it (restart checkpoints every `checkpoint_every` steps; the newest
is kept when the run dies, so `curve` still has its history), writes `results.csv` (one row
per staggered load step: gauge strain, engineering stress, max damage, nodes above the
0.8 threshold, alternations) with `summary.json`, and `fracture` applies the paper's failure
criterion to the damage snapshots: the first snapshot in which the nodes with d > 0.85 cut
the specimen in two (no path of intact nodes from the x- to the x+ grip face), with the
crack's location along the gage and its angle to the loading axis (`fracture.json`).
Example: `decks/ni/pfd/check_L1_h020.i`.
