#!/bin/bash
# Al2219 recheck (2026-11-16), Part I stages that follow the recalibration, in order on the
# compute host: parameter sensitivity of the PINN+FEM (linear) route, the 24 FE samples of the
# surrogate, parameter sensitivity of the grain-aware PINN-as-solver route.
# Run from the study root with the PINNACLE environment set; log in results/cp_study/bunge/chain.out.
set -u
python scripts/cp_parameter_study.py decks/cp_study/bunge/sens_h40_linear.i all
echo "$(date) stage 1 done"
python scripts/cp_parameter_study.py results/cp_study/bunge/surrogate/h40_data/study.i all
echo "$(date) stage 2 done"
python scripts/cp_parameter_study.py decks/cp_study/bunge/sens_h40_grain_solver.i all
echo "$(date) stage 3 done"
