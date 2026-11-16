#!/bin/bash
# Host chain of the Ni rerun after the Euler-convention fix (PINNACLE ADR 0036, 2026-11-16).
# Run from the study root on the ROCm host, after decks/ni/final/bunge/fem_15_h020.i was started:
#   nohup setsid bash scripts/tools/ni_bunge_chain.sh > results/ni_study/bunge/chain.out 2>&1 &
# Stage 1 waits for the 15 FEM runs and post-processes them; stage 2 runs PINN+FEM on the 15
# specimens and PINN on L.1; stage 3 runs the 30 synthetic specimens of the surrogate design and the
# h010 mesh check of L.1. Every run is a deck; rows in docs/run_ledger.md (ni-bunge-*).
export PATH=$HOME/venvs/pinnacle/bin:$PATH PYTHONPATH=$HOME/projects/pinnacle_bunge/src PYTHONUNBUFFERED=1
S=$HOME/projects/pinn-cp-pfd-benchmark
N=results/ni_study/bunge
cd "$S" || exit 1
post() { for d in "$1"/*/; do [ -f "$d/specimen.i" ] && (cd "$d" && python "$S/scripts/ni_specimen.py" specimen.i curve fracture > post.log 2>&1); done; }
until grep -q "summary.csv" $N/fem_15.out; do sleep 60; done
post $N/final/fem; echo "stage 1 done $(date)"
python scripts/ni_study.py decks/ni/final/bunge/hybrid_15_h020.i all > $N/hybrid_15.out 2>&1 &
python scripts/ni_specimen.py decks/ni/final/bunge/pinn_base_h020.i deck run curve fracture > $N/pinn_L1.out 2>&1 &
wait
post $N/final/hybrid; echo "stage 2 done $(date)"
python scripts/ni_study.py $N/final/surrogate/data/study.i all > $N/surrogate_fem_30.out 2>&1 &
python scripts/ni_specimen.py decks/ni/pfd/bunge/mesh_L1_h010.i deck run curve fracture > $N/mesh_h010.out 2>&1 &
wait
post $N/final/surrogate/data/runs; echo "stage 3 done $(date)"
