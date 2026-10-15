#!/bin/bash
# J2+PFD at W_c 200 (docs/10 §4.10): after the W_c 400 bisection (PID 2607143) ends,
# the Gc(l0) calibration at l0 10/20 and the 4-value sensitivity in parallel, then the
# l0 sweep and the mesh study (h20, h80) in parallel.
cd ~/projects/pinn-cp-pfd-benchmark && source ~/venvs/pinnacle/bin/activate
export TORCH_BLAS_PREFER_HIPBLASLT=0 TORCH_ROCM_AOTRITON_ENABLE_EXPERIMENTAL=1
while kill -0 2607143 2>/dev/null; do sleep 60; done
python3 scripts/pfd_calibrate.py decks/pfd_study/j2/cal_h40_wc200.i > results/pfd_study/j2/cal_h40_wc200_l10_l20.log 2>&1 &
P1=$!
python3 scripts/pfd_study.py decks/pfd_study/j2/sens_h40.i all > results/pfd_study/j2/sens_h40_wc200.log 2>&1 &
P2=$!
wait $P1 $P2
python3 scripts/pfd_study.py decks/pfd_study/j2/sweep_l_h40.i all > results/pfd_study/j2/sweep_l_h40_wc200.log 2>&1 &
P3=$!
python3 scripts/pfd_study.py decks/pfd_study/j2/sweep_h_l5.i all > results/pfd_study/j2/sweep_h_l5_wc200.log 2>&1 &
P4=$!
wait $P3 $P4
echo chain done $(date)
