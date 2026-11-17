#!/bin/bash
# Al2219 recheck (2026-11-16), Part II CP+PFD stages after the Gc(l0) calibration, in order on
# the compute host. Run from the study root with the PINNACLE environment set; log in
# results/pfd_study/cp/bunge/chain.out.
set -u
for s in neural_h40 sweep_l_h40 grains_h40 shape_h40 sens_h40 sweep_h_l5; do
  python scripts/pfd_study.py decks/pfd_study/cp/bunge/$s.i all
  echo "$(date) $s done"
done
