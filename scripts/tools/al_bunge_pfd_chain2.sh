#!/bin/bash
# Al2219 recheck (2026-11-16), last host stages, in order: the three neural-route cases of the Part I
# hardware benchmark on the torch backend, the Part II hardware benchmark (both alone on the host),
# the 24 CP+PFD surrogate samples, their collection and the two surrogate trainings.
# Run from the study root with the PINNACLE environment set; log in results/pfd_study/cp/bunge/chain2.out.
set -u
python scripts/hardware_benchmark.py decks/benchmark/bunge/hardware_h40.i hybrid_linear_t8 solver_grain_gpu solver_grain_cpu
echo "$(date) hardware part I done"
python scripts/hardware_benchmark.py decks/benchmark/bunge/hardware_pfd_h40.i
echo "$(date) hardware part II done"
python scripts/pfd_study.py results/pfd_study/cp/bunge/surrogate/h40_data/study.i all
echo "$(date) surrogate samples done"
python scripts/surrogate_design.py decks/pfd_study/cp/bunge/surrogate/design_h40.i collect
echo "$(date) surrogate collect done"
for t in curve field; do
  OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 pinnacle run decks/pfd_study/cp/bunge/surrogate/plate_h40_surrogate_$t.i > results/pfd_study/cp/bunge/surrogate/h40_$t.log 2>&1
  echo "$(date) surrogate $t done exit $?"
done
