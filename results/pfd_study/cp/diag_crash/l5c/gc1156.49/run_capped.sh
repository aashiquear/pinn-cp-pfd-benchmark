#!/bin/bash
# Diagnostic (PINNACLE e281696): the crashed implicit trial l0 5, Gc 1156.49 from the
# start under a 16 GB address-space cap, so a runaway allocation in the
# backward-Euler fallback raises MemoryError with a traceback instead of
# exhausting the host; restart checkpoints every 4 steps; RSS per minute.
cd "$(dirname "$0")"
source ~/venvs/pinnacle/bin/activate
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONFAULTHANDLER=1
ulimit -v 16000000
pinnacle run pfd.i > stdout.log 2> stderr.log &
pid=$!
while kill -0 $pid 2>/dev/null; do echo "$(date +%T) $(awk "/VmRSS/{print \$2}" /proc/$pid/status 2>/dev/null)" >> rss.log; sleep 60; done
wait $pid; echo "exit $?" > exit_status.txt
