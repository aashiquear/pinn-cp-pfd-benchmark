#!/bin/bash
# Diagnostic: resume the crashed trial from step 228 under a 16 GB address-space cap,
# so the fallback's large allocation raises MemoryError with a traceback.
cd "$(dirname "$0")"
source ~/venvs/pinnacle/bin/activate
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONFAULTHANDLER=1
ulimit -v 16000000
pinnacle run pfd.i > stdout.log 2> stderr.log
echo "exit $?" > exit_status.txt
