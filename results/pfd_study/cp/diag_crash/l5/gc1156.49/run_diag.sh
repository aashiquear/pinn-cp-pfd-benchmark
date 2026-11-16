#!/bin/bash
# Diagnostic rerun of the crashed implicit Gc trial (l0 5, Gc 1156.49): stderr kept,
# faulthandler on, exit status and peak RSS recorded, RSS sampled every 60 s.
cd "$(dirname "$0")"
source ~/venvs/pinnacle/bin/activate
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONFAULTHANDLER=1
/usr/bin/time -v pinnacle run pfd.i > stdout.log 2> stderr.log &
pid=$!
while kill -0 $pid 2>/dev/null; do
  c=$(pgrep -P $pid | head -1)
  [ -n "$c" ] && echo "$(date +%T) $(awk "/VmRSS/{print \$2}" /proc/$c/status 2>/dev/null)" >> rss.log
  sleep 60
done
wait $pid; echo "exit $?" > exit_status.txt
