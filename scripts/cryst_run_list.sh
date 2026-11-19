#!/bin/bash
# Run a list of decks, a few at once, and log exit code and wall time per deck.
# usage: scripts/cryst_run_list.sh <jobs> <logfile> deck1 deck2 ...   (run from the benchmark repo root)
jobs=$1; log=$2; shift 2
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}   # PINNACLE must be importable (pip install -e, or PYTHONPATH)
run_one() { d=$1; n=$(basename $d .i); t0=$(date +%s); python3 -W ignore -m pinnacle run $d > /dev/null 2> results/.err_$n.txt; rc=$?; echo "$n exit $rc wall $(( $(date +%s) - t0 )) s $(tail -1 results/.err_$n.txt | cut -c1-160)" >> $log; rm -f results/.err_$n.txt; }
export -f run_one; export log
printf '%s\n' "$@" | xargs -P $jobs -I{} bash -c 'run_one {}'
echo finished >> $log
