#!/usr/bin/env bash
# One Al2219 Gc-calibration run inside a Modal container.
#
# Inputs stage flat into the workdir root, so this rebuilds the layout the
# harness expects. PINNACLE ships as a tarball rather than an image layer: the
# repository is private, so the container cannot clone it, and staging the
# source keeps the running revision pinned to what was measured locally.
set -euo pipefail

DENSITY="${DENSITY:?}"
GC="${GC:?}"
L0="${L0:-1.0}"
LOAD_STEPS="${LOAD_STEPS:-20}"
MAX_ALT="${MAX_ALT:-3}"
PINNACLE_REV="${PINNACLE_REV:-unknown}"

mkdir -p ~/projects out
tar -xzf pinnacle_complete.tar.gz -C ~/projects

mkdir -p study/scripts study/decks/al2219 study/meshes/al2219
mv al2219_iso_001_d*.i        study/decks/al2219/  2>/dev/null || true
mv void_r5_L50_d*.e           study/meshes/al2219/ 2>/dev/null || true
mv run_gc_calibration.py damage_evidence.py study/scripts/

python -c "import numpy, scipy, meshio, netCDF4, matplotlib; print('deps ok', numpy.__version__)"
nproc

python study/scripts/run_gc_calibration.py "${DENSITY}" \
    --gc "${GC}" \
    --l0 "${L0}" \
    --load-steps "${LOAD_STEPS}" \
    --max-alternations "${MAX_ALT}" \
    --study-root study \
    --pinnacle-root ~/projects/pinnacle \
    --pinnacle-rev "${PINNACLE_REV}" \
    --out out

# Flatten the evidence directory so the named outputs transfer back.
find out -type f | sed 's|^|  produced: |'
for f in $(find out -name 'damage_evolution.png' -o -name 'damage_evolution.gif' -o -name 'damage_history.csv'); do
    base=$(basename "$f")
    cp "$f" "out/${DENSITY}_Gc${GC}_${base}"
done
ls out/
