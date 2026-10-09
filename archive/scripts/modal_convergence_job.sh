#!/usr/bin/env bash
# Run one Al2219 convergence level inside a Modal container.
#
# Inputs stage flat into the workdir root, so this rebuilds the directory
# layout the harness expects before invoking it. PINNACLE is staged as a
# tarball rather than installed in the image: the repository is private, so a
# sandbox cannot clone it, and staging per job keeps one image usable across
# framework revisions.
#
# Environment: DENSITY, LOAD_STEPS, PINNACLE_REV.
set -euo pipefail

DENSITY="${DENSITY:?DENSITY not set}"
LOAD_STEPS="${LOAD_STEPS:-20}"
PINNACLE_REV="${PINNACLE_REV:-unknown}"

mkdir -p ~/projects
tar xzf pinnacle_complete.tar.gz -C ~/projects

# PINNACLE is added to sys.path by the harness rather than pip-installed: the
# image already carries every runtime dependency, and an editable install here
# would spend container time re-resolving them for no behavioral difference.
mkdir -p study/decks/al2219 study/meshes/al2219 study/scripts out
cp "al2219_iso_001_d${DENSITY}.i" study/decks/al2219/
cp "void_r5_L50_d${DENSITY}.e"    study/meshes/al2219/
cp mesh_family.json               study/meshes/al2219/
cp run_convergence_level.py       study/scripts/

python -c "import numpy, scipy, meshio, netCDF4; print('deps ok', numpy.__version__, scipy.__version__)"
nproc  # recorded so the wall time in the record can be read against the core count

python study/scripts/run_convergence_level.py "${DENSITY}" \
    --load-steps "${LOAD_STEPS}" \
    --study-root study \
    --pinnacle-root ~/projects/pinnacle \
    --pinnacle-rev "${PINNACLE_REV}" \
    --out out
