#!/usr/bin/env python3
"""Run one pseudo-3D plate case and write the standard outputs.

Every plate run in the study goes through here, so that no result can exist
without its ``results.csv`` and full-field VTK series. The G20 constitutive
strain-increment limit is derived from the loading rather than passed in, so
a caller cannot under-resolve a run by forgetting a flag.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import warnings
from pathlib import Path

import sys

import numpy as np

warnings.filterwarnings("ignore")
# The harness is run by path from the study root, so its own directory is not
# necessarily on sys.path; the shared results writer lives beside it.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pinnacle.plugins import discover_components  # noqa: E402

discover_components(force=True)
from pinnacle.io.loader import load_input  # noqa: E402
from pinnacle.runner.runner import Runner  # noqa: E402

from run_results import write_run_results  # noqa: E402

STUDY = Path(__file__).resolve().parent.parent
MAX_STRAIN_INCREMENT = 8.0e-4      # per constitutive update; study gap G20
PLATE_LENGTH = 100.0               # um, the loaded (x) dimension

# Calibrated Al2219 CP set (plain-plate refit supersedes this in Phase 1).
CAL_AL2219 = {
    "TAUC0": "60.31771394231891", "RHO0": "2.393387495846754",
    "K1": "0.08114503403820844", "YC2": "0.0020581167904802733",
}


def build_deck(*, mesh: Path, orient: str, gc: float, l0: float, nsteps: int,
               mechsub: int, maxalt: int, strain: float, out_dir: Path,
               cal: dict[str, str]) -> Path:
    src = (STUDY / "decks/al2219/plate_template.i").read_text()
    sub = dict(cal, MESH=str(mesh), ORIENT=orient, GC=f"{gc:.9g}", L0=f"{l0:.9g}",
               NSTEPS=str(nsteps), MECHSUB=str(mechsub), MAXALT=str(maxalt),
               STRAIN=f"{strain:.9g}", OUTDIR=str(out_dir / "solver"))
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    leftover = set(re.findall(r"@([A-Z0-9_]+)@", re.sub(r"#.*", "", src)))
    if leftover:
        raise SystemExit(f"unsubstituted placeholders: {sorted(leftover)}")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "deck.i"
    path.write_text(src)
    return path


def _find(obj, key, depth=0):
    if depth > 7:
        return None
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            got = _find(v, key, depth + 1)
            if got is not None:
                return got
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--density", type=int, required=True)
    ap.add_argument("--geometry", choices=("hole", "plain"), default="hole")
    ap.add_argument("--orientation", default="[1.0, 0.0, 0.0]")
    ap.add_argument("--gc", type=float, default=10.0)
    ap.add_argument("--l0", type=float, default=1.0)
    ap.add_argument("--strain", type=float, default=0.12)
    # Denser by default than the smoke runs: at 20 steps to 12% strain the
    # first recorded point is already at 0.6%, past yield, so the elastic
    # branch is one point wide and fracture initiation is bracketed to a
    # whole load step. 60 steps puts the first point at 0.2%. Load steps
    # multiply damage solves but NOT constitutive work, since the sub-step
    # count per step falls as the increment shrinks.
    ap.add_argument("--load-steps", type=int, default=60)
    ap.add_argument("--alternations", type=int, default=2)
    ap.add_argument("--no-vtk", action="store_true")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", default="results/plate")
    args = ap.parse_args()

    # G20: the constitutive increment per update must stay at or below 8e-4.
    # Derived, not passed, so a run cannot be silently under-resolved.
    d_eps = args.strain / max(1, args.load_steps)
    mechsub = max(1, int(np.ceil(d_eps / MAX_STRAIN_INCREMENT)))

    tag = args.tag or f"{args.geometry}_d{args.density}_Gc{args.gc:g}_l0{args.l0:g}"
    out_dir = STUDY / args.out / tag
    mesh_npz = STUDY / f"meshes/plate/plate_{args.geometry}_d{args.density}.npz"
    mesh_e = STUDY / f"meshes/plate/plate_{args.geometry}_d{args.density}.e"
    z = np.load(mesh_npz)
    nodes, cells = z["nodes"], z["cells"]

    deck = build_deck(mesh=mesh_e, orient=args.orientation, gc=args.gc, l0=args.l0,
                      nsteps=args.load_steps, mechsub=mechsub,
                      maxalt=args.alternations, strain=args.strain,
                      out_dir=out_dir, cal=CAL_AL2219)

    t0 = time.perf_counter()
    result = Runner.from_config(load_input(str(deck)))
    wall = time.perf_counter() - t0

    snaps = _find(result, "snapshots") or []
    summary = write_run_results(
        out_dir, snaps, nodes, cells,
        boundary_nodes=np.asarray(z["set_xmax"], dtype=int),
        gauge_length=PLATE_LENGTH, load_axis=0, write_vtk=not args.no_vtk,
    )
    rec = {
        "tag": tag, "geometry": args.geometry, "density": args.density,
        "n_elements": int(cells.shape[0]), "n_nodes": int(nodes.shape[0]),
        "orientation": args.orientation, "Gc_MPa_um": args.gc, "l0_um": args.l0,
        "target_strain": args.strain, "load_steps": args.load_steps,
        "mech_substeps": mechsub, "alternations": args.alternations,
        "strain_increment_per_update": d_eps / mechsub,
        "wall_seconds": wall,
        "vtk_series": str(out_dir / "fields" / "fields.pvd"),
        **summary,
    }
    (out_dir / "run.json").write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: rec[k] for k in
                      ("tag", "n_elements", "mech_substeps", "n_steps",
                       "vtk_files", "wall_seconds")}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
