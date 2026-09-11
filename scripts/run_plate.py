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
               cal: dict[str, str], maxcutbacks: int = 15) -> Path:
    src = (STUDY / "decks/al2219/plate_template.i").read_text()
    sub = dict(cal, MESH=str(mesh), ORIENT=orient, GC=f"{gc:.9g}", L0=f"{l0:.9g}",
               NSTEPS=str(nsteps), MECHSUB=str(mechsub), MAXALT=str(maxalt),
               MAXCUTBACKS=str(maxcutbacks),
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
    ap.add_argument("--max-cutbacks", type=int, default=15,
                    help="Load-increment halvings before a step is abandoned. "
                         "The framework default of 6 bottoms out at 1/64 of a "
                         "load step, which is too coarse at the crack tip.")
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
                      out_dir=out_dir, cal=CAL_AL2219,
                      maxcutbacks=args.max_cutbacks)

    t0 = time.perf_counter()
    result = Runner.from_config(load_input(str(deck)))
    wall = time.perf_counter() - t0

    # Diagnostics kept for BOTH outcomes. A failed run is when the
    # alternation history is most wanted, and re-running to obtain it costs
    # what the run cost. ``converged``/``failed_steps`` are the top-level
    # flags the mechanics executioner now reports; ``staggered`` is the
    # per-load-step alternation record.
    diag = {
        "staggered": _find(result, "staggered"),
        "mech_converged": _find(result, "converged"),
        "mech_failed_steps": _find(result, "failed_steps"),
    }
    snaps = _find(result, "snapshots") or []
    summary = write_run_results(
        out_dir, snaps, nodes, cells,
        boundary_nodes=np.asarray(z["set_xmax"], dtype=int),
        gauge_length=PLATE_LENGTH, load_axis=0, write_vtk=not args.no_vtk,
    )
    _early = dict(diag)
    _early.update({"tag": tag, "n_elements": int(n_elem), "wall_seconds": wall,
                   "status": "pending_validation"})
    (out_dir / "run.json").write_text(json.dumps(_early, indent=1, default=str))

    # VERIFY THE BOUNDARY CONDITION WAS ACTUALLY APPLIED (study gap G25).
    #
    # Checked against the physical invariant rather than a status flag. A
    # mechanics step that hits the cutback floor is recorded
    # ``converged: False`` and then ABANDONED -- step_idx advances, the load
    # increment is never solved, and the run returns success at the top level
    # with a displacement field that does not satisfy its own Dirichlet
    # condition. Guarding on the reported flag does not work: with a single
    # alternation the staggered record is ``converged: False`` by
    # construction, so it cannot distinguish a lagged coupling from an
    # abandoned load step. The prescribed displacement is unambiguous.
    import csv as _csv

    with (out_dir / "results.csv").open() as fh:
        _rows = list(_csv.DictReader(fh))
    bad = []
    for _r in _rows:
        want = args.strain * float(_r["load_factor"])
        got = float(_r["boundary_strain"])
        if abs(got - want) > max(1e-6, 1e-3 * abs(want)):
            bad.append(int(_r["step"]))

    # Distinguish a run that FAILED from one that FINISHED FRACTURING.
    #
    # Once the specimen is severed it carries no load and is a mechanism:
    # there is no equilibrium state at a larger prescribed displacement, so
    # the solve legitimately has nowhere to go. That is a terminal condition
    # of the physics, not a defect, and the load steps before it are a valid
    # fracture curve that must not be thrown away. The discriminator is the
    # load-carrying capacity at the last good step, measured against the peak
    # the run reached: a genuine failure aborts while the specimen is still
    # carrying load.
    status = "complete"
    valid = _rows
    if bad:
        first_bad = min(bad)
        valid = [r for r in _rows if int(r["step"]) < first_bad]
        peak = max((abs(float(r["avg_stress_axial"])) for r in _rows), default=0.0)
        last_ok = abs(float(valid[-1]["avg_stress_axial"])) if valid else float("inf")
        severed = bool(valid) and peak > 0 and last_ok <= 0.05 * peak
        if severed:
            status = "fractured"
            with (out_dir / "results.csv").open("w", newline="") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(valid[0].keys()))
                w.writeheader()
                w.writerows(valid)
            (out_dir / "FRACTURED.txt").write_text(
                f"Run terminated at complete fracture after step {first_bad - 1}.\n"
                f"Axial stress fell to {last_ok:.4g} MPa from a peak of {peak:.4g} MPa "
                f"({100 * last_ok / peak:.2f}%).\n"
                f"Steps {first_bad}-{int(_rows[-1]['step'])} carry no equilibrium "
                "solution and were dropped; the retained steps are a valid "
                "fracture curve.\n"
            )
        else:
            (out_dir / "INVALID.txt").write_text(
                f"Run failed at step {first_bad} while still carrying "
                f"{last_ok:.4g} MPa ({100 * last_ok / max(peak, 1e-30):.1f}% of peak "
                f"{peak:.4g} MPa), i.e. before fracture completed.\n"
            )
            raise SystemExit(
                f"{tag}: the prescribed boundary displacement was not applied from "
                f"step {first_bad}, while the specimen was still carrying "
                f"{last_ok:.4g} MPa ({100 * last_ok / max(peak, 1e-30):.1f}% of peak). "
                "This is a failed solve, not completed fracture. See G25."
            )

    rec = {
        "tag": tag, "geometry": args.geometry, "density": args.density,
        "n_elements": int(cells.shape[0]), "n_nodes": int(nodes.shape[0]),
        "orientation": args.orientation, "Gc_MPa_um": args.gc, "l0_um": args.l0,
        "target_strain": args.strain, "load_steps": args.load_steps,
        "mech_substeps": mechsub, "alternations": args.alternations,
        "strain_increment_per_update": d_eps / mechsub,
        "wall_seconds": wall,
        "diagnostics": diag,
        "status": status,
        "steps_retained": len(valid),
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
