#!/usr/bin/env python3
"""Run one refinement level of the Al2219 mesh-convergence study.

One level per process, writing a self-contained JSON record, because the levels
are independent and the wall time of the family is then the cost of its largest
member rather than the sum. Each record carries enough provenance -- code
revision, deck, mesh, timings, per-load-step history -- to be interpreted
without the process that produced it.

Global and local metrics are recorded together. A convergence study that tracks
only a global quantity (peak load, total reaction) will show convergence long
before the damage band is resolved, because the global response is an integral
that averages over the very region the mesh is failing to capture. The local
metrics here -- peak damage, damaged volume fraction, band width, and the
location of the damage maximum -- are the ones expected to converge last.

Usage:
    run_convergence_level.py <density> [--load-steps N] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent
PINNACLE = Path("/mnt/extended/projects/pinnacle")
VOID_CENTRE = np.array([25.0, 25.0, 25.0])
VOID_RADIUS_UM = 5.0
L0_UM = 1.0


def _git_rev(repo: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return "unknown"


def _find(obj, key, depth=0):
    if depth > 7:
        return None
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for value in obj.values():
            found = _find(value, key, depth + 1)
            if found is not None:
                return found
    return None


def local_metrics(nodes: np.ndarray, damage: np.ndarray) -> dict:
    """Damage-field measures, the ones that converge last."""
    if damage.size != nodes.shape[0]:
        return {"note": f"damage size {damage.size} != n_nodes {nodes.shape[0]}"}
    peak = float(damage.max())
    imax = int(np.argmax(damage))
    loc = nodes[imax]
    radial = float(np.linalg.norm(loc - VOID_CENTRE))
    # Band width along the ray through the damage maximum: the extent over
    # which damage exceeds half its peak. Reported in units of l0, since that
    # is what the phase-field regularization sets it by.
    hot = damage > 0.5 * peak if peak > 0 else np.zeros_like(damage, dtype=bool)
    if hot.any():
        r_hot = np.linalg.norm(nodes[hot] - VOID_CENTRE, axis=1)
        band = float(r_hot.max() - r_hot.min())
    else:
        band = float("nan")
    return {
        "peak_damage": peak,
        "mean_damage": float(damage.mean()),
        "damaged_node_fraction_gt_0p65": float((damage > 0.65).mean()),
        "damage_max_location_um": [float(x) for x in loc],
        "damage_max_radius_um": radial,
        "damage_max_distance_from_void_surface_um": radial - VOID_RADIUS_UM,
        "half_peak_band_width_um": band,
        "half_peak_band_width_over_l0": band / L0_UM,
    }


def main() -> int:
    # Declared before any use: `global` must precede the first reference to the
    # name in the function, and the argparse defaults below read STUDY. This is
    # a compile-time (symtable) error, not a parse error, so `ast.parse` accepts
    # the bad form and only `compile()` rejects it -- which is why the check
    # before dispatch is now compile(), and why three containers failed in one
    # second each on the first attempt.
    global STUDY, PINNACLE

    ap = argparse.ArgumentParser()
    ap.add_argument("density", type=int)
    ap.add_argument("--load-steps", type=int, default=20)
    # Alternations per load step. Exposed mainly so the whole harness path can
    # be exercised cheaply before a long dispatch; production runs leave it at
    # the deck value, since a low cap reports the coupling as unconverged
    # rather than silently lagging it (ADR 0024).
    ap.add_argument("--max-alternations", type=int, default=None)
    ap.add_argument("--out", default=str(STUDY / "results" / "convergence"))
    # Roots are arguments rather than constants so the same harness runs
    # unchanged locally and inside a remote container, where the study tree and
    # the framework are unpacked wherever the job staged them. A run whose
    # provenance says "unknown" for the code revision is not reproducible, so
    # --pinnacle-rev lets the caller pass the revision when the container has
    # no git metadata to read.
    ap.add_argument("--study-root", default=None)
    ap.add_argument("--pinnacle-root", default=None)
    ap.add_argument("--pinnacle-rev", default=None)
    args = ap.parse_args()

    if args.study_root:
        STUDY = Path(args.study_root).resolve()
    if args.pinnacle_root:
        PINNACLE = Path(args.pinnacle_root).resolve()

    sys.path.insert(0, str(PINNACLE / "src"))
    from pinnacle.io.loader import load_input
    from pinnacle.plugins import discover_components
    from pinnacle.runner.runner import Runner

    discover_components(force=True)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    d = args.density

    deck_src = (STUDY / "decks" / "al2219" / f"al2219_iso_001_d{d}.i").read_text()
    if args.load_steps != 20:
        deck_src = deck_src.replace("  n_steps        = 20", f"  n_steps        = {args.load_steps}")
    if args.max_alternations is not None:
        deck_src = deck_src.replace(
            "  max_iterations = 10", f"  max_iterations = {args.max_alternations}"
        )
    run_dir = out_dir / f"d{d}"
    deck_src = deck_src.replace(f"./runs/al2219_iso_001_d{d}", str(run_dir))
    deck = out_dir / f"al2219_iso_001_d{d}_run.i"
    deck.write_text(deck_src)

    mesh_path = (STUDY / "meshes" / "al2219" / f"void_r5_L50_d{d}.e").resolve()
    if not mesh_path.is_file():
        raise SystemExit(f"mesh not found: {mesh_path}")
    # The deck names the mesh by a path relative to the process cwd, so pin it
    # to the resolved absolute path: the harness must give the same answer
    # whatever directory it is launched from, locally or in a container.
    deck_src = deck_src.replace(
        f"path   = meshes/al2219/void_r5_L50_d{d}.e", f"path   = {mesh_path}"
    )

    mesh_rec = next(
        r for r in json.loads((STUDY / "meshes" / "al2219" / "mesh_family.json").read_text())
        if r["density"] == d
    )

    t0 = time.time()
    result = Runner.from_config(load_input(str(deck)))
    wall = time.time() - t0

    fields = result["fields"]
    # Nodes come from the MESH, not from ``fields["coordinates"]``. That field
    # carries a different grid on the composed path -- 216 points against the
    # mesh's 6,096 at d=16 -- so using it would have placed every local metric
    # at the wrong coordinates. It happened to raise on the reshape here; had
    # the two counts been compatible it would have silently reported the damage
    # maximum at a fictitious location. The assertion below keeps that failure
    # loud if either shape ever changes again.
    from pinnacle import registry

    mesh_obj = registry.create(
        "mesh", "file", path=str(mesh_path), format="exodus",
        scale=1.0, name_mapping={}, verify=False,
    )
    nodes = np.asarray(mesh_obj.nodes, dtype=float).reshape(-1, 3)
    damage = np.asarray(fields.get("damage", []), dtype=float).ravel()
    disp_raw = np.asarray(fields["displacement"], dtype=float).ravel()
    assert disp_raw.size == nodes.size, (
        f"displacement has {disp_raw.size} entries but the mesh has {nodes.size} "
        "coordinates; the field-to-mesh correspondence is not what this harness assumes"
    )
    disp = disp_raw.reshape(nodes.shape)
    stress = np.asarray(fields.get("stress", []), dtype=float)
    history = _find(result, "staggered") or []

    record = {
        "density": d,
        "n_elements": mesh_rec["n_elements"],
        "first_layer_um": mesh_rec["first_layer_um"],
        "points_across_l0": mesh_rec["points_across_l0"],
        "under_resolved_at_l0": bool(mesh_rec["points_across_l0"] < 4.0),
        "n_nodes": int(nodes.shape[0]),
        "n_dofs": int(3 * nodes.shape[0]),
        "load_steps": args.load_steps,
        "max_alternations": args.max_alternations,
        "wall_seconds": wall,
        "history": history,
        "all_steps_converged": bool(history) and all(h["converged"] for h in history),
        "global": {
            "max_displacement_um": float(np.abs(disp).max()),
            "mean_stress": float(np.mean(stress)) if stress.size else None,
            "max_stress": float(np.max(stress)) if stress.size else None,
        },
        "local": local_metrics(nodes, damage),
        "provenance": {
            "pinnacle_rev": args.pinnacle_rev or _git_rev(PINNACLE),
            "study_rev": _git_rev(STUDY),
            "deck": deck.name,
            "mesh": f"meshes/al2219/void_r5_L50_d{d}.e",
            "python": platform.python_version(),
            "precision": "float64",
            "backend": "numpy",
        },
    }
    dest = out_dir / f"level_d{d}.json"
    dest.write_text(json.dumps(record, indent=1))
    print(f"d={d} wall {wall/3600:.2f} h -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
