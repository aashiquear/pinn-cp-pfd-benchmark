#!/usr/bin/env python3
"""One coupled void run at a prescribed (Gc, l0, mesh density), for calibration.

CALIBRATED CP SET. Uses the parameters fitted against the Fig. 3 experiment,
not the published ones. In this model the published set gives tau_c(0) =
176.71 MPa against the paper's reported 43.72 -- a factor sqrt(1+11q) = 4.05,
because PINNACLE forms the Taylor term from the latent-hardened forest sum --
so it over-predicts stress about threefold and the elastic driving energy about
ninefold. Calibrating Gc against that inflated driving force would fit the
error, not the material.

WHAT IS BEING CALIBRATED AGAINST. The source has no experimental
fracture-initiation strain for Al2219; Fig. 3 is a tensile curve with no
fracture, and the Fig. 8(a) initiation strains are the paper's own SIMULATION
output. So this script does not assert a single "correct" Gc. It measures the
Gc -> initiation-strain RELATION at fixed h and l0, which is the calibration
curve; selecting a point on it against the paper's reported 7.0% for isochoric
[001] is a simulation-to-simulation comparison and is labeled as such wherever
reported.

INITIATION CRITERION. The paper's own: the strain at which the damage variable
first reaches 0.65 anywhere in the domain (Sec. 3.2). Recorded here by load
step, so the reported strain carries the load-step resolution as its
uncertainty rather than being quoted as exact.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent
PINNACLE = Path("/mnt/extended/projects/pinnacle")

# Fitted against the nine Vaughan (2024) markers; see results/calibration/.
CALIBRATED_CP = {
    "tau_c0": 60.317714,
    "rho0": 2.393387,
    "k1": 0.081145,
    "annihilation_diameter": 0.002058,
}
D_THRESHOLD = 0.65
VOID_CENTRE = np.array([25.0, 25.0, 25.0])
VOID_RADIUS_UM = 5.0


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


def _sub(text: str, key: str, value, width: int = 21) -> str:
    """Replace `key = ...` on its own line, preserving the deck's alignment."""
    out = []
    hit = False
    for line in text.split("\n"):
        stripped = line.strip()
        if stripped.startswith(f"{key} ") or stripped.startswith(f"{key}="):
            if "=" in stripped and stripped.split("=")[0].strip() == key:
                indent = line[: len(line) - len(line.lstrip())]
                out.append(f"{indent}{key:<{width}} = {value}")
                hit = True
                continue
        out.append(line)
    if not hit:
        raise SystemExit(f"key {key!r} not found in deck; cannot set it silently")
    return "\n".join(out)


def build_deck(density: int, gc: float, l0: float, load_steps: int,
               max_alternations: int, strain: float, out_dir: Path) -> Path:
    src = (STUDY / "decks" / "al2219" / f"al2219_iso_001_d{density}.i").read_text()
    mesh = (STUDY / "meshes" / "al2219" / f"void_r5_L50_d{density}.e").resolve()
    if not mesh.is_file():
        raise SystemExit(f"mesh not found: {mesh}")
    src = src.replace(f"path   = meshes/al2219/void_r5_L50_d{density}.e",
                      f"path   = {mesh}")
    for key, val in CALIBRATED_CP.items():
        src = _sub(src, key, f"{val:.9g}", width=21)
    src = _sub(src, "Gc", f"{gc:.9g}", width=12)
    src = _sub(src, "length_scale", f"{l0:.9g}", width=12)
    src = _sub(src, "n_steps", str(load_steps), width=14)
    src = _sub(src, "max_iterations", str(max_alternations), width=14)
    src = src.replace("./runs/al2219_iso_001_d{}".format(density),
                      f"./runs/gc_d{density}")
    deck = out_dir / f"gc_d{density}_Gc{gc:g}_l0{l0:g}.i"
    deck.write_text(src)
    return deck


def main() -> int:
    # Declared before ANY use: `global` must precede the first reference to the
    # name in the function, and the argparse defaults below read STUDY. This is
    # a compile-time symtable error, not a parse error, so `ast.parse` accepts
    # it and only `compile()` rejects it -- which is why the pre-dispatch check
    # is compile().
    global STUDY, PINNACLE

    ap = argparse.ArgumentParser()
    ap.add_argument("density", type=int)
    ap.add_argument("--gc", type=float, required=True)
    ap.add_argument("--l0", type=float, default=1.0)
    ap.add_argument("--load-steps", type=int, default=20)
    ap.add_argument("--max-alternations", type=int, default=3)
    ap.add_argument("--strain", type=float, default=0.12)
    ap.add_argument("--out", default=str(STUDY / "results" / "gc_calibration"))
    ap.add_argument("--study-root", default=None)
    ap.add_argument("--pinnacle-root", default=None)
    ap.add_argument("--pinnacle-rev", default=None)
    args = ap.parse_args()

    if args.study_root:
        STUDY = Path(args.study_root).resolve()
    if args.pinnacle_root:
        PINNACLE = Path(args.pinnacle_root).resolve()
    sys.path.insert(0, str(PINNACLE / "src"))
    sys.path.insert(0, str(Path(__file__).resolve().parent))

    from pinnacle import registry
    from pinnacle.io.loader import load_input
    from pinnacle.plugins import discover_components
    from pinnacle.runner.runner import Runner

    discover_components(force=True)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    deck = build_deck(args.density, args.gc, args.l0, args.load_steps,
                      args.max_alternations, args.strain, out_dir)

    t0 = time.time()
    result = Runner.from_config(load_input(str(deck)))
    wall = time.time() - t0

    mesh_path = (STUDY / "meshes" / "al2219"
                 / f"void_r5_L50_d{args.density}.e").resolve()
    mesh_obj = registry.create("mesh", "file", path=str(mesh_path),
                               format="exodus", scale=1.0, name_mapping={},
                               verify=False)
    nodes = np.asarray(mesh_obj.nodes, dtype=float).reshape(-1, 3)

    snapshots = _find(result, "snapshots") or []
    history = _find(result, "staggered") or []

    # Initiation: FIRST load step at which d reaches the criterion anywhere.
    init_step = None
    per_step = []
    for snap in snapshots:
        d = np.asarray(snap["damage"], dtype=float).ravel()
        j = int(np.argmax(d))
        r = float(np.linalg.norm(nodes[j] - VOID_CENTRE))
        per_step.append({
            "step": int(snap["step"]),
            "load_factor": float(snap["step_factor"]),
            "strain": float(snap["step_factor"]) * args.strain,
            "peak_damage": float(d.max()),
            "mean_damage": float(d.mean()),
            "fraction_above_threshold": float((d >= D_THRESHOLD).mean()),
            "argmax_xyz_um": nodes[j].tolist(),
            "argmax_distance_from_void_surface_um": r - VOID_RADIUS_UM,
        })
        if init_step is None and d.max() >= D_THRESHOLD:
            init_step = per_step[-1]

    step_strain = args.strain / max(1, args.load_steps)
    record = {
        "density": args.density,
        "Gc_MPa_um": args.gc,
        "l0_um": args.l0,
        "load_steps": args.load_steps,
        "max_alternations": args.max_alternations,
        "applied_strain": args.strain,
        "wall_seconds": wall,
        "n_elements": int(len(mesh_obj.connectivity)),
        "n_nodes": int(nodes.shape[0]),
        "cp_parameters": CALIBRATED_CP,
        "cp_source": "fitted against the nine Vaughan (2024) Fig. 3 markers",
        "initiation_criterion_d": D_THRESHOLD,
        "initiation_strain": init_step["strain"] if init_step else None,
        "initiation_strain_resolution": step_strain,
        "initiation_step": init_step["step"] if init_step else None,
        "initiation_at_void_surface_um": (
            init_step["argmax_distance_from_void_surface_um"] if init_step else None),
        "initiated": init_step is not None,
        "per_step": per_step,
        "all_steps_converged": all(h.get("converged") for h in history) if history else None,
        "alternations_per_step": [int(h.get("iterations", 0)) for h in history],
    }

    if snapshots:
        from damage_evidence import write_evidence

        try:
            record["evidence"] = write_evidence(
                snapshots, nodes, out_dir / f"evidence_d{args.density}_Gc{args.gc:g}",
                label=f"Al2219 d={args.density}, Gc={args.gc:g} MPa.um, l0={args.l0:g} um",
                void_centre=VOID_CENTRE, void_radius=VOID_RADIUS_UM)
        except Exception as exc:
            record["evidence"] = {"error": f"{type(exc).__name__}: {exc}"}
    else:
        record["evidence"] = {"error": "no snapshots recorded"}

    rev = args.pinnacle_rev
    if rev is None:
        try:
            rev = subprocess.run(["git", "-C", str(PINNACLE), "rev-parse", "--short", "HEAD"],
                                 capture_output=True, text=True, timeout=20).stdout.strip()
        except Exception:
            rev = "unknown"
    record["provenance"] = {"pinnacle_rev": rev, "deck": deck.name,
                            "mesh": mesh_path.name, "backend": "numpy",
                            "precision": "float64"}

    dest = out_dir / f"gc_d{args.density}_Gc{args.gc:g}_l0{args.l0:g}.json"
    dest.write_text(json.dumps(record, indent=1))
    print(f"d={args.density} Gc={args.gc:g} l0={args.l0:g}: wall {wall/60:.1f} min, "
          f"initiated={record['initiated']}, "
          f"eps_init={record['initiation_strain']} -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
