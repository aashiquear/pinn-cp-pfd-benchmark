#!/usr/bin/env python3
"""Part V: run one PINNACLE process deck and keep its fields and field series.

usage: scripts/lpbf_run.py <deck.i> [key=value ...]

Runs the deck in this process (PINNACLE on PYTHONPATH) and writes, into its output_dir:

  fields.npz   nodes; the final nodal values of every variable; for decks with
               `snapshots_every` in [Executioner], `snapshot_times` and one stack
               `snapshot_<variable>` (frames x nodes x components, float32): the stack
               the animations are drawn from (scripts/lpbf_animate.py)
  history.csv  per step: step, time, dt, the largest value of every variable, and the
               values the services report
  wall.json    wall time, steps, mesh size, overrides

`key=value` replaces a deck key by its dotted path, e.g. executioner.end_time=2e-4.

With `global.ranks=N` or `global.rank_devices="cuda:0 cuda:1 ..."` (in the deck or as an
override) the mesh is cut into parts and every part runs in its own process on its own
device (PINNACLE ADR 0057, `pinnacle.parallel.run_decomposed`). The files are the same;
history.csv then holds the services of rank 0 only (one part of the mesh), and wall.json
lists the ranks with their devices, sizes, exchanges and wall times.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    from pinnacle.io.loader import load_input
    from pinnacle.plugins import discover_components
    from pinnacle.runner import Runner

    deck = Path(sys.argv[1]).resolve()
    discover_components(force=True)
    cfg = load_input(str(deck))
    applied = {}
    for item in sys.argv[2:]:
        dotted, value = item.split("=", 1)
        node = cfg
        parts = dotted.split(".")
        for part in parts[:-1]:
            node = node[part]
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError:
            parsed = value
        node[parts[-1]] = parsed
        applied[dotted] = parsed
    raw = Path(str(cfg["global"]["output_dir"]))
    target = raw if raw.is_absolute() else (deck.parent / raw).resolve()
    target.mkdir(parents=True, exist_ok=True)
    cfg["global"]["output_dir"] = str(target)
    os.chdir(deck.parent)
    names = [str(v) for v in cfg.get("variables", {})]
    decomposed = int(cfg["global"].get("ranks", 1)) > 1 or bool(cfg["global"].get("rank_devices"))
    ranks = None
    if decomposed:
        from pinnacle.parallel import run_decomposed

        tic = time.time()
        out = run_decomposed(cfg)
        total = time.time() - tic
        ranks = out["ranks"]
        wall = max(r["wall_s"] for r in ranks)
        n_elements = out["n_elements"]
        result = {"time": out["time"], "history": {"all": out["history"]}, "snapshots": out["snapshots"]}
        result.update(out["fields"])
        arrays = {"nodes": np.asarray(out["nodes"], dtype=float)}
    else:
        runner = Runner(cfg)
        runner.build()
        tic = time.time()
        result = runner.model.train()["composed"]
        wall = total = time.time() - tic
        n_elements = int(runner.model.composed.mesh.n_elements)
        arrays = {"nodes": np.asarray(runner.model.composed.mesh.nodes, dtype=float)}
    for name in names:
        if name in result:
            arrays[name] = np.asarray(result[name], dtype=np.float32)
    frames = result.get("snapshots", []) or []
    if frames:
        arrays["snapshot_times"] = np.array([float(f["time"]) for f in frames])
        for name in names:
            if name in frames[0]:
                arrays[f"snapshot_{name}"] = np.stack(
                    [np.asarray(f[name], dtype=np.float32) for f in frames]
                )
    np.savez_compressed(target / "fields.npz", **arrays)
    history = next(iter(result.get("history", {}).values()), [])
    if history:
        base = ["step", "t", "dt"] + [k for k in history[0] if k.endswith("_max")]
        service_keys = sorted(
            f"{s}.{k}" for s, d in history[0].get("services", {}).items() for k in d
        )
        with (target / "history.csv").open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(base + service_keys)
            for row in history:
                services = row.get("services", {})
                writer.writerow(
                    [row.get(k) for k in base]
                    + [
                        services.get(key.split(".")[0], {}).get(key.split(".", 1)[1])
                        for key in service_keys
                    ]
                )
    summary = {
        "exit": 0,
        "wall_s": round(wall, 1),
        "steps": len(history),
        "elements": n_elements,
        "nodes": int(arrays["nodes"].shape[0]),
        "frames": len(frames),
        "time_end": float(result.get("time", 0.0)),
        "overrides": applied,
        "backend": str(cfg["global"].get("backend")),
        "use_gpu": bool(cfg["global"].get("use_gpu")),
        "devices": cfg["global"].get("devices"),
    }
    if ranks is not None:
        # wall_s: the slowest rank's solve; with the start of the processes and the merge:
        summary["wall_with_start_and_merge_s"] = round(total, 1)
        summary["ranks"] = ranks
    (target / "wall.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
