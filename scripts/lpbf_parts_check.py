#!/usr/bin/env python3
"""Part V: run one deck as a single process and cut into parts, and compare the two.

usage: scripts/lpbf_parts_check.py <deck.i> <devices or number of ranks> [key=value ...]

`devices` is a quoted list such as "cuda:0 cuda:1 cuda:2 cuda:3" (one process per entry,
PINNACLE ADR 0057), or a number of ranks on the device of the deck. `key=value` replaces a
deck key by its dotted path. Options before the deck: --transport process|thread,
--partition slabs|rcb, --no-single (time the parts only).

Prints one JSON document: wall time of the single run and of the parts, per field the
shape, the largest value and the largest difference between the two runs (shapes are
checked before anything is subtracted), the last history row of both, and per rank its
device, sizes, exchanges, linear iterations, device memory and wall time.
Works for the thermal, flow and grain decks (`pf_fem_explicit`) and for the stress decks
(`cp_quasi_static`).
"""

from __future__ import annotations

import copy
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


def main() -> None:
    args = sys.argv[1:]
    transport, partition, single_too = "process", "slabs", True
    while args and args[0].startswith("--"):
        flag = args.pop(0)
        if flag == "--transport":
            transport = args.pop(0)
        elif flag == "--partition":
            partition = args.pop(0)
        elif flag == "--no-single":
            single_too = False
        else:
            raise SystemExit(__doc__)
    if len(args) < 2:
        raise SystemExit(__doc__)
    from pinnacle.io.loader import load_input
    from pinnacle.parallel import run_decomposed
    from pinnacle.plugins import discover_components
    from pinnacle.runner import Runner

    deck = Path(args[0]).resolve()
    devices = args[1].split() if not args[1].isdigit() else None
    ranks = int(args[1]) if args[1].isdigit() else len(devices)
    discover_components(force=True)
    cfg = load_input(str(deck))
    for item in args[2:]:
        dotted, value = item.split("=", 1)
        node = cfg
        keys = dotted.split(".")
        for key in keys[:-1]:
            node = node[key]
        try:
            node[keys[-1]] = json.loads(value)
        except json.JSONDecodeError:
            node[keys[-1]] = value
    os.chdir(deck.parent)
    report: dict = {"deck": deck.name, "ranks": ranks, "devices": devices, "transport": transport}
    single = None
    if single_too:
        runner = Runner(copy.deepcopy(cfg))
        runner.build()
        tic = time.time()
        single = runner.model.train()["composed"]
        report["single_wall_s"] = round(time.time() - tic, 2)
        report["elements"] = int(runner.model.composed.mesh.n_elements)
        del runner
    tic = time.time()
    out = run_decomposed(copy.deepcopy(cfg), ranks=ranks, devices=devices, partition=partition, transport=transport)
    report["parts_wall_with_start_s"] = round(time.time() - tic, 2)
    report["parts_wall_s"] = max(r["wall_s"] for r in out["ranks"])
    report["elements"] = out["n_elements"]
    if single is not None:
        fields = {}
        for name, values in out["fields"].items():
            if values is None or name not in single:
                continue
            reference = np.asarray(single[name], dtype=float)
            row = {"shape_single": list(reference.shape), "shape_parts": list(values.shape)}
            if reference.shape == values.shape:
                row["largest_value"] = float(np.abs(reference).max())
                row["largest_difference"] = float(np.abs(reference - values).max())
            fields[name] = row
        report["fields"] = fields
        history = single["history"]
        history = history if isinstance(history, list) else next(iter(history.values()))

        def numbers(row: dict) -> dict:
            flat = {k: v for k, v in row.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}
            for service, values in (row.get("services") or {}).items():
                flat.update({f"{service}.{k}": v for k, v in values.items() if isinstance(v, (int, float))})
            return flat

        report["steps"] = [len(history), len(out["history"])]
        report["last_row_single"] = numbers(history[-1])
        report["last_row_parts"] = numbers(out["history"][-1])
    report["rank_stats"] = [
        {k: r.get(k) for k in ("rank", "device", "device_index", "elements", "shared_nodes", "exchanges",
                               "linear_solves", "linear_iterations", "device_peak_bytes", "build_s", "wall_s")}
        for r in out["ranks"]
    ]
    print("PARTS-CHECK " + json.dumps(report))


if __name__ == "__main__":
    main()
