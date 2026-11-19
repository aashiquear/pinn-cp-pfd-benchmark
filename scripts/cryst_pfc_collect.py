#!/usr/bin/env python3
"""Part IV, stage D: collect the VTK series of a phase-field crystal run into frames.npz.

Usage: python3 scripts/cryst_pfc_collect.py [RUN ...]     (default: every run folder)
Reads results/crystallization/pfc/<run>/fields_step*.vtk, writes <run>/frames.npz (psi on the
structured node grid, x, y). NumPy only, so it runs on the compute host; the figure script
reads the cache.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk, series_files  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/pfc"


def load(run: str) -> dict:
    """Frames of a run on its structured node grid: psi (nt, ny, nx), t, x, y, history."""
    folder = RES / run
    cache = folder / "frames.npz"
    if not cache.exists() or cache.stat().st_mtime < (folder / "fields.pvd").stat().st_mtime:
        files = series_files(folder)
        first = read_vtk(files[0][1])
        nodes = first["nodes"]
        xs, ys = np.unique(np.round(nodes[:, 0], 8)), np.unique(np.round(nodes[:, 1], 8))
        ix, iy = np.searchsorted(xs, np.round(nodes[:, 0], 8)), np.searchsorted(ys, np.round(nodes[:, 1], 8))
        frames = []
        for _, path in files:
            point = read_vtk(path)["point"]
            grid = np.zeros((len(ys), len(xs)))
            grid[iy, ix] = point["field_psi" if "field_psi" in point else "psi"]
            frames.append(grid)
        np.savez_compressed(cache, psi=np.array(frames, dtype=np.float32), x=xs, y=ys)
    d = dict(np.load(cache))
    hist = list(csv.DictReader((folder / "solver_history.csv").open()))
    d["t_hist"] = np.array([float(r["t"]) for r in hist])
    d["energy"] = np.array([float(r["energy"]) for r in hist])
    d["newton"] = np.array([float(r["newton_iters"]) for r in hist])
    nt = len(d["psi"])
    # snapshots: step 0, every k-th step and the last one
    every = int(re.search(r"snapshots_every = (\d+)", (ROOT / f"decks/crystallization/pfc/{run}.i").read_text()).group(1))
    steps = [0] + list(range(every, len(hist) + 1, every))
    if steps[-1] != len(hist):
        steps.append(len(hist))
    d["t"] = np.array([0.0] + [d["t_hist"][s - 1] for s in steps[1:]])[:nt]
    d["area"] = float((d["x"][-1] - d["x"][0]) * (d["y"][-1] - d["y"][0]))
    return d



if __name__ == "__main__":
    names = sys.argv[1:] or sorted(p.name for p in RES.iterdir() if (p / "fields.pvd").exists())
    for name in names:
        d = load(name)
        print(name, d["psi"].shape, flush=True)
