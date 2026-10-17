#!/usr/bin/env python3
"""Plan-view field maps of Ni dog-bone runs, driven by one deck.

Reads ``decks/ni/figures/fields_*.i``: for every run in ``[Runs]`` and every time
in ``[Figures] times`` it picks the snapshot of the run's ``fields.pvd`` nearest
that time, reads the legacy-ASCII VTK (points, hex8 cells, point displacement,
the cell scalar ``field``), and draws the deformed bottom face of every element
in the x-y plane coloured by the scalar. One row per run, one column per time,
one colour scale per column; each panel names the run's gauge strain at that
time (from its results.csv).

Usage: ni_fields.py decks/ni/figures/fields_mesh.i
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_figures import _resolve, _str  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.labelsize": 9, "mathtext.fontset": "dejavusans"})


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    return Translator().translate(InputParser().parse(path.read_text()))


def snapshots(pvd: Path) -> list[tuple[float, Path]]:
    text = pvd.read_text()
    return [(float(t), pvd.parent / f) for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', text)]


def read_vtk(path: Path, field: str):
    lines = path.read_text().split("\n")
    i = 0
    pts = cells = disp = vals = None
    while i < len(lines):
        tok = lines[i].split()
        if tok[:1] == ["POINTS"]:
            n = int(tok[1])
            data = []
            i += 1
            while len(data) < 3 * n:
                data += lines[i].split()
                i += 1
            pts = np.array(data, dtype=float).reshape(n, 3)
            continue
        if tok[:1] == ["CELLS"]:
            n = int(tok[1])
            cells = [list(map(int, lines[i + 1 + k].split()[1:])) for k in range(n)]
            i += n + 1
            continue
        if tok[:2] == ["VECTORS", "displacement"]:
            n = pts.shape[0]
            data = []
            i += 1
            while len(data) < 3 * n:
                data += lines[i].split()
                i += 1
            disp = np.array(data, dtype=float).reshape(n, 3)
            continue
        if tok[:2] == ["SCALARS", field]:
            n = len(cells)
            data = []
            i += 2  # skip LOOKUP_TABLE
            while len(data) < n:
                data += lines[i].split()
                i += 1
            vals = np.array(data, dtype=float)
            continue
        i += 1
    if vals is None:
        raise SystemExit(f"{path}: no cell scalar {field!r}")
    return pts, cells, disp, vals


def gauge_at(results: Path, t: float) -> float:
    rows = list(csv.DictReader(open(results)))
    tt = np.array([float(r["time"]) for r in rows])
    ee = np.array([float(r["gauge_strain"]) for r in rows])
    return float(np.interp(t, tt, ee))


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    f, runs = cfg["figures"], cfg["runs"]
    field = _str(f["field"])
    times = [float(x) for x in _str(f["times"]).split()]
    x0, x1 = (float(x) for x in _str(f["x_range"]).split())
    maps = {}
    for name, r in runs.items():
        snaps = snapshots(_resolve(_str(r["pvd"]), deck))
        for t in times:
            ts, path = min(snaps, key=lambda s: abs(s[0] - t))
            pts, cells, disp, vals = read_vtk(path, field)
            xy = (pts + disp)[:, :2]
            # elements whose undeformed centroid lies in the x window
            keep = [k for k, c in enumerate(cells) if x0 <= pts[c, 0].mean() <= x1]
            polys = [xy[cells[k][:4]] for k in keep]
            maps[name, t] = (polys, vals[keep], ts, gauge_at(_resolve(_str(r["results"]), deck), ts))
    names = list(runs)
    fig, axes = plt.subplots(len(names), len(times), figsize=(4.2 * len(times), 2.0 * len(names) + 0.8),
                             dpi=int(f["dpi"]), squeeze=False)
    for j, t in enumerate(times):
        vmax = max(float(maps[n, t][1].max()) for n in names)
        for i, n in enumerate(names):
            ax = axes[i, j]
            polys, vals, ts, eg = maps[n, t]
            pc = PolyCollection(polys, array=vals, cmap="viridis", edgecolors="none")
            pc.set_clim(0.0, vmax)
            ax.add_collection(pc)
            allxy = np.concatenate(polys)
            ax.set_xlim(allxy[:, 0].min() - 0.2, allxy[:, 0].max() + 0.2)
            ax.set_ylim(-1.6, 1.6)
            ax.set_aspect("equal")
            ax.set_title("%s\nbar %.0f%%, gauge %.1f%%, max %.2f"
                         % (_str(runs[n]["label"]), 5.0 * ts, 100 * eg, vals.max()), fontsize=7)
            ax.tick_params(labelsize=6)
            if j == 0:
                ax.set_ylabel("y (mm)")
            if i == len(names) - 1:
                ax.set_xlabel("x (mm)")
        cb = fig.colorbar(pc, ax=axes[:, j], orientation="horizontal", fraction=0.06, pad=0.18)
        cb.set_label(_str(f["label"]), fontsize=7)
        cb.ax.tick_params(labelsize=6)
    out = _resolve(_str(f["out"]), deck)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
