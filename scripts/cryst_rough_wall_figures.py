#!/usr/bin/env python3
"""Case E5 (nucleation on a rough wall): where the crystals form, figures.

Usage: python3 scripts/cryst_rough_wall_figures.py
Reads results/crystallization/walls_rough/e5_rough_wall_std*_fem/ (VTK series on the mapped
grid of scripts/cryst_rough_wall_mesh.py); writes figures/crystallization/e5_rough_wall.png
(phase field at four early frames per run), e5_rough_wall_<run>.gif and
results/crystallization/walls_rough/summary.csv.

Measures at the first saved frame with a crystalline fraction above 1 %:
  crystals            connected regions with phi < 0.5 (on the logical grid)
  on the wall         crystals with a cell in the first two rows above the wall
  in a groove         crystals on the wall whose lowest point lies in the lower half of the
                      sawtooth (groove bottoms); "on a ridge": in the upper half
  expected on wall    share of the crystals that would touch the first two rows if the nuclei
                      were spread uniformly over the liquid (2 rows of ny): the reference for
                      "no preference"
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from scipy import ndimage  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk, series_files  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/walls_rough"
FIG = ROOT / "figures/crystallization"
DX = 6.25e-3


def logical(nodes: np.ndarray):
    """Index arrays (column, row) of the mapped grid: columns are vertical."""
    xs = np.unique(np.round(nodes[:, 0], 9))
    col = np.searchsorted(xs, np.round(nodes[:, 0], 9))
    row = np.zeros(len(nodes), dtype=int)
    for k in range(len(xs)):
        idx = np.where(col == k)[0]
        row[idx[np.argsort(nodes[idx, 1])]] = np.arange(len(idx))
    return col, row, len(xs), int(row.max()) + 1


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    runs = sorted(p for p in RES.glob("e5_rough_wall_std*_fem") if (p / "fields.pvd").exists())
    rows = []
    fig, axes = plt.subplots(len(runs), 4, figsize=(13.0, 1.9 * len(runs) + 0.4), squeeze=False)
    for i, run in enumerate(runs):
        std = float(re.search(r"_std([0-9.]+)_", run.name).group(1))
        files = series_files(run)
        first = read_vtk(files[0][1])
        nodes = first["nodes"]
        col, row, nc, nr = logical(nodes)
        X, Y = np.zeros((nr, nc)), np.zeros((nr, nc))
        X[row, col], Y[row, col] = nodes[:, 0] / DX, nodes[:, 1] / DX
        wall = Y[0]
        frames = []
        for _, path in files:
            point = read_vtk(path)["point"]
            a = np.zeros((nr, nc))
            a[row, col] = point["field_phi" if "field_phi" in point else "phi"]
            frames.append(a)
        fractions = [float((f < 0.5).mean()) for f in frames]
        k0 = next((k for k, f in enumerate(fractions) if f > 0.01), len(frames) - 1)
        solid = frames[k0] < 0.5
        labels, count = ndimage.label(solid)
        on_wall = groove = ridge = 0
        mid = 0.5 * (wall.min() + wall.max())
        for lab in range(1, count + 1):
            rr, cc = np.where(labels == lab)
            if rr.min() <= 1:
                on_wall += 1
                lowest = cc[np.argmin(rr)]
                if wall[lowest] < mid:
                    groove += 1
                else:
                    ridge += 1
        steps = np.linspace(0, len(list(csv.DictReader((run / "solver_history.csv").open()))), len(files)).round()
        rows.append({"noise_std": std, "frame_step": int(steps[k0]), "crystalline_fraction": fractions[k0],
                     "crystals": count, "on_wall": on_wall, "in_groove": groove, "on_ridge": ridge,
                     "share_on_wall": on_wall / max(count, 1), "expected_on_wall": 2.0 / nr})
        picks = [k for k in (k0 - 1, k0, k0 + 1, k0 + 2) if 0 <= k < len(frames)]
        for ax, k in zip(axes[i], picks):
            ax.pcolormesh(X, Y, frames[k], cmap="gray", vmin=0, vmax=1, shading="gouraud", rasterized=True)
            ax.fill_between(X[0], wall, wall.min() - 4, color="#c8553d", lw=0)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"noise {std:g}, step {int(steps[k])}, X = {fractions[k]:.2f}", fontsize=8)
        gfig, gax = plt.subplots(figsize=(5.2, 2.8))
        mesh = gax.pcolormesh(X, Y, frames[0], cmap="gray", vmin=0, vmax=1, shading="gouraud")
        gax.fill_between(X[0], wall, wall.min() - 4, color="#c8553d", lw=0)
        gax.set_aspect("equal"); gax.set_xticks([]); gax.set_yticks([])
        last = next((k for k, f in enumerate(fractions) if f > 0.98), len(frames) - 1)
        FuncAnimation(gfig, lambda k, mesh=mesh, frames=frames: (mesh.set_array(frames[k].ravel()),),
                      frames=last + 1).save(FIG / f"e5_rough_wall_{run.name}.gif", writer=PillowWriter(fps=3), dpi=80)
        plt.close(gfig)
    fig.tight_layout()
    fig.savefig(FIG / "e5_rough_wall.png", bbox_inches="tight")
    plt.close(fig)
    with (RES / "summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
    print("wrote figures/crystallization/e5_rough_wall.png")


if __name__ == "__main__":
    main()
