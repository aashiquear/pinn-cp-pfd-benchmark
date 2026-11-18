#!/usr/bin/env python3
"""Case F2b (eutectic grain of Ag-Cu with an orientation field): measures and figures.

Usage: python3 scripts/cryst_eutectic_grain_figures.py
Reads results/crystallization/eutectic/grain/f2b_grain_*_lam*_fem/ (VTK series); writes
figures/crystallization/f2b_grain.png (composition and orientation of every run at three
frames), f2b_<run>.gif and results/crystallization/eutectic/grain/summary.csv.

Measures at the last frame:
  solid fraction        share of the box with phi < 0.5 (phi = 1 liquid)
  Cu-rich share         share of the solid with c > 0.5 (lever rule at the stated T in the deck)
  reoriented share      share of the solid whose orientation differs from the seed's by more
                        than 0.05 (theta in [0, 1), periodic distance)
  orientation domains   connected solid regions of more than 20 cells in which neighboring
                        cells differ by less than 0.05
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
RES = ROOT / "results/crystallization/eutectic/grain"
FIG = ROOT / "figures/crystallization"


def pdist(a, b):
    d = np.abs(a - b) % 1.0
    return np.minimum(d, 1.0 - d)


def load(run: Path):
    files = series_files(run)
    first = read_vtk(files[0][1])
    n = first["nodes"]
    xs, ys = np.unique(np.round(n[:, 0], 6)), np.unique(np.round(n[:, 1], 6))
    ix, iy = np.searchsorted(xs, np.round(n[:, 0], 6)), np.searchsorted(ys, np.round(n[:, 1], 6))
    frames = []
    for _, path in files:
        point = read_vtk(path)["point"]
        out = []
        for key in ("phi", "c", "theta"):
            a = np.zeros((len(ys), len(xs)))
            a[iy, ix] = point[key if key in point else f"field_{key}"]
            out.append(a)
        frames.append(out)
    return frames


def measures(phi, c, theta, theta_seed):
    solid = phi < 0.5
    smooth = solid.copy()
    for axis in (0, 1):
        sharp = pdist(theta, np.roll(theta, -1, axis=axis)) > 0.05
        smooth &= ~sharp
        smooth &= ~np.roll(sharp, 1, axis=axis)
    labels, count = ndimage.label(smooth)
    sizes = ndimage.sum(np.ones_like(labels), labels, index=np.arange(1, count + 1)) if count else np.array([])
    return {"solid_fraction": float(solid.mean()), "cu_rich_share": float((c[solid] > 0.5).mean()),
            "reoriented_share": float((pdist(theta[solid], theta_seed) > 0.05).mean()),
            "orientation_domains": int(np.sum(sizes > 20))}


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    runs = sorted(p for p in RES.glob("f2b_grain_*_lam*_fem") if (p / "fields.pvd").exists())
    rows = []
    fig, axes = plt.subplots(2 * len(runs), 3, figsize=(8.4, 5.6 * len(runs)), squeeze=False)
    for i, run in enumerate(runs):
        s = float(re.search(r"_S([0-9.]+)_", run.name).group(1))
        chi = float(re.search(r"_chi([0-9.]+)_", run.name).group(1))
        frames = load(run)
        phi0, _, theta0 = frames[0]
        theta_seed = float(np.median(theta0[phi0 < 0.5]))
        rows.append({"run": run.name, "S": s, "chi_l": chi, "frames": len(frames),
                     **measures(*frames[-1], theta_seed)})
        for j, k in enumerate((0, len(frames) // 2, len(frames) - 1)):
            phi, c, theta = frames[k]
            axes[2 * i, j].imshow(c, origin="lower", cmap="RdBu_r", vmin=0, vmax=1)
            axes[2 * i, j].contour(phi, levels=[0.5], colors="k", linewidths=0.5)
            axes[2 * i, j].set_title(f"$S$ = {s:g}, $\\chi_l$ = {chi:g}, frame {k} of {len(frames) - 1}", fontsize=8)
            axes[2 * i + 1, j].imshow(np.ma.masked_where(phi > 0.5, theta), origin="lower", cmap="hsv", vmin=0, vmax=1)
            axes[2 * i + 1, j].set_facecolor("0.85")
        gfig, gax = plt.subplots(1, 2, figsize=(6.0, 3.1))
        im1 = gax[0].imshow(frames[0][1], origin="lower", cmap="RdBu_r", vmin=0, vmax=1)
        im2 = gax[1].imshow(np.ma.masked_where(frames[0][0] > 0.5, frames[0][2]), origin="lower", cmap="hsv", vmin=0, vmax=1)
        for a in gax:
            a.set_xticks([]); a.set_yticks([]); a.set_facecolor("0.85")

        def draw(k, frames=frames, im1=im1, im2=im2):
            im1.set_data(frames[k][1])
            im2.set_data(np.ma.masked_where(frames[k][0] > 0.5, frames[k][2]))
            return im1, im2

        FuncAnimation(gfig, draw, frames=len(frames)).save(FIG / f"f2b_{run.name}.gif", writer=PillowWriter(fps=3), dpi=80)
        plt.close(gfig)
    for a in axes.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    fig.savefig(FIG / "f2b_grain.png", bbox_inches="tight")
    plt.close(fig)
    with (RES / "summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
    print("wrote figures/crystallization/f2b_grain.png")


if __name__ == "__main__":
    main()
