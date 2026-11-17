#!/usr/bin/env python3
"""Figures and animations of the binary-model runs of stage C (case C2: Fig. 4 a, b, d, e).

Usage: python3 scripts/cryst_binary_figures.py <tag> <run dir> [<run dir> ...]
Each run directory holds frames.npz (phi, c, steps, dx, dt, c_inf, c_l, c_s) and history.csv,
as written by scripts/cryst_binary_reference.py (the PINNACLE runs are converted to the same
layout by scripts/cryst_binary_collect.py). Writes figures/crystallization/<tag>_maps.png,
<tag>_growth.png and one GIF of the composition map per run.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures/crystallization"
plt.rcParams.update({"font.size": 10, "axes.grid": False, "savefig.dpi": 200, "axes.linewidth": 0.8,
                     "legend.frameon": False, "mathtext.fontset": "cm"})
XI_UM = 2.1  # micrometres per unit length


def main() -> None:
    tag = sys.argv[1]
    runs = [ROOT / a for a in sys.argv[2:]]
    data = [np.load(r / "frames.npz") for r in runs]
    fig, axes = plt.subplots(2, len(runs), figsize=(3.1 * len(runs), 6.0), squeeze=False)
    for j, (run, d) in enumerate(zip(runs, data)):
        n = d["phi"].shape[-1]
        size = n * float(d["dx"]) * XI_UM
        # last frame whose crystal has not reached the periodic boundary
        solid = [(1.0 - f).mean() for f in d["phi"]]
        edge = [max((1 - f[0]).max(), (1 - f[:, 0]).max()) for f in d["phi"]]
        ok = [i for i in range(len(solid)) if edge[i] < 0.5]
        i = ok[-1] if ok else len(solid) - 1
        im0 = axes[0, j].imshow(d["c"][i], origin="lower", extent=(0, size, 0, size), cmap="viridis",
                                vmin=float(d["c_s"]) - 0.005, vmax=float(d["c_l"]) + 0.005)
        axes[1, j].imshow(1.0 - d["phi"][i], origin="lower", extent=(0, size, 0, size), cmap="Greys",
                          vmin=0, vmax=1)
        axes[0, j].set_title(f"{run.name}\nstep {int(d['steps'][i])}, solid {solid[i]:.2f}", fontsize=8)
        for a in axes[:, j]:
            a.set_xticks([]); a.set_yticks([])
    axes[0, 0].set_ylabel("composition $c$ (Cu fraction)")
    axes[1, 0].set_ylabel(r"crystal $1-\phi$")
    fig.colorbar(im0, ax=axes[0, :], shrink=0.8, pad=0.01)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{tag}_maps.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.5))
    for run in runs:
        rows = list(csv.DictReader((run / "history.csv").open()))
        t = np.array([float(r["t"]) for r in rows])
        rad = np.array([float(r["equivalent_radius_dx"]) for r in rows])
        keep = rad < 0.4 * data[0]["phi"].shape[-1]
        ax.loglog(t[keep], rad[keep], "o-", ms=3, lw=1.0, label=run.name)
        if keep.sum() > 4:
            slope = np.polyfit(np.log(t[keep][1:]), np.log(rad[keep][1:]), 1)[0]
            print(f"{run.name}: growth exponent d ln R / d ln t = {slope:.3f} ({keep.sum()} points)")
    ax.set_xlabel(r"$t\,D_l/\xi^2$"); ax.set_ylabel(r"equivalent radius ($\Delta x$)")
    ax.legend(fontsize=7)
    fig.savefig(FIG / f"{tag}_growth.png", bbox_inches="tight")
    plt.close(fig)

    for run, d in zip(runs, data):
        fig, ax = plt.subplots(figsize=(3.6, 3.6))
        im = ax.imshow(d["c"][0], origin="lower", cmap="viridis", vmin=float(d["c_s"]) - 0.005,
                       vmax=float(d["c_l"]) + 0.005)
        ax.set_xticks([]); ax.set_yticks([])
        title = ax.set_title("")

        def draw(i, d=d, im=im, title=title):
            im.set_data(d["c"][i])
            title.set_text(f"step {int(d['steps'][i])}")
            return im, title

        FuncAnimation(fig, draw, frames=len(d["steps"])).save(FIG / f"{tag}_{run.name}.gif",
                                                              writer=PillowWriter(fps=4), dpi=70)
        plt.close(fig)
    print("wrote", f"figures/crystallization/{tag}_maps.png, {tag}_growth.png and {len(runs)} GIFs")


if __name__ == "__main__":
    main()
