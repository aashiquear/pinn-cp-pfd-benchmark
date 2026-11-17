#!/usr/bin/env python3
"""Figures of case B1 (heterogeneous nucleus on a wall; review Fig. 6).

Usage: python3 scripts/cryst_wall_figures.py
Reads results/crystallization/walls/ (fields.vtk of the nucleus decks, wall_table.csv from
scripts/cryst_wall_collect.py) and data/digitized/cryst_review_fig6_barriers.csv; writes
figures/crystallization/b1_wall_nuclei.png and b1_wall_barriers.png.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from vtk_legacy import read_vtk  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.grid": False, "figure.dpi": 150, "savefig.dpi": 200,
                     "axes.linewidth": 0.8, "legend.frameon": False, "mathtext.fontset": "cm"})
RES = ROOT / "results/crystallization/walls"
FIG = ROOT / "figures/crystallization"
DELTA_NM = 2.0 / (2.0 * np.log(9.0))  # nm per scaled length


def field(name: str):
    m = read_vtk(RES / name / "fields.vtk")
    x, y, u = m["nodes"][:, 0], m["nodes"][:, 1], m["point"]["phi"]
    xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
    grid = np.zeros((len(ys), len(xs)))
    ix = np.searchsorted(xs, np.round(x, 9))
    iy = np.searchsorted(ys, np.round(y, 9))
    grid[iy, ix] = u
    # mirror about the symmetry line x = 0
    return (np.concatenate([-xs[:0:-1], xs]) * DELTA_NM, ys * DELTA_NM,
            np.concatenate([grid[:, :0:-1], grid], axis=1))


def main() -> None:
    shown = (20, 40, 80)
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 4.6), sharex=True, sharey=True)
    levels = [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6]
    for row, (label, pattern) in enumerate((("Model I", "b1_wall_model1_dT{:03d}_fem"),
                                            ("Model II", "b1_wall_model2_dT{:03d}_nucleus_fem"))):
        for col, d_t in enumerate(shown):
            ax = axes[row, col]
            x, y, g = field(pattern.format(d_t))
            im = ax.pcolormesh(x, y, g, cmap="Greys", vmin=0, vmax=1, shading="gouraud", rasterized=True)
            ax.contour(x, y, g, levels=levels, colors="black", linewidths=0.6)
            ax.set_xlim(-6, 6); ax.set_ylim(0, 6); ax.set_aspect("equal")
            ax.set_title(rf"{label}, $\Delta T = {d_t}$ K", fontsize=9)
            if row == 1:
                ax.set_xlabel(r"$x$ (nm)")
            if col == 0:
                ax.set_ylabel(r"$y$ (nm)")
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label(r"$\phi$ (0 liquid, 1 crystal)")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "b1_wall_nuclei.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote figures/crystallization/b1_wall_nuclei.png")

    rows = list(csv.DictReader((RES / "wall_table.csv").open()))
    dig = [ln.strip().split(",") for ln in (ROOT / "data/digitized/cryst_review_fig6_barriers.csv").open()
           if ln.startswith("model")]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for model, col, name in (("model_I", "#1f4e79", "Model I"), ("model_II", "#c0392b", "Model II")):
        px = [float(r[1]) for r in dig if r[0] == model]
        py = [float(r[2]) for r in dig if r[0] == model]
        ax.plot(px, py, ".", ms=1.5, color="#999999")
        rr = [r for r in rows if r["model"] == model]
        t = [float(r["undercooling_K"]) for r in rr]
        ax.plot(t, [float(r["ratio_gamma0"]) for r in rr], "o", ms=5, mfc="none", color=col,
                label=rf"{name}, PINNACLE FE, cap with $\gamma(T_m)$")
        ax.plot(t, [float(r["ratio_gammaT"]) for r in rr], "x", ms=4, color=col, alpha=0.7,
                label=rf"{name}, PINNACLE FE, cap with $\gamma(T)$")
    ax.plot([], [], ".", ms=4, color="#999999", label="published (Fig. 6, digitized)")
    ax.set_xlim(0, 100); ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"$T_m - T$ (K)"); ax.set_ylabel(r"$W / W_{\mathrm{cap}}$")
    ax.legend(fontsize=7, loc="lower left")
    fig.savefig(FIG / "b1_wall_barriers.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote figures/crystallization/b1_wall_barriers.png")


if __name__ == "__main__":
    main()
