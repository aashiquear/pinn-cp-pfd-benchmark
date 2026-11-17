#!/usr/bin/env python3
"""Composition and orientation maps of the orientation-field runs of stage C.

Usage: python3 scripts/cryst_orientation_figures.py <tag> <run dir> [<run dir> ...]
Each run directory holds frames.npz with phi, c, theta (phi = 1 liquid). Top row: composition;
bottom row: orientation theta (cyclic color map) in the crystal, liquid black. One GIF of the
orientation map per run.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from scipy import ndimage  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures/crystallization"
plt.rcParams.update({"font.size": 10, "axes.grid": False, "savefig.dpi": 200, "axes.linewidth": 0.8,
                     "mathtext.fontset": "cm"})
CMAP = plt.get_cmap("hsv")


def orientation_image(phi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    rgb = CMAP(theta % 1.0)[..., :3]
    solid = np.clip(1.0 - phi, 0.0, 1.0)[..., None]
    return rgb * solid  # liquid black


def grain_count(phi: np.ndarray, theta: np.ndarray, bins: int = 12) -> int:
    """Connected crystal regions of like orientation (theta binned), larger than 20 cells."""
    solid = phi < 0.5
    total = 0
    binned = np.floor((theta % 1.0) * bins).astype(int)
    for b in range(bins):
        lab, n = ndimage.label(solid & (binned == b))
        if n:
            sizes = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
            total += int(np.sum(sizes > 20))
    return total


def main() -> None:
    tag = sys.argv[1]
    runs = [ROOT / a for a in sys.argv[2:]]
    data = [np.load(r / "frames.npz") for r in runs]
    fig, axes = plt.subplots(2, len(runs), figsize=(3.2 * len(runs), 6.2), squeeze=False)
    for j, (run, d) in enumerate(zip(runs, data)):
        edge = [max((1 - f[0]).max(), (1 - f[:, 0]).max()) for f in d["phi"]]
        ok = [i for i in range(len(edge)) if edge[i] < 0.5]
        i = ok[-1] if ok else len(edge) - 1
        im = axes[0, j].imshow(d["c"][i], origin="lower", cmap="viridis", vmin=float(d["c_s"]) - 0.005,
                               vmax=float(d["c_l"]) + 0.005)
        axes[1, j].imshow(orientation_image(d["phi"][i], d["theta"][i]), origin="lower")
        n_grains = grain_count(d["phi"][i], d["theta"][i])
        axes[0, j].set_title(f"{run.name}\nstep {int(d['steps'][i])}, solid {(1 - d['phi'][i]).mean():.2f}, "
                             f"{n_grains} orientation domains", fontsize=7)
        for a in axes[:, j]:
            a.set_xticks([]); a.set_yticks([])
        print(run.name, "frame", i, "step", int(d["steps"][i]), "orientation domains", n_grains)
    axes[0, 0].set_ylabel("composition $c$")
    axes[1, 0].set_ylabel(r"orientation $\theta$ (crystal)")
    fig.colorbar(im, ax=axes[0, :], shrink=0.8, pad=0.01)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{tag}_maps.png", bbox_inches="tight")
    plt.close(fig)
    for run, d in zip(runs, data):
        fig, ax = plt.subplots(figsize=(3.6, 3.6))
        im = ax.imshow(orientation_image(d["phi"][0], d["theta"][0]), origin="lower")
        ax.set_xticks([]); ax.set_yticks([])
        title = ax.set_title("")

        def draw(i, d=d, im=im, title=title):
            im.set_data(orientation_image(d["phi"][i], d["theta"][i]))
            title.set_text(f"step {int(d['steps'][i])}")
            return im, title

        FuncAnimation(fig, draw, frames=len(d["steps"])).save(FIG / f"{tag}_{run.name}.gif",
                                                              writer=PillowWriter(fps=4), dpi=70)
        plt.close(fig)
    print("wrote", f"figures/crystallization/{tag}_maps.png")


if __name__ == "__main__":
    main()
