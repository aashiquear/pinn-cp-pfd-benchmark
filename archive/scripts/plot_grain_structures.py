#!/usr/bin/env python3
"""Grain maps of the Voronoi seeding and the relaxed structures at each horizon.

One panel per structure, coloured by grain id with a fixed categorical cycle,
boundaries drawn in dark ink; the title of each panel carries the number
that settles 'natural': grain count and the triple-junction angle spread.

Usage: plot_grain_structures.py --labels name=label_dir [...] --seeds seeds.npz --out fig.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree


def voronoi_grid(seeds, xs, ys):
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    _, idx = cKDTree(seeds).query(np.column_stack([X.ravel(), Y.ravel()]))
    return idx.reshape(X.shape)


def panel(ax, gid, xs, ys, title):
    rng = np.random.default_rng(3)
    n = int(gid.max()) + 1
    lut = rng.permutation(n)
    cmap = plt.get_cmap("tab20", 20)
    img = cmap(lut[gid] % 20)
    edge = np.zeros_like(gid, dtype=bool)
    edge[1:, :] |= gid[1:, :] != gid[:-1, :]
    edge[:, 1:] |= gid[:, 1:] != gid[:, :-1]
    img[edge] = (0.12, 0.12, 0.12, 1.0)
    ax.imshow(np.transpose(img, (1, 0, 2)), origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], interpolation="nearest")
    ax.set_title(title, fontsize=9, loc="left"); ax.set_xticks([]); ax.set_yticks([])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="+", required=True, help="name=dir with grains.npz and label.json")
    ap.add_argument("--seeds", required=True, help="npz with 'seeds' (n, 2)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cases = [s.split("=", 1) for s in a.labels]
    fig, axes = plt.subplots(1, len(cases) + 1, figsize=(4.0 * (len(cases) + 1), 2.6), dpi=150)
    z0 = np.load(Path(cases[0][1]) / "grains.npz")
    xs, ys = z0["xs"], z0["ys"]
    seeds = np.load(a.seeds)["seeds"]
    panel(axes[0], voronoi_grid(seeds, xs, ys), xs, ys, f"Voronoi seeding, {len(seeds)} grains")
    for ax, (name, d) in zip(axes[1:], cases):
        z = np.load(Path(d) / "grains.npz"); rec = json.loads((Path(d) / "label.json").read_text())
        m = rec["morphology"]
        spread = m.get("junction_rms_dev_from_120") if isinstance(m, dict) else None
        sp = f", junction rms dev {spread:.1f} deg" if isinstance(spread, (int, float)) else ""
        panel(ax, z["gid"], xs, ys, f"{name}: {rec['n_grains']} grains{sp}")
    fig.tight_layout(); Path(a.out).parent.mkdir(parents=True, exist_ok=True); fig.savefig(a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
