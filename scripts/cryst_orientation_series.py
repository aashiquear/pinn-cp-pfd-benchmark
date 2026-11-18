#!/usr/bin/env python3
"""Parameter series of the orientation-field model: maps, measures and supporting plots.

Usage:
  python3 scripts/cryst_orientation_series.py <tag> <symmetry k> <label> <value>=<run dir> [...]

  tag       prefix of the output files (figures/crystallization/<tag>_*.png,
            results/crystallization/orientation/<tag>_measures.csv)
  k         symmetry index of the crystal (theta is the crystal angle in units of 2 pi / k)
  label     axis label of the varied parameter (mathtext allowed, quote it)
  value=dir parameter value and the run directory (holding frames.npz from
            cryst_binary_collect.py or cryst_orientation_reference.py), relative to the repo

Measures, at the last frame in which the crystal has not reached the box edge:
  equivalent radius            sqrt(solid area / pi), in cells
  maximum radius               largest distance of a crystal cell from the seed, in cells
  compactness                  solid area / (pi maximum radius^2): 1 for a disc, small for needles
  growth exponent              slope of log(equivalent radius) against log(time), second half of the run
  orientation domains          connected crystal regions of like orientation (12 bins, > 20 cells)
  mean domain area             solid area / orientation domains, in cells
  mean misorientation          mean |delta theta| between neighboring crystal cells (0.25: random)
  boundary fraction            fraction of crystal neighbor pairs with |delta theta| > 0.05
  radial alignment             <cos(k (alpha - beta))>: alpha the fast growth direction of the cell
                               (2 pi theta / k), beta the direction from the seed to the cell;
                               1: fast directions point outward everywhere, 0: no relation
  trapped solute               mean (c - c_S) / (c_L - c_S) in the crystal: 0 full partitioning at
                               the solidus, (1 - S) complete trapping of the liquid composition
Plots: <tag>_maps.png (composition and orientation), <tag>_measures.png (measures against the
parameter), <tag>_misorientation.png (distribution of neighbor misorientation in the crystal),
<tag>_growth.png (equivalent radius against time).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cryst_orientation_figures import FIG, ROOT, grain_count, orientation_image  # noqa: E402

COLORS = ["#1b3a5c", "#2a78a8", "#3aa39b", "#8fb84e", "#e0a030", "#c8553d"]


def wrapped(d: np.ndarray) -> np.ndarray:
    return np.abs(d - np.round(d))


def last_free_frame(phi: np.ndarray) -> int:
    edge = [max((1 - f[0]).max(), (1 - f[:, 0]).max()) for f in phi]
    ok = [i for i in range(len(edge)) if edge[i] < 0.5]
    return ok[-1] if ok else len(edge) - 1


def measures(d: np.lib.npyio.NpzFile, k: int) -> tuple[dict[str, float], np.ndarray]:
    i = last_free_frame(d["phi"])
    phi, theta, c = d["phi"][i].astype(float), d["theta"][i].astype(float), d["c"][i].astype(float)
    n = phi.shape[0]
    solid = phi < 0.5
    area = float(solid.sum())
    yy, xx = np.mgrid[0:n, 0:n]
    cx = cy = n / 2 - 0.5
    r = np.hypot(xx - cx, yy - cy)
    r_max = float(r[solid].max()) if area else 0.0
    # neighbor misorientation inside the crystal (both cells solid)
    pairs = []
    for axis in (0, 1):
        dth = wrapped(np.diff(theta, axis=axis))
        both = solid[1:, :] & solid[:-1, :] if axis == 0 else solid[:, 1:] & solid[:, :-1]
        pairs.append(dth[both])
    mis = np.concatenate(pairs) if pairs else np.zeros(0)
    # radial alignment, outside the seed (r > 20 cells)
    beta = np.arctan2(yy - cy, xx - cx)
    alpha = 2.0 * np.pi * theta / k
    ring = solid & (r > 20)
    align = float(np.cos(k * (alpha - beta))[ring].mean()) if ring.any() else float("nan")
    # growth exponent from the second half of the frames before the box edge
    steps = d["steps"][: i + 1].astype(float)
    rad = np.array([np.sqrt((f < 0.5).sum() / np.pi) for f in d["phi"][: i + 1]])
    half = slice(max(1, (i + 1) // 2), i + 1)
    expo = float(np.polyfit(np.log(steps[half]), np.log(rad[half]), 1)[0]) if i >= 4 else float("nan")
    domains = grain_count(phi, theta)
    c_l, c_s = float(d["c_l"]), float(d["c_s"])
    out = {
        "step": float(d["steps"][i]),
        "solid_fraction": area / n**2,
        "equivalent_radius": float(np.sqrt(area / np.pi)),
        "maximum_radius": r_max,
        "compactness": area / (np.pi * r_max**2) if r_max else float("nan"),
        "growth_exponent": expo,
        "orientation_domains": float(domains),
        "mean_domain_area": area / max(domains, 1),
        "mean_misorientation": float(mis.mean()) if mis.size else float("nan"),
        "boundary_fraction": float((mis > 0.05).mean()) if mis.size else float("nan"),
        "radial_alignment": align,
        "trapped_solute": float((c[solid].mean() - c_s) / (c_l - c_s)) if area else float("nan"),
    }
    return out, mis


def main() -> None:
    tag, k, label = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    items = [a.split("=", 1) for a in sys.argv[4:]]
    values = [v for v, _ in items]
    data = [np.load(ROOT / p / "frames.npz") for _, p in items]
    rows, hists = [], []
    for v, d in zip(values, data):
        m, mis = measures(d, k)
        rows.append({"value": v, **m})
        hists.append(mis)
    out = ROOT / "results/crystallization/orientation"
    out.mkdir(parents=True, exist_ok=True)
    with (out / f"{tag}_measures.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{key}={val:.4g}" if isinstance(val, float) else f"{key}={val}" for key, val in r.items()))
    plt.rcParams.update({"font.size": 9, "axes.grid": False, "savefig.dpi": 200, "axes.linewidth": 0.8,
                         "mathtext.fontset": "cm"})
    # ---- maps
    fig, axes = plt.subplots(2, len(data), figsize=(2.6 * len(data), 5.4), squeeze=False)
    for j, (v, d) in enumerate(zip(values, data)):
        i = last_free_frame(d["phi"])
        im = axes[0, j].imshow(d["c"][i], origin="lower", cmap="viridis", vmin=float(d["c_s"]) - 0.005,
                               vmax=float(d["c_l"]) + 0.005)
        axes[1, j].imshow(orientation_image(d["phi"][i], d["theta"][i]), origin="lower")
        axes[0, j].set_title(f"{label} = {v}\nstep {int(d['steps'][i])}", fontsize=9)
        for a in axes[:, j]:
            a.set_xticks([]); a.set_yticks([])
    axes[0, 0].set_ylabel("composition $c$")
    axes[1, 0].set_ylabel(r"orientation $\theta$")
    fig.colorbar(im, ax=axes[0, :], shrink=0.8, pad=0.01)
    fig.savefig(FIG / f"{tag}_maps.png", bbox_inches="tight")
    plt.close(fig)
    # ---- measures against the parameter
    try:
        x = np.array([float(v) for v in values])
        names = None
    except ValueError:  # named cases instead of a numeric parameter
        x, names = np.arange(len(values), dtype=float), values
    panels = [("equivalent_radius", "equivalent radius (cells)"), ("compactness", "compactness"),
              ("growth_exponent", "growth exponent"), ("orientation_domains", "orientation domains"),
              ("mean_misorientation", r"mean misorientation ($2\pi/k$)"), ("radial_alignment", "radial alignment"),
              ("trapped_solute", "trapped solute"), ("boundary_fraction", "boundary fraction")]
    fig, axes = plt.subplots(2, 4, figsize=(11.5, 5.0))
    for ax, (key, name) in zip(axes.ravel(), panels):
        ax.plot(x, [r[key] for r in rows], "o-", color=COLORS[1], lw=1.2, ms=4)
        ax.set_xlabel(label); ax.set_ylabel(name)
        if names:
            ax.set_xticks(x); ax.set_xticklabels(names, fontsize=7)
        if key == "orientation_domains":
            ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(FIG / f"{tag}_measures.png", bbox_inches="tight")
    plt.close(fig)
    # ---- misorientation distribution
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    bins = np.linspace(0.0, 0.5, 51)
    for j, (v, mis) in enumerate(zip(values, hists)):
        if mis.size:
            h, _ = np.histogram(mis, bins=bins, density=True)
            ax.plot(0.5 * (bins[1:] + bins[:-1]), np.maximum(h, 1e-4), color=COLORS[j % len(COLORS)], lw=1.2,
                    label=f"{label} = {v}")
    ax.set_yscale("log")
    ax.set_xlabel(r"neighbor misorientation $|\Delta\theta|$ ($2\pi/k$)"); ax.set_ylabel("probability density")
    ax.legend(frameon=False, fontsize=7)
    fig.savefig(FIG / f"{tag}_misorientation.png", bbox_inches="tight")
    plt.close(fig)
    # ---- growth
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    for j, (v, d) in enumerate(zip(values, data)):
        i = last_free_frame(d["phi"])
        rad = [np.sqrt((f < 0.5).sum() / np.pi) for f in d["phi"][1 : i + 1]]
        ax.loglog(d["steps"][1 : i + 1], rad, color=COLORS[j % len(COLORS)], lw=1.2, label=f"{label} = {v}")
    ax.set_xlabel("step"); ax.set_ylabel("equivalent radius (cells)")
    ax.legend(frameon=False, fontsize=7)
    fig.savefig(FIG / f"{tag}_growth.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote", f"figures/crystallization/{tag}_{{maps,measures,misorientation,growth}}.png")


if __name__ == "__main__":
    main()
