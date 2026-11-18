#!/usr/bin/env python3
"""Case C6, heat treatments (Fig. 9(d), (e) of granasy 2005 spherulites): figures and measures.

Usage: python3 scripts/cryst_heat_treatment_figures.py
Reads results/crystallization/orientation/heat_treatment/c6_fig9{d_direct,e_short,e_mild,e_deep}_fem/
(frames.npz from cryst_binary_collect.py, solver_history.csv); writes
figures/crystallization/c6_heat_treatment.png, c6_heat_treatment_<run>.gif and
results/crystallization/orientation/heat_treatment/summary.csv.

Measures (time counted from the end of the hold at 1595 K, step 3000):
  crystals        largest number of connected regions with phi < 0.5 over the saved frames with
                  a crystalline fraction below 0.6 (periodic images are not merged)
  t_10, t_50      time at which the crystalline fraction passes 0.1 and 0.5
  domains         orientation domains of the final frame (grain_count of
                  cryst_orientation_figures.py), and the mean domain area in cells
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
from scipy import ndimage  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cryst_orientation_figures import FIG, ROOT, grain_count, orientation_image  # noqa: E402

RES = ROOT / "results/crystallization/orientation/heat_treatment"
RUNS = {"c6_fig9d_direct_fem": "direct quench: 1595 K, then 1574 K",
        "c6_fig9e_short_fem": "20 steps at 1550 K: 1595 K, 1550 K, then 1574 K",
        "c6_fig9e_mild_fem": "100 steps at 1550 K: 1595 K, 1550 K, then 1574 K",
        "c6_fig9e_deep_fem": "100 steps at 1520 K: 1595 K, 1520 K, then 1574 K"}
RUNS = {k: v for k, v in RUNS.items() if (RES / k / "frames.npz").exists()}
T_HOLD = 3000 * 4.75e-6


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    fig = plt.figure(figsize=(12.5, 3.2 * len(RUNS)))
    gs = fig.add_gridspec(len(RUNS), 5, width_ratios=[1, 1, 1, 1, 1.25], wspace=0.08, hspace=0.18)
    ax_x = fig.add_subplot(gs[0, 4])
    ax_n = fig.add_subplot(gs[len(RUNS) - 1, 4])
    colors = ["#1b3a5c", "#8a6fb0", "#3aa39b", "#c8553d"]
    rows = []
    for i, (name, label) in enumerate(RUNS.items()):
        run = RES / name
        d = np.load(run / "frames.npz")
        phi, theta, steps = d["phi"], d["theta"], d["steps"]
        frac = (phi < 0.5).mean(axis=(1, 2))
        count = np.array([ndimage.label(p < 0.5)[1] for p in phi])
        hist = list(csv.DictReader((run / "solver_history.csv").open()))
        t = np.array([float(h["t"]) for h in hist]) - T_HOLD
        x = 1.0 - np.array([float(h["phi_mean"]) for h in hist])
        x = np.clip((x - np.median(x[: len(x) // 2])) / (x[-1] - np.median(x[: len(x) // 2])), 0, 1)
        t10, t50 = (float(t[np.argmax(x > v)]) for v in (0.1, 0.5))
        open_frames = (frac > 0.005) & (frac < 0.6)
        domains = grain_count(phi[-1], theta[-1])
        rows.append({"run": name, "crystals": int(count[open_frames].max()) if open_frames.any() else 0,
                     "t_10": t10, "t_50": t50, "final_fraction": float(frac[-1]), "domains": domains,
                     "mean_domain_area_cells": phi[-1].size / max(domains, 1)})
        first = int(np.argmax(frac > 0.01))  # first saved frame with crystal
        picks = [max(first - 1, 0)] + [int(np.argmin(np.abs(frac - v))) for v in (0.3, 0.7)] + [len(steps) - 1]
        for j, k in enumerate(picks):
            ax = fig.add_subplot(gs[i, j])
            ax.set_xticks([]); ax.set_yticks([])
            if j in (1, 2) and not 0.05 < frac[k] < 0.95:  # no saved frame during the transformation
                ax.set_frame_on(False)
                continue
            ax.imshow(orientation_image(phi[k], theta[k]), origin="lower")
            ax.set_title(f"step {int(steps[k])}, X = {frac[k]:.2f}", fontsize=8)
            if j == 0:
                ax.set_ylabel(label.split(":")[0], fontsize=9)
        ax_x.plot(t * 1e3, x, color=colors[i], lw=1.3, label=label.split(":")[0])
        ax_n.plot((steps * float(d["dt"]) - T_HOLD)[open_frames] * 1e3, count[open_frames], "o-", color=colors[i],
                  ms=3, lw=1.1, label=label.split(":")[0])
        frames = [orientation_image(phi[k], theta[k]) for k in range(len(steps))]
        gfig, gax = plt.subplots(figsize=(3.2, 3.2))
        gax.set_xticks([]); gax.set_yticks([])
        im = gax.imshow(frames[0], origin="lower")
        FuncAnimation(gfig, lambda k, im=im, frames=frames: (im.set_data(frames[k]),), frames=len(frames)).save(
            FIG / f"c6_heat_treatment_{name}.gif", writer=PillowWriter(fps=4), dpi=80)
        plt.close(gfig)
    ax_x.set_xlim(-0.5, None)
    ax_x.set_ylabel("crystalline fraction $X$")
    ax_x.legend(frameon=False, fontsize=7)
    ax_n.set_xlim(ax_x.get_xlim())
    ax_n.set_xlabel("time after the hold at 1595 K ($10^{-3}$)")
    ax_n.set_ylabel("separate crystals ($X$ < 0.6)")
    for a in (ax_x, ax_n):
        a.yaxis.tick_right(); a.yaxis.set_label_position("right")
    fig.savefig(FIG / "c6_heat_treatment.png", bbox_inches="tight")
    plt.close(fig)
    with (RES / "summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
    print("wrote figures/crystallization/c6_heat_treatment.png")


if __name__ == "__main__":
    main()
