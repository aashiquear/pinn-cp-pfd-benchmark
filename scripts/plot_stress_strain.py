#!/usr/bin/env python3
"""Stress-strain from the elastic regime through plasticity to damage initiation.

One curve per mesh density plus the digitized experimental markers, with the
initiation point on each curve marked. The point of the figure is the
comparison a reader has to make to judge convergence: do consecutive
densities lie on top of each other, and does initiation land at the
experimental strain-to-damage?
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

STUDY = Path(__file__).resolve().parent.parent


def _curve(path: Path):
    eps, sig = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            eps.append(float(row["boundary_strain"]))
            sig.append(float(row["avg_stress_axial"]))
    order = np.argsort(eps)
    return np.asarray(eps)[order] * 100.0, np.asarray(sig)[order]


def _markers(path: Path):
    x, y = [], []
    with open(path) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            x.append(float(row["strain_pct"]))  # already percent
            y.append(float(row["stress_MPa"]))
    return np.asarray(x), np.asarray(y)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="label=path/to/results.csv")
    ap.add_argument("--markers", default="data/digitized/al2219_fig3_experiment_markers.csv")
    ap.add_argument("--title", default="Calibrated CP set reproduces the measured curve")
    ap.add_argument("--out", default="figures/stress_strain_convergence.png")
    args = ap.parse_args()

    apply_figure_style()
    fig, ax = plt.subplots(figsize=(5.2, 3.9))

    xe, ye = _markers(STUDY / args.markers)
    ax.plot(xe, ye, "o", ms=6, mfc="none", mec="#111111", mew=1.4, zorder=5,
            label="experiment (Vaughan 2024)")

    # Focal series is the finest mesh; coarser ones recede. Threading is by
    # density, so the same density keeps its color in every later figure.
    cmap = plt.get_cmap("viridis")
    specs = [s.split("=", 1) for s in args.runs]
    n = len(specs)
    for i, (label, rel) in enumerate(specs):
        p = STUDY / rel
        if not p.exists():
            continue
        eps, sig = _curve(p)
        focal = i == n - 1
        ax.plot(eps, sig, "-", lw=2.0 if focal else 1.1,
                color=cmap(0.15 + 0.7 * i / max(1, n - 1)),
                alpha=1.0 if focal else 0.75, zorder=4 if focal else 3,
                label=label)
        init = p.parent / "initiation.json"
        if init.exists():
            rec = json.load(open(init))
            ax.plot(rec["boundary_strain"] * 100.0, rec["avg_stress_axial_MPa"],
                    marker="v", ms=7, color=cmap(0.15 + 0.7 * i / max(1, n - 1)),
                    mec="white", mew=0.8, ls="none", zorder=6)

    # The experimental strain-to-damage: the test was stopped here by
    # catastrophic fracture, so it is the target for Gc, not a fitted value.
    ax.axvline(xe[-1], color="#B4462F", lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.annotate("fracture in test", xy=(xe[-1], ax.get_ylim()[0]),
                xytext=(-6, 14), textcoords="offset points",
                rotation=90, ha="right", va="bottom", color="#B4462F", fontsize=6)

    ax.set_xlabel("engineering strain (%)")
    ax.set_ylabel("volume-averaged axial stress (MPa)")
    ax.set_title(args.title)
    ax.margins(0.04)
    ax.legend(frameon=False, loc="lower right", fontsize=6)
    fig.tight_layout()

    out = STUDY / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    print(json.dumps({"figure": str(out), "series": [s[0] for s in specs],
                      "n_markers": int(len(xe))}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
