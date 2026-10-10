#!/usr/bin/env python3
"""Nominal stress-strain of the polycrystal plate against the two references.

Series, in fixed colour order: PINNACLE polycrystal plate(s) (blue, then aqua
and yellow for further runs), the MOOSE single-crystal reference run
(orange), and the digitized Vaughan (2024) experiment markers (black).
One axis; direct labels at the curve ends; no number on every point.

Usage: plot_plate_curve.py --runs label=results.csv [...] --out figures/x.png
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY = Path(__file__).resolve().parent.parent
SERIES = ["#2a78d6", "#1baf7a", "#eda100", "#e87ba4"]
MOOSE = "#eb6834"


def read(path, xcol, ycol):
    rows = list(csv.DictReader(open(path)))
    return [float(r[xcol]) for r in rows], [float(r[ycol]) for r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="label=path/to/results.csv")
    ap.add_argument("--moose", default="", help="optional MOOSE reference CSV; off by default")
    ap.add_argument("--markers", default=str(STUDY / "data/digitized/al2219_fig3_experiment_markers.csv"))
    ap.add_argument("--xmax", type=float, default=0.0)
    ap.add_argument("--title", default="Polycrystal plate, crystal plasticity only")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    for i, spec in enumerate(a.runs):
        label, path = spec.split("=", 1)
        x, y = read(path, "nominal_strain", "eng_stress")
        ax.plot([100 * v for v in x], y, color=SERIES[i % len(SERIES)], lw=2, label=label)
        short = label.split("(")[0].strip()
        ax.annotate(short, (100 * x[-1], y[-1]), xytext=(-4, 6), textcoords="offset points",
                    fontsize=8, ha="right", va="bottom", color="#333", clip_on=True)
    if a.moose:
        x, y = read(a.moose, "nominal_strain_yy", "mean_stress_yy")
        ax.plot([100 * v for v in x], y, color=MOOSE, lw=2, label="MOOSE single crystal (0,0,0), reference deck")
        ax.annotate("MOOSE single crystal", (100 * x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                    fontsize=8, va="center", color="#333")
    if Path(a.markers).exists():
        lines = [ln for ln in open(a.markers) if not ln.startswith("#")]
        rows = list(csv.DictReader(lines))
        ex = [float(r["strain_pct"]) for r in rows]; ey = [float(r["stress_MPa"]) for r in rows]
        if a.xmax:
            ex, ey = zip(*[(u, v) for u, v in zip(ex, ey) if u <= 100 * a.xmax]) if any(u <= 100 * a.xmax for u in ex) else ([], [])
        ax.plot(ex, ey, "o", ms=6, color="#222", mfc="white", mew=1.5, label="Experiment (Vaughan 2024, Fig. 3, digitized)")
    if a.xmax:
        ax.set_xlim(0, 100 * a.xmax)
    ax.set_xlabel("Nominal strain along the pull axis (%)")
    ax.set_ylabel("Engineering axial stress (MPa)")
    ax.set_title(a.title, fontsize=10, loc="left")
    ax.grid(True, color="#e6e6e6", lw=0.6); ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout(); Path(a.out).parent.mkdir(parents=True, exist_ok=True); fig.savefig(a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
