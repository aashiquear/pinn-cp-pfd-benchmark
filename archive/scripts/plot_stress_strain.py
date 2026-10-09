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
    """Strain (%) and axial stress, truncated at the first INVALID step.

    A run that loses the boundary-strain invariant writes rows past that point
    that are not solutions -- strain stops tracking the prescribed value and
    the stress can go negative or absurd. Reading past it produces a curve
    that looks like physics and is not; that mistake has already been made
    once in this study. Validity is inferred from the DATA (strain must track
    load_factor times the program) rather than from run.json, because the
    partial record a stalled run writes omits the strain program -- metadata
    is missing exactly when it matters.
    """
    rows = list(csv.DictReader(open(path)))
    eps = np.array([float(r["boundary_strain"]) for r in rows])
    sig = np.array([float(r["avg_stress_axial"]) for r in rows])
    lam = np.array([float(r["load_factor"]) for r in rows])
    good = np.ones(len(rows), dtype=bool)
    pos = lam > 1.0e-12
    if pos.any():
        target = float(np.median(eps[pos] / lam[pos]))
        dev = np.abs(eps - lam * target)
        good = dev <= np.maximum(1.0e-6, 1.0e-3 * np.abs(lam * target))
    if not good.all():
        first_bad = int(np.argmin(good))
        eps, sig = eps[:first_bad], sig[:first_bad]
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
    ap.add_argument("--inset", type=float, default=0.0,
                    help="Upper strain limit (%) of the elastic-regime inset; "
                         "0 disables it. Graded steps make this the only way "
                         "to read the elastic branch on a 20% program.")
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
    # Elastic-regime inset. Graded load steps put the elastic samples between
    # 0.004% and 0.2% strain, which on an axis spanning 20% is a single pixel
    # against the origin -- the elastic branch and the knee are present in the
    # data but invisible. The inset is the only way to actually read them.
    axi = None
    if args.inset > 0.0:
        axi = ax.inset_axes((0.56, 0.13, 0.40, 0.38))
        axi.set_title("elastic regime", fontsize=7, pad=2)
        axi.tick_params(labelsize=6)
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
        if axi is not None:
            k = eps <= args.inset
            axi.plot(eps[k], sig[k], "-o", ms=2.6, lw=1.2,
                     color=cmap(0.15 + 0.7 * i / max(1, n - 1)))
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

    if axi is not None:
        axi.set_xlabel("strain (%)", fontsize=6, labelpad=1)
        axi.set_ylabel("stress (MPa)", fontsize=6, labelpad=1)
        axi.grid(alpha=0.25, lw=0.4)

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
