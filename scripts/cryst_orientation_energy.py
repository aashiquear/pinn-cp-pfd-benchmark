#!/usr/bin/env python3
"""Case C1: the orientational free energy of (granasy 2005 spherulites), Eq. 2.1 and Fig. 3.

    f_ori = [1 - p(phi)] (H T / 2 xi0) { x F0 + (1 - x) F1 },   u = xi0 |grad theta|
    F0 = |sin(2 pi m u)|  for u < 3/(4m), else 1
    F1 = |sin(2 pi n u)|  for u < 1/(4n), else 1

theta in [0, 1] is the normalized orientation; for k-fold symmetry the misorientation angle
between neighboring cells is u * 360 deg / k. With n = 1/2 the second term has its only
minimum at zero misorientation; the first adds a metastable minimum at u = 1/(2m), i.e. at
360 / (2 m k) degrees: 30 degrees for k = 2, m = 3 (branching with a fixed misorientation).

Usage: python3 scripts/cryst_orientation_energy.py
Writes figures/crystallization/c1_orientation_energy.png and
results/crystallization/orientation/c1_features.csv (reference values, no PINNACLE run).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.size": 10, "axes.grid": False, "savefig.dpi": 200, "axes.linewidth": 0.8,
                     "legend.frameon": False, "mathtext.fontset": "cm"})


def f_ori(u: np.ndarray, m: float, n: float, x: float) -> np.ndarray:
    """Orientational free energy in units of H T / (2 xi0), in the liquid-side factor 1."""
    f0 = np.where(u < 3.0 / (4.0 * m), np.abs(np.sin(2 * np.pi * m * u)), 1.0)
    f1 = np.where(u < 1.0 / (4.0 * n), np.abs(np.sin(2 * np.pi * n * u)), 1.0)
    return x * f0 + (1.0 - x) * f1


def main() -> None:
    k, m, n, x = 2, 3.0, 0.5, 0.2  # the parameters of the paper's Fig. 3
    u = np.linspace(0.0, 0.5, 20001)
    angle = u * 360.0 / k
    f = f_ori(u, m, n, x)
    # features named in the caption of Fig. 3
    lo = (angle > 5) & (angle < 29)
    i_max = np.where(lo)[0][np.argmax(f[lo])]
    mid = (angle > 25) & (angle < 35)
    i_min = np.where(mid)[0][np.argmin(f[mid])]
    out = ROOT / "results/crystallization/orientation"
    out.mkdir(parents=True, exist_ok=True)
    rows = [("local maximum", angle[i_max], f[i_max], "about 20 deg (caption)"),
            ("metastable minimum", angle[i_min], f[i_min], "30 deg (caption)"),
            ("value at 90 deg", 90.0, f[-1], "1 (axis of the figure)"),
            ("value at 0 deg", 0.0, f[0], "0")]
    with (out / "c1_features.csv").open("w") as fh:
        fh.write("feature,angle_deg,f_ori_over_HT_2xi0,published\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]:.3f},{r[2]:.5f},{r[3]}\n")
            print(r)
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    ax.plot(angle, f, color="black", lw=1.2, label=r"$k=2,\ m=3,\ n=1/2,\ x=0.2$")
    for mm, col in ((1.0, "#1f4e79"), (2.0, "#c0392b")):
        ax.plot(angle, f_ori(u, mm, n, x), color=col, lw=0.9, ls="--", label=rf"$m={mm:g}$")
    ax.set_xlim(0, 90); ax.set_ylim(0, 1.2)
    ax.set_xlabel(r"$\xi_0|\nabla\theta|\,2\pi/k$ (deg)")
    ax.set_ylabel(r"$f_{ori}$ $[HT/2\xi_0]$")
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(ROOT / "figures/crystallization/c1_orientation_energy.png", bbox_inches="tight")
    print("wrote figures/crystallization/c1_orientation_energy.png")


if __name__ == "__main__":
    main()
