#!/usr/bin/env python3
"""Case E1 (planar front, solute trapping): measures and figures.

Usage: python3 scripts/cryst_trapping_figures.py
Reads results/crystallization/trapping/e1_planar_s*_fem/; writes summary.csv there and
figures/crystallization/e1_trapping.png.

Measures per run (dx = 6.25e-3, dt = 4.75e-6):
  growth exponent     slope of log(front position - start) against log(time), second half
  front speed         cells per step at the end of the run
  crystal composition mean c of the crystal formed in the second half of the run
  trapped fraction    (c_crystal - c_S) / (c_inf - c_S): 0 at the solidus (full partitioning),
                      1 when the liquid composition is frozen in
  solute peak         largest c ahead of the front at the end, as (c_peak - c_inf)/(c_L - c_inf)
  partition coefficient  c_crystal / c_peak (equilibrium: c_S / c_L = 0.856)
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/trapping"
FIG = ROOT / "figures/crystallization"
C_L, C_S, DX = 0.466219, 0.399112, 6.25e-3
COLORS = ["#1b3a5c", "#2a78a8", "#3aa39b", "#8fb84e", "#e0a030", "#c8553d", "#7a3b69"]


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    runs = sorted(RES.glob("e1_planar_s*_fem"), key=lambda p: float(re.search(r"_s([0-9.]+)_", p.name).group(1)))
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2))
    for k, run in enumerate(runs):
        color = plt.get_cmap("viridis")(k / max(len(runs) - 1, 1))
        s = float(re.search(r"_s([0-9.]+)_", run.name).group(1))
        c_inf = C_L - s * (C_L - C_S)
        hist = list(csv.DictReader((run / "solver_history.csv").open()))
        step = np.array([float(r["step"]) for r in hist]) + 1
        data = read_vtk(run / "fields.vtk")
        order = np.argsort(data["nodes"][:, 0])
        x = data["nodes"][order, 0] / DX
        point = data["point"]
        phi = point["field_phi" if "field_phi" in point else "phi"][order]
        c = point["field_c" if "field_c" in point else "c"][order]
        n_cells = x[-1]
        front = (1.0 - np.array([float(r["phi_mean"]) for r in hist])) * n_cells
        half = len(step) // 2
        expo = float(np.polyfit(np.log(step[half:]), np.log(front[half:] - front[0] + 1e-9), 1)[0])
        speed = float((front[-1] - front[-len(step) // 10]) / (len(step) // 10))
        grown = (x > max(front[half], front[0] + 12)) & (x < front[-1] - 10) & (phi < 0.01)
        c_crystal = float(c[grown].mean()) if grown.any() else float("nan")
        ahead = x > front[-1] - 15
        c_peak = float(c[ahead].max())
        rows.append({"S": s, "c_inf": c_inf, "front_cells": float(front[-1]), "growth_exponent": expo,
                     "speed_cells_per_step": speed, "c_crystal": c_crystal,
                     "trapped_fraction": (c_crystal - C_S) / (c_inf - C_S) if abs(c_inf - C_S) > 1e-9 else float("nan"),
                     "solute_peak": (c_peak - c_inf) / (C_L - c_inf), "partition_coefficient": c_crystal / c_peak})
        axes[0].loglog(step, np.maximum(front - front[0], 1e-2), color=color, lw=1.2, label=f"S = {s:g}")
        axes[1].plot(x - front[-1], c, color=color, lw=1.1)
    axes[0].set_xlabel("step"); axes[0].set_ylabel("front advance (cells)"); axes[0].legend(frameon=False, fontsize=6, ncol=2)
    axes[1].set_xlim(-150, 150); axes[1].axhline(C_S, color="0.5", ls=":", lw=0.8); axes[1].axhline(C_L, color="0.5", ls=":", lw=0.8)
    axes[1].set_xlabel("distance from the front (cells)"); axes[1].set_ylabel("composition $c$")
    s_values = [r["S"] for r in rows]
    axes[2].plot(s_values, [r["partition_coefficient"] for r in rows], "o-", color="#1b3a5c", ms=4, lw=1.1,
                 label=r"partition coefficient $c_{cr}/c_{peak}$")
    axes[2].plot(s_values, [r["solute_peak"] for r in rows], "s-", color="#c8553d", ms=4, lw=1.1,
                 label="solute peak (0: none, 1: liquidus)")
    axes[2].axhline(C_S / C_L, color="0.5", ls=":", lw=0.8)
    axes[2].set_xlabel("supersaturation $S$"); axes[2].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "e1_trapping.png", bbox_inches="tight")
    plt.close(fig)
    with (RES / "summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={v:.4g}" for k, v in r.items()))
    print("wrote figures/crystallization/e1_trapping.png")


if __name__ == "__main__":
    main()
