#!/usr/bin/env python3
"""Case F1 (hard-sphere nucleus): collect the PINNACLE runs, compare with the reference, plot.

Usage: python3 scripts/cryst_hard_sphere_figures.py
Reads results/crystallization/hard_sphere/{reference.csv, reference_profiles.npz,
f1_nucleus_*_fem/}; writes summary.csv there and figures/crystallization/f1_*.png.
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
from cryst_binary_collect import read_vtk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/hard_sphere"
FIG = ROOT / "figures/crystallization"
COLORS = {"average": "#1b3a5c", "111": "#3aa39b", "110": "#c8553d", "100": "#e0a030"}
# Monte Carlo barriers of Auer and Frenkel as read by eye from Fig. 5 of the paper's preprint
# (`digitized`, about +/- 1 kT): volume fraction, W*/kT.
MONTE_CARLO = [(0.5207, 42.0), (0.5277, 26.0), (0.5343, 19.0)]


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    ref = list(csv.DictReader((RES / "reference.csv").open()))
    profiles = np.load(RES / "reference_profiles.npz")
    rows = []
    for r in ref:
        run = RES / f"f1_nucleus_{r['interface']}_phi{float(r['phi_inf']):.4f}_fem"
        hist = list(csv.DictReader((run / "solver_history.csv").open()))[-1]
        data = read_vtk(run / "fields.vtk")
        order = np.argsort(data["nodes"][:, 0])
        x, m, phi = data["nodes"][order, 0], data["point"]["m"][order], data["point"]["phi"][order]
        rr, mr, pr = profiles[f"{r['interface']}_{float(r['phi_inf']):.4f}"]
        inside = x <= rr[-1]
        rows.append({"interface": r["interface"], "phi_inf": float(r["phi_inf"]), "W_pinnacle": float(hist["energy"]),
                     "W_reference": float(r["W_pft"]), "W_cnt": float(r["W_cnt"]),
                     "rel_diff": float(hist["energy"]) / float(r["W_pft"]) - 1,
                     "newton_iters": int(float(hist["newton_iters"])), "m_center": float(m[0]),
                     "phi_center": float(phi[0]),
                     "max_profile_diff_m": float(np.max(np.abs(np.interp(x[inside], rr, mr) - m[inside]))),
                     "radius_half": float(np.interp(0.5 * (1 + m[0]), m, x)),
                     "excess_particles": float(np.trapezoid(4 * np.pi * x**2 * (phi - float(r["phi_inf"])) * 6 / np.pi, x))})
        rows[-1]["_profile"] = (x, m, phi)
    with (RES / "summary.csv").open("w", newline="") as fh:
        keys = [k for k in rows[0] if not k.startswith("_")]
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(f"{r['interface']:8s} {r['phi_inf']:.4f} W {r['W_pinnacle']:8.3f} ref {r['W_reference']:8.3f} "
              f"({r['rel_diff']:+.1e}) CNT {r['W_cnt']:7.3f} m0 {r['m_center']:.4f} R_half {r['radius_half']:.2f} "
              f"dN {r['excess_particles']:.1f}")
    # ---- barrier against volume fraction
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.3))
    ax = axes[0]
    for name, color in COLORS.items():
        sel = [r for r in rows if r["interface"] == name]
        x = [r["phi_inf"] for r in sel]
        ax.plot(x, [r["W_reference"] for r in sel], "-", color=color, lw=1.1)
        ax.plot(x, [r["W_pinnacle"] for r in sel], "o", color=color, ms=4, mfc="none",
                label=f"{name}" if name == "average" else f"({name})")
    sel = [r for r in rows if r["interface"] == "average"]
    ax.plot([r["phi_inf"] for r in sel], [r["W_cnt"] for r in sel], "--", color="0.4", lw=1.0, label="classical, average")
    ax.plot(*zip(*MONTE_CARLO), "s", color="k", ms=5, label="Monte Carlo (digitized)")
    ax.plot([0.5207], [40.63], "*", color="#c8553d", ms=9, label="paper, 40.63")
    ax.set_yscale("log"); ax.set_xlabel(r"liquid volume fraction $\phi_\infty$"); ax.set_ylabel(r"$W^*/kT$")
    ax.legend(frameon=False, fontsize=6.5)
    ax = axes[1]
    ax.plot([r["phi_inf"] for r in sel], [r["W_pinnacle"] / r["W_cnt"] for r in sel], "o-", color=COLORS["average"],
            ms=4, lw=1.1, label=r"$W^*/W^*_{\rm classical}$")
    ax.plot([r["phi_inf"] for r in sel], [r["m_center"] for r in sel], "s-", color=COLORS["110"], ms=4, lw=1.1,
            label=r"$m(0)$")
    ax.set_xlabel(r"liquid volume fraction $\phi_\infty$"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "f1_barrier.png", bbox_inches="tight"); plt.close(fig)
    # ---- profiles at the three Monte Carlo states (average interface)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.0))
    for (phi_inf, _), color in zip(MONTE_CARLO, ("#1b3a5c", "#3aa39b", "#c8553d")):
        r = next(q for q in rows if q["interface"] == "average" and abs(q["phi_inf"] - phi_inf) < 1e-9)
        x, m, phi = r["_profile"]
        rr, mr, pr = profiles[f"average_{phi_inf:.4f}"]
        axes[0].plot(rr, 1 - mr, "-", color=color, lw=2.2, alpha=0.35)
        axes[0].plot(x, 1 - m, "-", color=color, lw=1.0, label=rf"$\phi_\infty = {phi_inf}$")
        axes[1].plot(rr, pr, "-", color=color, lw=2.2, alpha=0.35)
        axes[1].plot(x, phi, "-", color=color, lw=1.0)
    for ax, lab in zip(axes, (r"crystallinity $1 - m$", r"volume fraction $\phi$")):
        ax.set_xlim(0, 20); ax.set_xlabel(r"$r/\sigma$"); ax.set_ylabel(lab)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "f1_profiles.png", bbox_inches="tight"); plt.close(fig)
    print("wrote figures/crystallization/f1_barrier.png, f1_profiles.png")


if __name__ == "__main__":
    main()
