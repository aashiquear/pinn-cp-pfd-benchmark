#!/usr/bin/env python3
"""Transformation kinetics of case B2 (noise-induced nucleation, growth, impingement).

Usage: python3 scripts/cryst_noise_figures.py [run directory name]
Reads results/crystallization/noise/<run>/ (the fields_step*.vtk series and
solver_history.csv); writes
  results/crystallization/noise/<run>_kinetics.csv   t, X (threshold), X (integral), grains
  results/crystallization/noise/<run>_jmak.json      fit of X = 1 - exp{-[(t - t0)/tau]^p}
  figures/crystallization/b2_jmak_snapshots.png, b2_jmak_kinetics.png, b2_jmak.gif

X is the area fraction with phi > 0.5 (nodal values weighted by the nodal area). The fit uses
0.01 < X < 0.95, as in (granasy 2005 spherulites), section III.B. The uncertainty of p is the
standard error of the nonlinear least-squares fit; the spread over realizations needs more
seeds and is reported separately when they are run.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from scipy import ndimage  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from vtk_legacy import read_vtk, series_files  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.grid": False, "figure.dpi": 150, "savefig.dpi": 200,
                     "axes.linewidth": 0.8, "legend.frameon": False, "mathtext.fontset": "cm"})
run_name = sys.argv[1] if len(sys.argv) > 1 else "b2_jmak_2d_fem"
RUN = ROOT / "results/crystallization/noise" / run_name
FIG = ROOT / "figures/crystallization"
tag = run_name.replace("_2d_fem", "")


def jmak(t, t0, tau, p):
    s = np.clip((t - t0) / tau, 0.0, None)
    return 1.0 - np.exp(-(s**p))


def main() -> None:
    times = [float(t) for t in re.findall(r'timestep="([^"]+)"', (RUN / "fields.pvd").read_text())]
    files = series_files(RUN)
    first = read_vtk(files[0][1])
    x, y = first["nodes"][:, 0], first["nodes"][:, 1]
    xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
    ix, iy = np.searchsorted(xs, np.round(x, 9)), np.searchsorted(ys, np.round(y, 9))
    # nodal area weights of the structured mesh (half cells on the boundary)
    wx = np.gradient(xs); wx[[0, -1]] *= 0.5
    wy = np.gradient(ys); wy[[0, -1]] *= 0.5
    area = np.outer(wy, wx)
    grids, frac, grains = [], [], []
    for (_step, path), _t in zip(files, times):
        u = read_vtk(path)["point"]["phi"]
        g = np.zeros((len(ys), len(xs)))
        g[iy, ix] = u
        grids.append(g.astype(np.float32))
        solid = g > 0.5
        frac.append(float(np.sum(area[solid]) / np.sum(area)))
        grains.append(int(ndimage.label(solid)[1]))
    t = np.array(times)
    X = np.array(frac)
    hist = list(csv.DictReader((RUN / "solver_history.csv").open()))
    th = np.array([float(r["t"]) for r in hist])
    xi = np.array([float(r["mass"]) for r in hist]) / float(np.sum(area))
    with (RUN.parent / f"{run_name}_kinetics.csv").open("w") as fh:
        fh.write("t,X_threshold,X_integral,grains\n")
        for a, b, c in zip(t, X, grains):
            fh.write(f"{a:.6g},{b:.6f},{np.interp(a, th, xi) if a > 0 else 0.0:.6f},{c}\n")

    sel = (X > 0.01) & (X < 0.95)
    popt, pcov = curve_fit(jmak, t[sel], X[sel], p0=[t[sel][0] * 0.5, t[sel].mean(), 3.0],
                           bounds=([0.0, 1e-6, 0.5], [t[sel][0], 10 * t[-1], 8.0]))
    perr = np.sqrt(np.diag(pcov))
    t0, tau, p = popt
    fit = {"t0": t0, "tau": tau, "p": p, "t0_se": perr[0], "tau_se": perr[1], "p_se": perr[2],
           "n_points": int(sel.sum()), "max_grains": int(max(grains)),
           "rms_residual": float(np.sqrt(np.mean((jmak(t[sel], *popt) - X[sel]) ** 2)))}
    (RUN.parent / f"{run_name}_jmak.json").write_text(json.dumps(fit, indent=2))
    print(json.dumps(fit, indent=2))

    # ---- snapshots
    picks = [int(np.argmin(np.abs(X - v))) for v in (0.02, 0.15, 0.4, 0.7, 0.92)]
    fig, axes = plt.subplots(1, len(picks), figsize=(12.5, 2.9))
    for ax, k in zip(axes, picks):
        im = ax.imshow(grids[k], origin="lower", extent=(xs[0], xs[-1], ys[0], ys[-1]), cmap="cividis",
                       vmin=0, vmax=1)
        ax.set_title(rf"$t = {t[k]:.1f}$, $X = {X[k]:.2f}$", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.01)
    cb.set_label(r"$\phi$")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{tag}_snapshots.png", bbox_inches="tight")
    plt.close(fig)

    # ---- kinetics
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.2))
    tt = np.linspace(0, t[-1], 400)
    ax[0].plot(t, X, "o", ms=3.5, mfc="none", color="#1f4e79", label="PINNACLE, FE")
    ax[0].plot(tt, jmak(tt, *popt), color="black", lw=1.0,
               label=rf"JMAK fit, $p = {p:.2f} \pm {perr[2]:.2f}$")
    ax[0].set_xlim(0, t[-1]); ax[0].set_ylim(0, 1.02)
    ax[0].set_xlabel(r"$t\,M_\phi w$"); ax[0].set_ylabel(r"$X$"); ax[0].legend(loc="lower right")
    s2 = sel & (t > t0)
    ax[1].plot(np.log(t[s2] - t0), np.log(-np.log(1 - X[s2])), "o", ms=3.5, mfc="none", color="#1f4e79")
    lx = np.log(t[s2] - t0)
    ax[1].plot(lx, p * (lx - np.log(tau)), color="black", lw=1.0, label=rf"slope $p = {p:.2f}$")
    ax[1].plot(lx, 3.0 * (lx - np.log(tau)), "--", color="#888888", lw=0.9, label=r"slope $1 + d = 3$")
    ax[1].set_xlabel(r"$\ln(t - t_0)$"); ax[1].set_ylabel(r"$\ln[-\ln(1 - X)]$"); ax[1].legend()
    ax[2].plot(t, grains, color="#c0392b", lw=1.2)
    ax[2].set_xlim(0, t[-1]); ax[2].set_ylim(0, None)
    ax[2].set_xlabel(r"$t\,M_\phi w$"); ax[2].set_ylabel(r"separate crystal regions ($\phi > 0.5$)")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / f"{tag}_kinetics.png", bbox_inches="tight")
    plt.close(fig)

    # ---- animation
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    im = ax.imshow(grids[0], origin="lower", extent=(xs[0], xs[-1], ys[0], ys[-1]), cmap="cividis", vmin=0, vmax=1)
    ax.set_xticks([]); ax.set_yticks([])
    title = ax.set_title("")

    def draw(i):
        im.set_data(grids[i])
        title.set_text(rf"$t = {t[i]:.1f}$   $X = {X[i]:.2f}$")
        return im, title

    FuncAnimation(fig, draw, frames=len(grids)).save(FIG / f"{tag}.gif", writer=PillowWriter(fps=8), dpi=80)
    plt.close(fig)
    print("wrote", f"figures/crystallization/{tag}_snapshots.png, {tag}_kinetics.png, {tag}.gif")


if __name__ == "__main__":
    main()
