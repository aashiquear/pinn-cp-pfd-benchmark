#!/usr/bin/env python3
"""Case C5 (nucleation and growth of spherulites): transformation kinetics.

Usage: python3 scripts/cryst_kinetics_figures.py <tag> [t_start=<time of the quench>] <run dir> [<run dir> ...]
  run dirs relative to the repo; each holds history.csv (step, t, solid_fraction) and frames.npz
  (from cryst_orientation_reference.py, or cryst_binary_collect.py for a PINNACLE run; a
  PINNACLE run's solver_history.csv is used for the fine time series).

The crystalline fraction X(t) is averaged over the runs (realizations of the noise) and fitted
with the Johnson-Mehl-Avrami-Kolmogorov law X = 1 - exp{-[(t - t0)/tau]^p} for 0.01 < X < 0.95,
as in (granasy 2005 spherulites). Also reported: the fit of every single run, the local
exponent d ln[-ln(1 - X)]/d ln(t - t0) against X, and the number of crystals (connected
regions of the last frame in which the fraction is below 0.6).
Writes results/crystallization/orientation/<tag>_kinetics.csv and
figures/crystallization/<tag>_kinetics.png, <tag>_maps.png.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cryst_orientation_figures import orientation_image  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures/crystallization"
PAPER = {"p": 3.04, "tau": 0.0106, "t0": 0.00178}


def series(run: Path) -> tuple[np.ndarray, np.ndarray]:
    """(t, X): fraction above the background of the noisy liquid."""
    solver = run / "solver_history.csv"
    if solver.exists():
        rows = list(csv.DictReader(solver.open()))
        t = np.array([float(r["t"]) for r in rows])
        x = 1.0 - np.array([float(r["phi_mean"]) for r in rows])
    else:
        rows = list(csv.DictReader((run / "history.csv").open()))
        t = np.array([float(r["t"]) for r in rows])
        x = np.array([float(r["solid_fraction"]) for r in rows])
    background = float(np.median(x[: max(3, len(x) // 50)]))
    full = float(x[-1]) if x[-1] > 0.9 else 1.0
    return t, np.clip((x - background) / (full - background), 0.0, 1.0)


def jmak(t, t0, tau, p):
    return 1.0 - np.exp(-(np.clip(t - t0, 1e-12, None) / tau) ** p)


def fit(t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sel = (x > 0.01) & (x < 0.95)
    t_half = float(np.interp(0.5, x, t))
    start = [0.3 * t[sel][0], t_half - 0.3 * t[sel][0], 3.0]
    values, cov = curve_fit(jmak, t[sel], x[sel], p0=start,
                            bounds=([0.0, 1e-6, 0.5], [t[sel][0], 10 * t_half, 12.0]), maxfev=20000)
    return values, np.sqrt(np.diag(cov))


def crystals(run: Path) -> tuple[int, int]:
    """(step, count) of connected crystals in the last frame with a fraction below 0.6."""
    frames = np.load(run / "frames.npz")
    best = (0, 0)
    for i, step in enumerate(frames["steps"]):
        solid = frames["phi"][i] < 0.5
        if 0.02 < solid.mean() < 0.6:
            best = (int(step), int(ndimage.label(solid)[1]))
    return best


def main() -> None:
    tag = sys.argv[1]
    t_start = next((float(a.split("=")[1]) for a in sys.argv[2:] if a.startswith("t_start=")), 0.0)
    runs = [ROOT / a for a in sys.argv[2:] if "=" not in a]
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    data = []
    for r in runs:  # time counted from the quench (t_start = duration of a pre-heat stage)
        ti, xi = series(r)
        keep = ti > t_start
        ti, xi = ti[keep] - t_start, xi[keep]
        base = float(np.median(xi[: max(3, len(xi) // 50)]))  # level of the noisy liquid after the quench
        data.append((ti, np.clip((xi - base) / (1.0 - base), 0.0, 1.0)))
    t = data[0][0]
    mean = np.mean([np.interp(t, ti, xi) for ti, xi in data], axis=0)
    rows = []
    for run, (ti, xi) in zip(runs, data):
        (t0, tau, p), err = fit(ti, xi)
        step, count = crystals(run)
        rows.append({"run": run.name, "t0": t0, "tau": tau, "p": p, "p_error": err[2],
                     "crystals": count, "counted_at_step": step})
    (t0, tau, p), err = fit(t, mean)
    rows.append({"run": "ensemble mean", "t0": t0, "tau": tau, "p": p, "p_error": err[2],
                 "crystals": sum(r["crystals"] for r in rows), "counted_at_step": 0})
    out = ROOT / "results/crystallization/orientation"
    with (out / f"{tag}_kinetics.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(f"{r['run']:42s} t0 {r['t0']:.5f} tau {r['tau']:.5f} p {r['p']:.2f} +/- {r['p_error']:.2f} "
              f"crystals {r['crystals']}")
    print(f"paper: t0 {PAPER['t0']}, tau {PAPER['tau']}, p {PAPER['p']} +/- 0.02 (5000 x 5000)")
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.1))
    for ti, xi in data:
        axes[0].plot(ti, xi, color="0.7", lw=0.8)
    axes[0].plot(t, mean, color="#1b3a5c", lw=1.4, label="mean of the runs")
    axes[0].plot(t, jmak(t, t0, tau, p), "--", color="#c8553d", lw=1.1, label=f"fit, $p$ = {p:.2f}")
    axes[0].set_xlabel("time"); axes[0].set_ylabel("crystalline fraction $X$"); axes[0].legend(frameon=False, fontsize=7)
    sel = (mean > 0.005) & (mean < 0.995) & (t > t0 + 1e-9)
    lt, ly = np.log(t[sel] - t0), np.log(-np.log(1.0 - mean[sel]))
    axes[1].plot(lt, ly, color="#1b3a5c", lw=1.3)
    axes[1].plot(lt, p * (lt - np.log(tau)), "--", color="#c8553d", lw=1.0)
    axes[1].set_xlabel(r"$\ln(t - t_0)$"); axes[1].set_ylabel(r"$\ln[-\ln(1 - X)]$")
    width = max(5, sel.sum() // 12)
    local = np.array([np.polyfit(lt[i:i + width], ly[i:i + width], 1)[0] for i in range(0, len(lt) - width)])
    axes[2].plot(mean[sel][width // 2: width // 2 + len(local)], local, color="#1b3a5c", lw=1.3)
    axes[2].axhline(3.0, color="0.5", ls=":", lw=0.9)
    axes[2].set_ylim(0, 6); axes[2].set_xlabel("crystalline fraction $X$"); axes[2].set_ylabel("local exponent $p$")
    fig.tight_layout()
    fig.savefig(FIG / f"{tag}_kinetics.png", bbox_inches="tight")
    plt.close(fig)
    frames = np.load(runs[0] / "frames.npz")
    picks = [i for i in range(len(frames["steps"])) if 0.02 < (frames["phi"][i] < 0.5).mean()][:3] or [len(frames["steps"]) - 1]
    fig, axes = plt.subplots(1, len(picks), figsize=(3.3 * len(picks), 3.5), squeeze=False)
    for ax, i in zip(axes[0], picks):
        ax.imshow(orientation_image(frames["phi"][i], frames["theta"][i]), origin="lower")
        ax.set_title(f"step {int(frames['steps'][i])}, X = {(frames['phi'][i] < 0.5).mean():.2f}", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(FIG / f"{tag}_maps.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote", f"figures/crystallization/{tag}_kinetics.png, {tag}_maps.png")


if __name__ == "__main__":
    main()
