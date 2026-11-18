#!/usr/bin/env python3
"""Case F3 (ternary eutectic growth at a fixed driving force): measures and figures.

Usage: python3 scripts/cryst_ternary_figures.py [long]   (long: only the long runs)
Reads results/crystallization/ternary/f3_ternary_fs*_fem/; writes summary.csv there and
figures/crystallization/f3_maps.png, f3_measures.png, f3_<run>.gif. The long runs
(f3_long_fs*_fem, 600,000 steps) give f3_long.png and long_summary.csv: maps at the end, the
front position against time with the fitted growth exponent (front - seed ~ t^n over the last
80 % of the run; n = 1 is steady growth, n = 1/2 is growth limited by diffusion), and the
profile of c3 along the growth direction at the end.

Measures per run (dx = 1, dt = 0.0025, as in the paper's dimensionless units):
  front speed          slope of the solid fraction over the last 40 % of the run times the box
                       length (cells per unit time); "steady": relative change between the last
                       two fifths
  lamella pairs        number of pairs across the box in the solid 12 cells behind the front
                       (sign changes of c1 - c2 divided by two); spacing = width / pairs
  contrast             mean |c1 - c2| of the solid grown in the run (0.735 for the two solid
                       solutions of the start, 0 for a homogeneous solid)
  c3 in the solid      mean of the solid grown in the run, and its ratio to the liquid value
                       (0.09): 1 is complete trapping
  c3 peak              largest c3 in the liquid ahead of the front
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
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk, series_files  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/ternary"
FIG = ROOT / "figures/crystallization"
SEED = 40.0


def grid(path: Path):
    d = read_vtk(path)
    n = d["nodes"]
    xs, ys = np.unique(np.round(n[:, 0], 6)), np.unique(np.round(n[:, 1], 6))
    ix, iy = np.searchsorted(xs, np.round(n[:, 0], 6)), np.searchsorted(ys, np.round(n[:, 1], 6))
    out = []
    for key in ("phi", "c1", "c2"):
        a = np.zeros((len(ys), len(xs)))
        a[iy, ix] = d["point"][key if key in d["point"] else f"field_{key}"]
        out.append(a)
    return xs, ys, out[0], out[1], out[2]


def analyse(run: Path) -> dict:
    fs = float(re.search(r"_fs(-[0-9.]+)_", run.name).group(1))
    rows = list(csv.DictReader((run / "solver_history.csv").open()))
    t = np.array([float(r["t"]) for r in rows])
    solid = 1.0 - np.array([float(r["phi_mean"]) for r in rows])
    xs, ys, phi, c1, c2 = grid(run / "fields.vtk")
    length, width = float(ys[-1]), float(xs[-1])
    n = len(t)
    full = np.argmax(solid > 0.97) if (solid > 0.97).any() else n  # the front reached the top
    n_eff = max(int(full), 10)

    def slope(a: float, b: float) -> float:
        i, j = int(a * n_eff), int(b * n_eff) - 1
        return float((solid[j] - solid[i]) / (t[j] - t[i]) * length)

    front = float(solid[-1] * length)
    c3 = 1.0 - c1 - c2
    grown = (ys[:, None] > SEED + 15.0) & (ys[:, None] < front - 12.0) & (phi < 0.05)
    row = int(np.clip(np.searchsorted(ys, front - 12.0), 0, len(ys) - 1))
    diff = (c1 - c2)[row, :-1]
    pairs = int(np.sum(np.sign(diff) != np.sign(np.roll(diff, 1))) // 2)
    liquid = phi > 0.5
    return {"f_s": fs, "t_end": float(t[min(n_eff, n) - 1]), "front": front, "reached_top": bool(full < n),
            "speed": slope(0.6, 1.0), "steady": abs(slope(0.8, 1.0) / slope(0.6, 0.8) - 1.0),
            "lamella_pairs": pairs, "spacing": width / pairs if pairs else float("nan"),
            "contrast": float(np.abs(c1 - c2)[grown].mean()) if grown.any() else float("nan"),
            "c3_solid": float(c3[grown].mean()) if grown.any() else float("nan"),
            "c3_trapped": float(c3[grown].mean() / 0.09) if grown.any() else float("nan"),
            "c3_peak": float(c3[liquid].max()) if liquid.any() else float("nan"),
            "_fields": (xs, ys, phi, c1, c2), "_run": run}


def long_runs() -> None:
    runs = sorted((p for p in RES.glob("f3_long_fs*_fem") if (p / "solver_history.csv").exists()),
                  key=lambda p: -float(re.search(r"_fs(-[0-9.]+)_", p.name).group(1)))
    if not runs:
        return
    rows = [analyse(r) for r in runs]
    n = len(rows)
    fig = plt.figure(figsize=(1.25 * 2 * n + 7.4, 5.0))
    gs = fig.add_gridspec(1, 2 * n + 4, width_ratios=[1.0] * (2 * n) + [0.5, 2.6, 0.7, 2.6], wspace=0.1)
    ax_t, ax_c = fig.add_subplot(gs[0, 2 * n + 1]), fig.add_subplot(gs[0, 2 * n + 3])
    colors = ["#1b3a5c", "#c8553d", "#3aa39b", "#8a6fb0"]
    for j, r in enumerate(rows):
        xs, ys, phi, c1, c2 = r["_fields"]
        hist = list(csv.DictReader((r["_run"] / "solver_history.csv").open()))
        t = np.array([float(h["t"]) for h in hist])
        front = (1.0 - np.array([float(h["phi_mean"]) for h in hist])) * float(ys[-1])
        sel = (t > 0.2 * t[-1]) & (front - SEED > 1.0)
        r["exponent"] = float(np.polyfit(np.log(t[sel]), np.log(front[sel] - SEED), 1)[0])
        a1, a2 = fig.add_subplot(gs[0, 2 * j]), fig.add_subplot(gs[0, 2 * j + 1])
        a1.imshow(np.where(phi < 0.5, c1 - c2, np.nan), origin="lower", cmap="RdBu_r", vmin=-0.8, vmax=0.8,
                  extent=[0, xs[-1], 0, ys[-1]])
        a2.imshow(1.0 - c1 - c2, origin="lower", cmap="magma", vmin=0.0, vmax=0.2, extent=[0, xs[-1], 0, ys[-1]])
        a1.set_title(f"$f_s$ = {r['f_s']:g}\n$c_1 - c_2$", fontsize=8)
        a2.set_title("$c_3$", fontsize=8)
        for a in (a1, a2):
            a.set_xticks([]); a.set_yticks([])
        ax_t.plot(t[::200], front[::200], color=colors[j], lw=1.2,
                  label=f"$f_s$ = {r['f_s']:g}, $n$ = {r['exponent']:.2f}")
        ax_c.plot(ys, (1.0 - c1 - c2).mean(axis=1), color=colors[j], lw=1.2, label=f"$f_s$ = {r['f_s']:g}")
    ax_t.set_xlabel("time"); ax_t.set_ylabel("front position (cells)")
    ax_t.legend(frameon=False, fontsize=7)
    ax_c.axhline(0.09, color="0.5", ls=":", lw=0.9)
    ax_c.set_xlabel("position along the growth direction (cells)")
    ax_c.set_ylabel("$c_3$, mean across the box")
    ax_c.legend(frameon=False, fontsize=7)
    fig.savefig(FIG / "f3_long.png", bbox_inches="tight")
    plt.close(fig)
    keys = [k for k in rows[0] if not k.startswith("_")]
    with (RES / "long_summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print("long", " ".join(f"{k}={r[k]:.4g}" if isinstance(r[k], float) else f"{k}={r[k]}" for k in keys))
    print("wrote figures/crystallization/f3_long.png")


def main() -> None:
    if "long" in sys.argv[1:]:
        plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
        long_runs()
        return
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    runs = sorted((p for p in RES.glob("f3_ternary_fs*_fem") if (p / "solver_history.csv").exists()),
                  key=lambda p: -float(re.search(r"_fs(-[0-9.]+)_", p.name).group(1)))
    rows = [analyse(r) for r in runs]
    keys = [k for k in rows[0] if not k.startswith("_")]
    with (RES / "summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={r[k]:.4g}" if isinstance(r[k], float) else f"{k}={r[k]}" for k in keys))
    fig, axes = plt.subplots(2, len(rows), figsize=(1.5 * len(rows) + 1.2, 7.0), squeeze=False)
    for j, r in enumerate(rows):
        xs, ys, phi, c1, c2 = r["_fields"]
        im1 = axes[0, j].imshow(np.where(phi < 0.5, c1 - c2, np.nan), origin="lower", cmap="RdBu_r", vmin=-0.8,
                                vmax=0.8, extent=[0, xs[-1], 0, ys[-1]])
        im2 = axes[1, j].imshow(1.0 - c1 - c2, origin="lower", cmap="magma", vmin=0.0, vmax=0.3,
                                extent=[0, xs[-1], 0, ys[-1]])
        axes[0, j].set_title(f"$f_s$ = {r['f_s']:g}", fontsize=8)
        for a in axes[:, j]:
            a.set_xticks([]); a.set_yticks([])
    axes[0, 0].set_ylabel("$c_1 - c_2$ in the solid")
    axes[1, 0].set_ylabel("$c_3$")
    fig.colorbar(im1, ax=axes[0, :], shrink=0.8, pad=0.01)
    fig.colorbar(im2, ax=axes[1, :], shrink=0.8, pad=0.01)
    fig.savefig(FIG / "f3_maps.png", bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.0))
    fs = [-r["f_s"] for r in rows]
    axes[0].semilogy(fs, [r["speed"] for r in rows], "o-", color="#1b3a5c", ms=4, lw=1.1)
    axes[0].axhline(0.35, color="#c8553d", ls=":", lw=1.0)
    axes[0].set_ylabel("front speed (cells per unit time)")
    axes[1].plot(fs, [r["c3_trapped"] for r in rows], "o-", color="#1b3a5c", ms=4, lw=1.1, label="$c_3$ solid / liquid")
    axes[1].plot(fs, [r["contrast"] / 0.735 for r in rows], "s-", color="#c8553d", ms=4, lw=1.1,
                 label="phase contrast (1: two solids)")
    axes[1].legend(frameon=False, fontsize=7)
    axes[2].plot(fs, [r["c3_peak"] for r in rows], "o-", color="#3aa39b", ms=4, lw=1.1)
    axes[2].set_ylabel("$c_3$ peak ahead of the front")
    for ax in axes:
        ax.set_xlabel("driving force $-f_s$")
    fig.tight_layout()
    fig.savefig(FIG / "f3_measures.png", bbox_inches="tight")
    plt.close(fig)
    for r in rows:
        if abs(r["f_s"] + 0.86) > 1e-9:
            continue
        frames = [grid(path) for _, path in series_files(r["_run"])]
        fig, ax = plt.subplots(figsize=(2.4, 4.6))
        xs, ys = frames[0][0], frames[0][1]
        im = ax.imshow(frames[0][3] - frames[0][4], origin="lower", cmap="RdBu_r", vmin=-0.8, vmax=0.8,
                       extent=[0, xs[-1], 0, ys[-1]])
        ax.set_xticks([]); ax.set_yticks([])

        def draw(i):
            im.set_data(frames[i][3] - frames[i][4])
            return (im,)

        FuncAnimation(fig, draw, frames=len(frames)).save(FIG / f"f3_{r['_run'].name}.gif",
                                                          writer=PillowWriter(fps=4), dpi=80)
        plt.close(fig)
    print("wrote figures/crystallization/f3_maps.png, f3_measures.png")


if __name__ == "__main__":
    main()
