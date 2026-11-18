#!/usr/bin/env python3
"""Case F2 (lamellar eutectic growth): front speed, grown structure and figures.

Usage: python3 scripts/cryst_eutectic_figures.py
Reads results/crystallization/eutectic/f2_lamellar_*_fem/ (solver_history.csv, fields.vtk and the
VTK series); writes summary.csv there and figures/crystallization/f2_*.png and f2_*.gif.

Measures per run:
  front speed         slope of the solid fraction over the last 40 % of the run times the box
                      length (nm per time unit; 1 time unit = 7.143e-10 s), also in m/s
  steady              relative change of that slope between the last two fifths of the run
  beta fraction       share of the Cu-rich phase (c > 0.5) in the solid grown during the run
  c_alpha, c_beta     mean compositions of the two grown solids (equilibrium values from the deck)
  lambda^2 v          Jackson-Hunt product in m^3/s
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
from cryst_binary_collect import read_vtk, series_files  # noqa: E402
from cryst_eutectic_thermo import eutectic, groups  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/eutectic"
FIG = ROOT / "figures/crystallization"
TIME_UNIT, SEED_HEIGHT = 7.142857e-10, 10.0


def grid(path: Path):
    d = read_vtk(path)
    n = d["nodes"]
    xs, ys = np.unique(np.round(n[:, 0], 6)), np.unique(np.round(n[:, 1], 6))
    ix, iy = np.searchsorted(xs, np.round(n[:, 0], 6)), np.searchsorted(ys, np.round(n[:, 1], 6))
    out = {}
    for key in ("phi", "c"):
        a = np.zeros((len(ys), len(xs)))
        a[iy, ix] = d["point"][key]
        out[key] = a
    return xs, ys, out["phi"], out["c"]


def analyse(run: Path) -> dict:
    m = re.match(r"f2_lamellar_T([0-9.]+)_c([0-9.]+)_lam([0-9.]+)_fem", run.name)
    t_k, c0, lam = float(m.group(1)), float(m.group(2)), float(m.group(3))
    rows = list(csv.DictReader((run / "solver_history.csv").open()))
    t = np.array([float(r["t"]) for r in rows])
    solid = 1.0 - np.array([float(r["phi_mean"]) for r in rows])
    xs, ys, phi, c = grid(run / "fields.vtk")
    length = float(ys[-1])
    n = len(t)

    def slope(a: float, b: float) -> float:
        i, j = int(a * n), int(b * n) - 1
        return float((solid[j] - solid[i]) / (t[j] - t[i]) * length)

    v = slope(0.6, 1.0)
    front = float(solid[-1] * length)
    grown = (ys[:, None] > SEED_HEIGHT + 2.0) & (ys[:, None] < front - 2.0) & (phi < 0.1)
    beta = grown & (c > 0.5)
    alpha = grown & (c <= 0.5)
    gr = groups(t_k)
    return {"T": t_k, "undercooling": eutectic()[3] - t_k, "c0": c0, "lambda_nm": lam, "t_end": float(t[-1]),
            "front_nm": front, "speed_nm_per_unit": v, "speed_m_s": v * 1e-9 / TIME_UNIT,
            "steady": abs(slope(0.8, 1.0) / slope(0.6, 0.8) - 1.0),
            "lambda2_v_m3_s": (lam * 1e-9) ** 2 * v * 1e-9 / TIME_UNIT,
            "beta_fraction": float(beta.sum() / max(grown.sum(), 1)),
            "c_alpha_grown": float(c[alpha].mean()) if alpha.any() else float("nan"),
            "c_beta_grown": float(c[beta].mean()) if beta.any() else float("nan"),
            "c_alpha_eq": gr["c_alpha"], "c_beta_eq": gr["c_beta"],
            "grown_rows": int(grown.any(axis=1).sum()), "_fields": (xs, ys, phi, c)}


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    runs = sorted(p for p in RES.glob("f2_lamellar_*_fem") if (p / "solver_history.csv").exists())
    rows = [analyse(r) for r in runs]
    keys = [k for k in rows[0] if not k.startswith("_")]
    with (RES / "summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={r[k]:.4g}" for k in keys))

    def panel(sel: list[dict], label, tag: str) -> None:
        fig, axes = plt.subplots(1, len(sel), figsize=(1.0 + 1.5 * len(sel), 3.6), squeeze=False)
        for ax, r in zip(axes[0], sel):
            xs, ys, phi, c = r["_fields"]
            shown = np.where(phi > 0.5, np.nan, c)  # liquid blank
            reps = max(1, int(round(48.0 / xs[-1])))
            tiled = np.tile(shown[:, :-1], (1, reps))
            im = ax.imshow(tiled, origin="lower", extent=[0, reps * xs[-1], 0, ys[-1]], cmap="viridis", vmin=0, vmax=1,
                           aspect="equal")
            ax.contour(np.linspace(0, reps * xs[-1], tiled.shape[1]), ys, np.tile(phi[:, :-1], (1, reps)),
                       levels=[0.5], colors="k", linewidths=0.6)
            ax.set_title(label(r), fontsize=8)
            ax.set_xlabel("nm")
        axes[0, 0].set_ylabel("growth direction (nm)")
        fig.colorbar(im, ax=axes[0, :], shrink=0.7, pad=0.02, label="Cu fraction $c$ (solid)")
        fig.savefig(FIG / f"{tag}.png", bbox_inches="tight")
        plt.close(fig)

    t_ref = 955.5
    spacing = sorted((r for r in rows if abs(r["T"] - t_ref) < 1e-6), key=lambda r: r["lambda_nm"])
    cooling = sorted((r for r in rows if abs(r["lambda_nm"] - 16.0) < 1e-6), key=lambda r: -r["T"])
    if spacing:
        panel(spacing, lambda r: rf"$\lambda$ = {r['lambda_nm']:g} nm", "f2_spacing_maps")
        fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.0))
        lam = [r["lambda_nm"] for r in spacing]
        axes[0].plot(lam, [r["speed_m_s"] for r in spacing], "o-", color="#1b3a5c", ms=4, lw=1.1)
        axes[0].set_ylabel("front speed (m/s)")
        axes[1].plot(lam, [r["lambda2_v_m3_s"] * 1e17 for r in spacing], "o-", color="#1b3a5c", ms=4, lw=1.1,
                     label="PINNACLE")
        axes[1].axhline(0.41, color="#c8553d", ls=":", lw=1.0, label="simulation, Rátkai et al. 2017")
        axes[1].axhline(1.4, color="0.4", ls="--", lw=1.0, label="experiment (cited there)")
        axes[1].set_ylabel(r"$\lambda^2 v$ ($10^{-17}$ m$^3$/s)")
        axes[1].set_yscale("log")
        axes[1].legend(frameon=False, fontsize=7)
        axes[2].plot(lam, [r["beta_fraction"] for r in spacing], "o-", color="#3aa39b", ms=4, lw=1.1)
        axes[2].set_ylabel("Cu-rich fraction of the grown solid")
        for ax in axes:
            ax.set_xlabel(r"lamellar spacing $\lambda$ (nm)")
        fig.tight_layout()
        fig.savefig(FIG / "f2_spacing_measures.png", bbox_inches="tight")
        plt.close(fig)
    if len(cooling) > 1:
        panel(cooling, lambda r: rf"$\Delta T$ = {r['undercooling']:.0f} K", "f2_undercooling_maps")
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
        dt = [r["undercooling"] for r in cooling]
        axes[0].plot(dt, [r["speed_m_s"] for r in cooling], "o-", color="#1b3a5c", ms=4, lw=1.1)
        axes[0].set_ylabel("front speed (m/s)")
        axes[1].plot(dt, [r["c_alpha_grown"] for r in cooling], "o-", color="#1b3a5c", ms=4, lw=1.1, label="Ag-rich, grown")
        axes[1].plot(dt, [r["c_alpha_eq"] for r in cooling], "--", color="#1b3a5c", lw=1.0, label="Ag-rich, equilibrium")
        axes[1].plot(dt, [1 - r["c_beta_grown"] for r in cooling], "s-", color="#c8553d", ms=4, lw=1.1,
                     label="Cu-rich, grown ($1-c$)")
        axes[1].plot(dt, [1 - r["c_beta_eq"] for r in cooling], "--", color="#c8553d", lw=1.0,
                     label="Cu-rich, equilibrium ($1-c$)")
        axes[1].set_ylabel("minority fraction in each solid")
        axes[1].legend(frameon=False, fontsize=7)
        for ax in axes:
            ax.set_xlabel(r"undercooling $\Delta T$ (K)")
        fig.tight_layout()
        fig.savefig(FIG / "f2_undercooling_measures.png", bbox_inches="tight")
        plt.close(fig)
    # animation of the reference case
    ref = RES / "f2_lamellar_T955.5_c0.345_lam16_fem"
    files = series_files(ref) if ref.exists() else []
    if files:
        frames = [grid(path) for _, path in files]
        xs, ys = frames[0][0], frames[0][1]
        fig, ax = plt.subplots(figsize=(3.2, 3.4))
        im = ax.imshow(np.tile(frames[0][3][:, :-1], (1, 3)), origin="lower", extent=[0, 3 * xs[-1], 0, ys[-1]],
                       cmap="viridis", vmin=0, vmax=1)
        ax.set_xlabel("nm"); ax.set_ylabel("nm")
        title = ax.set_title("")

        def draw(i):
            im.set_data(np.tile(frames[i][3][:, :-1], (1, 3)))
            title.set_text(f"frame {i}")
            return im, title

        FuncAnimation(fig, draw, frames=len(frames)).save(FIG / "f2_lamellar_T955.5_lam16.gif",
                                                          writer=PillowWriter(fps=4), dpi=80)
        plt.close(fig)
    # long run: is the lamellar growth steady?
    long = RES / "f2_long100_T955.5_c0.345_lam16_fem"
    if (long / "solver_history.csv").exists():
        hist = list(csv.DictReader((long / "solver_history.csv").open()))
        t = np.array([float(r["t"]) for r in hist])
        xs, ys, phi, c = grid(long / "fields.vtk")
        front = (1.0 - np.array([float(r["phi_mean"]) for r in hist])) * ys[-1]
        solid_rows = (phi.mean(axis=1) < 0.05) & (ys > SEED_HEIGHT + 2.0)
        beta_y = (c[:, :-1] > 0.5).mean(axis=1)
        half = len(t) // 2
        v = float(np.polyfit(t[half:], front[half:], 1)[0])
        v_first = float(np.polyfit(t[len(t) // 10: half], front[len(t) // 10: half], 1)[0])
        gr = groups(955.5)
        lever = (0.345 - gr["c_alpha"]) / (gr["c_beta"] - gr["c_alpha"])
        print(f"long run: front {front[0]:.1f} -> {front[-1]:.1f} nm, speed {v:.5f} nm per unit = "
              f"{v * 1e-9 / TIME_UNIT:.4f} m/s (first half {v_first * 1e-9 / TIME_UNIT:.4f}), Cu-rich fraction "
              f"{beta_y[solid_rows].mean():.3f} +/- {beta_y[solid_rows].std():.3f} (lever rule {lever:.3f}), "
              f"lambda^2 v = {(16e-9) ** 2 * v * 1e-9 / TIME_UNIT:.3e} m^3/s")
        with (RES / "long_summary.csv").open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["front_start_nm", "front_end_nm", "speed_m_s", "speed_first_half_m_s", "beta_fraction",
                         "beta_fraction_std", "lever_rule", "lambda2_v_m3_s"])
            wr.writerow([front[0], front[-1], v * 1e-9 / TIME_UNIT, v_first * 1e-9 / TIME_UNIT,
                         beta_y[solid_rows].mean(), beta_y[solid_rows].std(), lever,
                         (16e-9) ** 2 * v * 1e-9 / TIME_UNIT])
        fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), gridspec_kw={"width_ratios": [1.0, 1.2, 1.2]})
        shown = np.where(phi > 0.5, np.nan, c)
        im = axes[0].imshow(np.tile(shown[:, :-1], (1, 3)), origin="lower", extent=[0, 3 * xs[-1], 0, ys[-1]],
                            cmap="viridis", vmin=0, vmax=1, aspect="equal")
        axes[0].set_xlabel("nm"); axes[0].set_ylabel("growth direction (nm)")
        fig.colorbar(im, ax=axes[0], shrink=0.8, label="Cu fraction $c$")
        axes[1].plot(t * TIME_UNIT * 1e6, front, color="#1b3a5c", lw=1.3)
        axes[1].set_xlabel(r"time ($\mu$s)"); axes[1].set_ylabel("front position (nm)")
        axes[2].plot(ys[solid_rows], beta_y[solid_rows], color="#3aa39b", lw=1.2)
        axes[2].axhline(lever, color="0.4", ls="--", lw=1.0, label="lever rule")
        axes[2].set_ylim(0, 0.5); axes[2].set_xlabel("height in the grown solid (nm)")
        axes[2].set_ylabel("Cu-rich fraction across the box"); axes[2].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG / "f2_long.png", bbox_inches="tight")
        plt.close(fig)
    print("wrote figures/crystallization/f2_*.png")


if __name__ == "__main__":
    main()
