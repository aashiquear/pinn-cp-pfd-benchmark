#!/usr/bin/env python3
"""Case B1, network routes against the FE saddle: table and figures.

Usage: python3 scripts/cryst_wall_network_figures.py
Reads results/crystallization/walls/<run>/ of the decks written by
scripts/cryst_wall_network_decks.py and of the FE decks they derive from; writes
results/crystallization/walls/network_table.csv,
figures/crystallization/b1_wall_network.png (40 K, both models: contour lines of the FE
saddle, of the starting cap and of the trained network; error map; training history) and
figures/crystallization/b1_wall_newton_control.png (Model II: where plain Newton ends from
the cap, with and without the network).

Columns of the table: `energy` is free_energy (+ wall_energy for Model I) of the half domain;
`d_energy_rel` its relative difference from the FE saddle deck; `field_max`, `field_rms` the
nodal difference from the FE saddle field; `pred_max` (hybrid) the distance of the network
prediction from the corrected field; `newton` the corrector's (or the control's) iterations;
`wall_s` the run time recorded by the run list.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/walls"
FIG = ROOT / "figures/crystallization"
DELTA_NM = 2.0 / (2.0 * np.log(9.0))
LEVELS = [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6]


def energy(run: str) -> float:
    row = next(csv.DictReader((RES / run / "postprocessors.csv").open()))
    return float(row["free_energy"]) + float(row.get("wall_energy", 0.0) or 0.0)


def nodal(run: str, key: str = "phi"):
    m = read_vtk(RES / run / "fields.vtk")
    return m["nodes"], np.asarray(m["point"][key], dtype=float)


def grid(nodes: np.ndarray, u: np.ndarray):
    x, y = nodes[:, 0], nodes[:, 1]
    xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
    g = np.zeros((len(ys), len(xs)))
    g[np.searchsorted(ys, np.round(y, 9)), np.searchsorted(xs, np.round(x, 9))] = u
    return xs * DELTA_NM, ys * DELTA_NM, g


def cap(run_stem: str, nodes: np.ndarray) -> np.ndarray:
    """The initial condition of the deck (its `cap` function) on the nodes."""
    text = (ROOT / f"decks/crystallization/walls/{run_stem}_fem.i").read_text()
    expr = re.search(r"\[\./cap\]\s*\n\s*kind\s*=\s*parsed\s*\n\s*expression\s*=\s*(\S+)", text).group(1)
    return eval(expr, {"x": nodes[:, 0], "y": nodes[:, 1], "tanh": np.tanh, "sqrt": np.sqrt, "exp": np.exp})


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    walls = {}
    for log in (RES / "runs_network.txt", RES / "runs_network_20k.txt"):
        for ln in log.read_text().splitlines() if log.exists() else []:
            m = re.match(r"(\S+) exit (\d+) wall (\d+) s", ln)
            if m:
                walls[m.group(1)] = (int(m.group(2)), int(m.group(3)))
    rows = []
    stems = sorted({re.sub(r"_(pinn_seed\d|hybrid\d*_seed\d|newton)$", "", p.name)
                    for p in RES.glob("b1_wall_model*") if re.search(r"_(pinn_seed\d|hybrid\d*_seed\d|newton)$", p.name)})
    for stem in stems:
        fe = f"{stem}_fem"
        e_fe = energy(fe)
        nodes, u_fe = nodal(fe)
        guess = cap(stem, nodes)
        model = "I" if "model1" in stem else "II"
        d_t = int(re.search(r"dT(\d+)", stem).group(1))
        inv = json.load((RES / fe / "summary.json").open())["phase_field"]["invariants"]
        rows.append({"model": model, "dT_K": d_t, "route": "FE saddle search", "seed": "", "exit": 0,
                     "energy": e_fe, "d_energy_rel": 0.0, "field_max": 0.0, "field_rms": 0.0, "pred_max": "",
                     "newton": inv.get("max_newton_iters", ""), "wall_s": ""})
        rows.append({"model": model, "dT_K": d_t, "route": "starting cap", "seed": "", "exit": "",
                     "energy": "", "d_energy_rel": "", "field_max": float(np.abs(guess - u_fe).max()),
                     "field_rms": float(np.sqrt(np.mean((guess - u_fe) ** 2))), "pred_max": "", "newton": "", "wall_s": ""})
        for run in sorted(p.name for p in RES.glob(f"{stem}_*") if re.search(r"_(pinn_seed\d|hybrid\d*_seed\d|newton)$", p.name)):
            route = ("PINN" if "_pinn_" in run else "PINN+FEM" if "_hybrid_" in run
                     else "PINN+FEM (150 iterations)" if "_hybrid150_" in run else "FE plain Newton")
            seed = (re.search(r"seed(\d)", run) or [None, ""])[1]
            code, wall = walls.get(run, ("", ""))
            if not (RES / run / "postprocessors.csv").exists():
                rows.append({"model": model, "dT_K": d_t, "route": route, "seed": seed, "exit": code, "energy": "",
                             "d_energy_rel": "", "field_max": "", "field_rms": "", "pred_max": "", "newton": "", "wall_s": wall})
                continue
            e = energy(run)
            _, u = nodal(run)
            inv = json.load((RES / run / "summary.json").open())["phase_field"]["invariants"]
            rows.append({"model": model, "dT_K": d_t, "route": route, "seed": seed, "exit": code, "energy": e,
                         "d_energy_rel": (e - e_fe) / abs(e_fe), "field_max": float(np.abs(u - u_fe).max()),
                         "field_rms": float(np.sqrt(np.mean((u - u_fe) ** 2))),
                         "pred_max": inv.get("hybrid_prediction_max_abs", ""),
                         "newton": inv.get("hybrid_corrector_newton_iters_first_step", inv.get("max_newton_iters", "")),
                         "wall_s": wall})
    with (RES / "network_table.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))

    # ---- figure 1: the trained network against the FE saddle at 40 K ----
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 5.4))
    for i, (label, stem) in enumerate((("Model I", "b1_wall_model1_dT040"), ("Model II", "b1_wall_model2_dT040_nucleus"))):
        nodes, u_fe = nodal(f"{stem}_fem")
        x, y, g_fe = grid(nodes, u_fe)
        _, _, g_cap = grid(nodes, cap(stem, nodes))
        _, u_net = nodal(f"{stem}_pinn_seed0")
        _, _, g_net = grid(nodes, u_net)
        ax = axes[i, 0]
        ax.contour(x, y, g_cap, levels=LEVELS, colors="0.6", linewidths=0.7, linestyles="dotted")
        ax.contour(x, y, g_fe, levels=LEVELS, colors="black", linewidths=0.9)
        ax.contour(x, y, g_net, levels=LEVELS, colors="#c8553d", linewidths=0.9, linestyles="dashed")
        ax.set_xlim(0, 9); ax.set_ylim(0, 6); ax.set_aspect("equal")
        ax.set_ylabel(r"$y$ (nm)"); ax.set_title(f"{label}, 40 K: contour lines", fontsize=9)
        ax.plot([], [], color="0.6", ls="dotted", label="starting cap")
        ax.plot([], [], color="black", label="FE saddle")
        ax.plot([], [], color="#c8553d", ls="dashed", label="PINN")
        ax.legend(fontsize=7, loc="upper right", frameon=False)
        ax = axes[i, 1]
        im = ax.pcolormesh(x, y, g_net - g_fe, cmap="RdBu_r", shading="gouraud", rasterized=True,
                           vmin=-np.abs(g_net - g_fe).max(), vmax=np.abs(g_net - g_fe).max())
        ax.set_xlim(0, 9); ax.set_ylim(0, 6); ax.set_aspect("equal")
        ax.set_title(r"$\phi_{PINN} - \phi_{FE}$", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
        ax = axes[i, 2]
        for seed, color in zip((0, 1, 2), ("#1f4e79", "#c8553d", "#3a7d44")):
            run = RES / f"{stem}_pinn_seed{seed}" / "solver_history.csv"
            if not run.exists():
                continue
            h = np.genfromtxt(run, delimiter=",", names=True)
            ax.semilogy(np.arange(1, len(h) + 1), h["loss"], color=color, lw=1.0, label=f"seed {seed}")
        ax.set_title("training history", fontsize=9); ax.set_ylabel("loss (relative to the cap's residual)")
        ax.legend(fontsize=7, frameon=False)
        if i == 1:
            axes[i, 0].set_xlabel(r"$x$ (nm)"); axes[i, 1].set_xlabel(r"$x$ (nm)")
            ax.set_xlabel("history row (regression, then L-BFGS)")
    fig.tight_layout()
    fig.savefig(FIG / "b1_wall_network.png", bbox_inches="tight")
    plt.close(fig)

    # ---- figure 2: Model II, plain Newton from the cap with and without the network ----
    shown = [d for d in (20, 40, 60, 80) if (RES / f"b1_wall_model2_dT{d:03d}_nucleus_newton/fields.vtk").exists()]
    if shown:
        fig, axes = plt.subplots(2, len(shown), figsize=(3.3 * len(shown), 4.6), squeeze=False)
        for j, d_t in enumerate(shown):
            stem = f"b1_wall_model2_dT{d_t:03d}_nucleus"
            nodes, u_fe = nodal(f"{stem}_fem")
            x, y, g_fe = grid(nodes, u_fe)
            for i, (run, label) in enumerate(((f"{stem}_newton", "plain Newton from the cap"),
                                              (f"{stem}_hybrid_seed0", "Newton from the network"))):
                ax = axes[i, j]
                if (RES / run / "fields.vtk").exists():
                    _, u = nodal(run)
                    _, _, g = grid(nodes, u)
                    ax.pcolormesh(x, y, g, cmap="Greys", vmin=0, vmax=1, shading="gouraud", rasterized=True)
                    ax.contour(x, y, g_fe, levels=[0.5], colors="#c8553d", linewidths=0.9)
                else:
                    ax.text(0.5, 0.5, "no converged state", transform=ax.transAxes, ha="center", fontsize=8)
                ax.set_aspect("equal"); ax.set_title(f"{d_t} K: {label}", fontsize=8)
                if j == 0:
                    ax.set_ylabel(r"$y$ (nm)")
                if i == 1:
                    ax.set_xlabel(r"$x$ (nm)")
        fig.tight_layout()
        fig.savefig(FIG / "b1_wall_newton_control.png", bbox_inches="tight")
        plt.close(fig)
    print("wrote figures/crystallization/b1_wall_network.png, b1_wall_newton_control.png")


if __name__ == "__main__":
    main()
