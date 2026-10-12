#!/usr/bin/env python3
"""Compare the fields of neural-route runs with the FE reference on the same mesh.

Deck-driven (PINNACLE ``.i`` syntax):

  [Compare]
    figures    = ../../figures/report
    prefix     = neural_h40
    strains    = '2 8 15 20'        # nominal strains (%) of the field snapshots
    strain_max = 0.20
    vmax       = 0.6                 # colour scale of the equivalent plastic strain
    markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  [Runs]
    [./fe]      dir = ../../results/cp_study/cal/gn2/base      label = 'FE (cold Newton)' []
    [./hybrid]  dir = ../../results/cp_study/hybrid/h40_physics label = 'PINN+FEM' []
    [./pinn]    dir = ../../results/cp_study/hybrid/h40_pinn_solver label = 'PINN as solver' []
  []

The first run is the reference. For every run: ``results.csv`` (the curve)
and ``solver/fields_step*.vtk`` (the series, every snapshot the deck wrote).
Writes

* ``<prefix>_curves.png``: engineering stress-strain of every run and the
  experiment markers (report style);
* ``<prefix>_fields_<strain>pct.png``: one row per run, the element-averaged
  equivalent plastic strain on the deformed plate at that strain, plus the
  difference to the reference in the last row;
* ``<prefix>_field_errors.csv``: per run and strain, the relative L2 error
  of the displacement, the equivalent plastic strain and the Cauchy stress
  sigma_xx against the reference (same mesh, element by element).

Usage: field_compare.py decks/report/neural_compare.i
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plastic_strain_animation import read_vtk, series_files, top_quads  # noqa: E402


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("compare", "runs"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def _list(v) -> list[float]:
    return [float(x) for x in (v if isinstance(v, (list, tuple)) else str(v).split())]


def read_curve(d: Path):
    rows = list(csv.DictReader((d / "results.csv").open()))
    e = np.array([100.0 * float(r["nominal_strain"]) for r in rows]); s = np.array([float(r["eng_stress"]) for r in rows])
    return np.concatenate([[0.0], e]), np.concatenate([[0.0], s])


def read_markers(path: Path):
    lines = [ln for ln in open(path) if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    return np.array([float(r["strain_pct"]) for r in rows]), np.array([float(r["stress_MPa"]) for r in rows])


def frame_at(d: Path, strain_pct: float, strain_max: float) -> dict:
    files = series_files(d / "solver")
    n = len(files)
    j = int(np.argmin([abs(100.0 * strain_max * (k + 1) / n - strain_pct) for k in range(n)]))
    v = read_vtk(files[j][1])
    return {"nodes": v["nodes"], "cells": v["cells"], "disp": v["point"]["displacement"],
            "peeq": v["cell"]["equivalent_plastic_strain"], "stress": v["cell"].get("stress"),
            "strain_pct": 100.0 * strain_max * (j + 1) / n}


def rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1.0e-300))


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck); comp = cfg["compare"]
    fig_dir = _resolve(str(comp["figures"]), deck); fig_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(comp.get("prefix", "neural")); strains = _list(comp.get("strains", "2 8 15 20"))
    strain_max = float(comp.get("strain_max", 0.20)); vmax = float(comp.get("vmax", 0.6))
    runs = {n: (_resolve(str(r["dir"]), deck), _str(r.get("label", n))) for n, r in cfg["runs"].items()}
    runs = {n: v for n, v in runs.items() if (v[0] / "results.csv").exists()}
    names = list(runs); ref = names[0]
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "mathtext.fontset": "dejavusans",
                         "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    colors = ["#1f4e9c", "#d1495b", "#00798c", "#edae49"]
    # curves
    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=200)
    for i, n in enumerate(names):
        e, s = read_curve(runs[n][0]); ax.plot(e, s, lw=1.6 if i else 2.0, color=colors[i % 4], label=runs[n][1])
    mk = _resolve(str(comp.get("markers", "")), deck) if comp.get("markers") else None
    if mk and mk.exists():
        ex, ey = read_markers(mk); ax.plot(ex, ey, "o", ms=5.5, color="#111", mfc="white", mew=1.3, label="Experiment (Vaughan 2024, Fig. 3 of Rezwan 2025)", zorder=5)
    smax = max(float(np.max(read_curve(runs[n][0])[1])) for n in names)
    top = 100.0 * (np.floor(smax / 100.0) + 1)
    ax.set_xlim(0, 100 * strain_max + 0.5); ax.set_ylim(0, max(700.0, top)); ax.grid(False)
    ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)"); ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)")
    ax.legend(loc="lower right", frameon=False); fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_curves.png"); plt.close(fig)
    # fields
    rows_out = []
    cmap = plt.get_cmap("viridis")
    for eps in strains:
        frames = {n: frame_at(runs[n][0], eps, strain_max) for n in names}
        fr0 = frames[ref]
        nrow = len(names) + (1 if len(names) > 1 else 0)
        fig, axes = plt.subplots(nrow, 1, figsize=(6.4, 2.6 * nrow), dpi=200)
        axes = np.atleast_1d(axes)
        for i, n in enumerate(names):
            fr = frames[n]
            quads = top_quads(fr["nodes"], fr["cells"], fr["disp"], 1.0)
            pc = PolyCollection(quads, array=fr["peeq"], cmap=cmap, edgecolors="face", linewidths=0.3); pc.set_clim(0.0, vmax)
            axes[i].add_collection(pc); axes[i].set_xlim(-2, 102 + 100 * strain_max + 2); axes[i].set_ylim(-2, 52); axes[i].set_aspect("equal")
            axes[i].set_title(r"%s, $\bar\varepsilon^p$ at $\varepsilon$ = %.1f%%" % (runs[n][1], fr["strain_pct"]), fontsize=9, loc="left")
            axes[i].set_xticks([0, 50, 100]); axes[i].set_yticks([0, 25, 50])
            fig.colorbar(pc, ax=axes[i], shrink=0.9, pad=0.02)
            if n != ref:
                rows_out.append({"run": n, "strain_pct": round(fr["strain_pct"], 3),
                                 "rel_l2_displacement": rel_l2(fr["disp"], fr0["disp"]),
                                 "rel_l2_peeq": rel_l2(fr["peeq"], fr0["peeq"]),
                                 "rel_l2_sigma_xx": rel_l2(fr["stress"][:, 0], fr0["stress"][:, 0]) if fr["stress"] is not None else np.nan,
                                 "max_abs_peeq_diff": float(np.max(np.abs(fr["peeq"] - fr0["peeq"])))})
        if len(names) > 1:
            n = names[-1]; fr = frames[n]
            quads = top_quads(fr0["nodes"], fr0["cells"], fr0["disp"], 1.0)
            diff = fr["peeq"] - fr0["peeq"]; lim = max(float(np.max(np.abs(diff))), 1e-12)
            pc = PolyCollection(quads, array=diff, cmap="RdBu_r", edgecolors="face", linewidths=0.3); pc.set_clim(-lim, lim)
            axes[-1].add_collection(pc); axes[-1].set_xlim(-2, 102 + 100 * strain_max + 2); axes[-1].set_ylim(-2, 52); axes[-1].set_aspect("equal")
            axes[-1].set_title(r"%s minus %s" % (runs[n][1], runs[ref][1]), fontsize=9, loc="left")
            axes[-1].set_xticks([0, 50, 100]); axes[-1].set_yticks([0, 25, 50]); fig.colorbar(pc, ax=axes[-1], shrink=0.9, pad=0.02)
        fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_fields_{eps:g}pct.png"); plt.close(fig)
    if rows_out:
        with (fig_dir.parent.parent / "results" / "cp_study" / "hybrid" / f"{prefix}_field_errors.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_out[0])); w.writeheader(); w.writerows(rows_out)
        for r in rows_out:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", fig_dir / f"{prefix}_curves.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
