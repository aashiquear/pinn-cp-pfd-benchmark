#!/usr/bin/env python3
"""Parameter sensitivity of a neural route against the FE sensitivity study.

Deck-driven (PINNACLE ``.i`` syntax):

  [Compare]
    fe_root   = ../../results/cp_study/sens/h40           # FE cases (sens_h40.i)
    fe_base   = ../../results/cp_study/cal/gn2/base        # FE base case
    out_dir   = ../../results/cp_study/sens_neural
    figures   = ../../figures/report
    prefix    = sens_neural_h40
    markers   = ../../data/digitized/al2219_fig3_experiment_markers.csv
  [Routes]
    [./linear]  root = ../../results/cp_study/sens_neural/h40_linear       base = ../../results/cp_study/hybrid/h40_linear          label = 'PINN+FEM (linear extrapolation)' []
    [./grain]   root = ../../results/cp_study/sens_neural/h40_grain_solver base = ../../results/cp_study/hybrid/h40_fem_solver_grain label = 'PINN as solver, grain-aware (verified)' []
  []

For every case of the FE study (``<fe_root>/<case>/results.csv``, plus the
base) and every route (``<root>/<case>/`` with ``results.csv`` and
``solver/solver_history.csv``), the script reports:

* the maximum and the 20%-strain deviation of the route's engineering
  stress from the FE case (the accuracy under the parameter change);
* Newton iterations, steps solved without Newton, predictor time and total
  solve time of the route, and the FE case's iterations (the cost);
* the FE sensitivity itself (stress at 20% against the base) so the
  route's deviation can be read against the effect it is tracking.

Writes ``<out_dir>/<prefix>.csv`` and two figures: the route-vs-FE
deviation per case, and the cost per case.

Usage: neural_sensitivity.py decks/report/neural_sensitivity.i
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("compare", "routes"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def read_curve(d: Path):
    p = d / "results.csv"
    if not p.exists():
        return None
    rows = list(csv.DictReader(p.open()))
    return np.array([float(r["nominal_strain"]) for r in rows]), np.array([float(r["eng_stress"]) for r in rows])


def read_history(d: Path) -> list[dict]:
    p = d / "solver" / "solver_history.csv"
    if not p.exists():
        return []
    out = []
    for r in csv.DictReader(p.open()):
        conv = str(r.get("converged", "1")).strip().lower()
        ok = conv in ("true",) or (conv not in ("", "false") and float(conv) == 1.0)
        if ok:
            out.append({k: (float(v) if v not in ("", None, "nan") else np.nan) for k, v in r.items()})
    return out


def cost(h: list[dict]) -> dict:
    it = np.array([r.get("iterations", np.nan) for r in h])
    nt = np.nan_to_num(np.array([r.get("newton_time_s", 0.0) for r in h]))
    pt = np.nan_to_num(np.array([r.get("predictor_time_s", 0.0) for r in h]))
    return {"steps": len(h), "iterations": int(np.nansum(it)), "accepted_without_newton": int(np.sum(it == 0)),
            "newton_time_s": float(nt.sum()), "predictor_time_s": float(pt.sum()), "total_time_s": float(nt.sum() + pt.sum())}


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck); comp = cfg["compare"]
    fe_root = _resolve(str(comp["fe_root"]), deck); fe_base = _resolve(str(comp["fe_base"]), deck)
    out_dir = _resolve(str(comp["out_dir"]), deck); out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = _resolve(str(comp.get("figures", comp["out_dir"])), deck); fig_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(comp.get("prefix", "sens_neural"))
    routes = {n: (_resolve(str(r["root"]), deck), _resolve(str(r["base"]), deck), _str(r.get("label", n))) for n, r in cfg["routes"].items()}
    cases = sorted(p.name for p in fe_root.iterdir() if (p / "results.csv").exists() and p.name != "base")
    fe_curves = {"base": read_curve(fe_base), **{c: read_curve(fe_root / c) for c in cases}}
    fe_hist = {"base": read_history(fe_base), **{c: read_history(fe_root / c) for c in cases}}
    rows = []
    for c in ["base", *cases]:
        e_fe, s_fe = fe_curves[c]
        fe_20 = float(s_fe[-1]); base_20 = float(fe_curves["base"][1][-1])
        params = {}
        if c != "base" and (fe_root / c / "case.json").exists():
            params = json.loads((fe_root / c / "case.json").read_text()).get("params", {})
        for n, (root, base, label) in routes.items():
            d = base if c == "base" else root / c
            cur = read_curve(d); h = read_history(d)
            if cur is None or not h:
                continue
            e, s = cur
            dev = np.interp(e_fe, e, s) - s_fe
            rec = {"route": n, "case": c, "parameter": next(iter(params), "") if params else "", "value": next(iter(params.values()), np.nan) if params else np.nan,
                   "fe_stress_20": fe_20, "fe_effect_20_pct": 100.0 * (fe_20 - base_20) / base_20,
                   "route_stress_20": float(s[-1]), "max_abs_dev_MPa": float(np.max(np.abs(dev))), "dev_20_pct": 100.0 * (float(s[-1]) - fe_20) / fe_20,
                   "fe_iterations": int(np.nansum([r.get("iterations", np.nan) for r in fe_hist[c]])) if fe_hist[c] else -1, **cost(h)}
            rows.append(rec)
    if not rows:
        raise SystemExit("no route results found")
    with (out_dir / f"{prefix}.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k in ("route", "case", "fe_effect_20_pct", "max_abs_dev_MPa", "dev_20_pct", "iterations", "accepted_without_newton", "total_time_s")})
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    colors = ["#1f4e9c", "#d1495b", "#00798c", "#edae49"]
    order = ["base", *cases]
    x = np.arange(len(order))
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 6.2), dpi=200, sharex=True)
    width = 0.8 / max(len(routes), 1)
    for i, n in enumerate(routes):
        rr = {r["case"]: r for r in rows if r["route"] == n}
        axes[0].bar(x + (i - (len(routes) - 1) / 2) * width, [rr[c]["max_abs_dev_MPa"] if c in rr else np.nan for c in order], width, color=colors[i % 4], label=routes[n][2])
        axes[1].bar(x + (i - (len(routes) - 1) / 2) * width, [rr[c]["total_time_s"] if c in rr else np.nan for c in order], width, color=colors[i % 4], label=routes[n][2])
    axes[0].set_ylabel(r"max $|\sigma_{route} - \sigma_{FE}|$ (MPa)"); axes[1].set_ylabel("solve + predictor time (s)")
    axes[1].set_xticks(x); axes[1].set_xticklabels(order, rotation=60, ha="right", fontsize=7)
    for ax in axes:
        ax.grid(False); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_cases.png"); plt.close(fig)
    print("wrote", out_dir / f"{prefix}.csv", fig_dir / f"{prefix}_cases.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
