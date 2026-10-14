#!/usr/bin/env python3
"""Accuracy, time and hardware cost of the plate surrogate against the other routes.

Deck-driven (PINNACLE ``.i`` syntax):

  [Compare]
    out_dir   = ../../results/cp_study/surrogate
    figures   = ../../figures/report
    prefix    = surrogate_h40
    design    = ../../results/cp_study/surrogate/h40_data        # samples.csv (FE time per sample)
    markers   = ../../data/digitized/al2219_fig3_experiment_markers.csv
    routes    = ../../results/cp_study/hybrid/hybrid_h40_summary.csv  # cost of the solver routes (hybrid_compare.py), or a study summary (case, wall_s)
    mesh_vtk  = ../../results/pfd_study/cp/surrogate/h40_field/fields.vtk  # optional: nodes and cells for the held-out field maps
  [Surrogates]
    [./curve]  dir = ../../results/cp_study/surrogate/h40_curve  label = 'curve surrogate (GPU)' []
    [./field]  dir = ../../results/cp_study/surrogate/h40_field  label = 'field surrogate (GPU)' []
    [./curve_cpu] dir = ../../results/cp_study/surrogate/h40_curve_cpu label = 'curve surrogate (CPU)' []
  []

For every surrogate run (``surrogate_summary.json``, ``surrogate_holdout.csv``,
``holdout/<i>.csv``, ``pinnacle.log`` for the per-query time) the script
reports the held-out error, the data cost (FE time of the training
samples from the design's ``samples.csv``), the training time, the time
per query, and the break-even number of queries against every solver
route of the summary CSV (the queries after which data + training is
cheaper than one solver run per query). Figures: held-out curves (truth
against prediction) and a cost-per-query chart.

Usage: surrogate_compare.py decks/report/surrogate_compare.i
"""
from __future__ import annotations

import csv
import json
import re
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
    for key in ("compare", "surrogates"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def query_time(log: Path) -> float:
    """Seconds per query from the executioner's log line, or nan."""
    if not log.exists():
        return float("nan")
    m = re.findall(r"(\d+) queries predicted in ([\d.]+) s", log.read_text())
    if not m:
        return float("nan")
    n, t = m[-1]
    return float(t) / max(int(n), 1)


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); cfg = read_deck(deck); comp = cfg["compare"]
    out_dir = _resolve(str(comp["out_dir"]), deck); out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = _resolve(str(comp.get("figures", comp["out_dir"])), deck); fig_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(comp.get("prefix", "surrogate"))
    design = _resolve(str(comp["design"]), deck)
    samples = list(csv.DictReader((design / "samples.csv").open()))
    fe_time = {int(r["sample"]): float(r["run_time_s"]) for r in samples}
    routes = list(csv.DictReader(_resolve(str(comp["routes"]), deck).open())) if comp.get("routes") else []
    for r in routes:  # a study summary (pfd_study.py collect) names the run `case` and its cost `wall_s`
        r.setdefault("run", r.get("case", "")); r.setdefault("label", r["run"]); r.setdefault("total_solve_time_s", r.get("wall_s", "nan"))
    mesh_vtk = _resolve(str(comp["mesh_vtk"]), deck) if comp.get("mesh_vtk") else None
    rows = []
    for name, sp in cfg["surrogates"].items():
        d = _resolve(str(sp["dir"]), deck); label = _str(sp.get("label", name))
        sfile = d / "surrogate_summary.json"
        if not sfile.exists():
            print("skipping", name, "(no summary)"); continue
        s = json.loads(sfile.read_text())
        hold = list(csv.DictReader((d / "surrogate_holdout.csv").open()))
        train_ids = [int(r["sample"]) for r in hold if r["split"] == "train"]
        data_cost = float(sum(fe_time.get(i, 0.0) for i in train_ids))
        tq = query_time(d / "pinnacle.log")
        rec = {"surrogate": name, "label": label, "target": s["target"], "train": s["train"], "holdout": s["holdout"],
               "holdout_rel_l2_mean": s["holdout_rel_l2_mean"], "holdout_rel_l2_max": s["holdout_rel_l2_max"],
               "data_cost_s": data_cost, "train_time_s": s["train_time_s"], "query_s": tq, "device": s.get("device", ""),
               "network_parameters": s.get("network_parameters", 0)}
        for r in routes:
            per_run = float(r["total_solve_time_s"])
            rec[f"break_even_vs_{r['run']}"] = (data_cost + s["train_time_s"]) / max(per_run - (tq if np.isfinite(tq) else 0.0), 1e-9)
        rows.append(rec)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()})
        # held-out curves figure (reduced target)
        if s["target"] == "reduced" and (d / "holdout").exists():
            plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 7, "mathtext.fontset": "dejavusans",
                                 "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
            fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=200)
            files = sorted((d / "holdout").glob("*.csv"))
            cmap = plt.get_cmap("viridis")
            for k, f in enumerate(files):
                rr = list(csv.DictReader(f.open()))
                t = np.array([100 * float(x[next(iter(x))]) for x in rr]); yt = np.array([float(x[[c for c in x if c.endswith("_true")][0]]) for x in rr]); yp = np.array([float(x[[c for c in x if c.endswith("_pred")][0]]) for x in rr])
                c = cmap(k / max(len(files) - 1, 1))
                ax.plot(np.concatenate([[0], t]), np.concatenate([[0], yt]), color=c, lw=1.6, label=f"FE, held-out sample {f.stem}" if k < 6 else None)
                ax.plot(np.concatenate([[0], t]), np.concatenate([[0], yp]), color=c, lw=1.0, ls="--", label="surrogate" if k == 0 else None)
            top = 100.0 * (np.floor(max(float(np.max(yt)), 1.0) / 100.0) + 1)
            ax.set_xlim(0, 20.5); ax.set_ylim(0, max(700.0, top)); ax.grid(False)
            ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)"); ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)"); ax.legend(loc="lower right", frameon=False)
            fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_{name}_holdout.png"); plt.close(fig)
        # held-out field maps (field target): truth, prediction and their difference on the mesh
        if s["target"] == "field" and (d / "holdout").exists() and mesh_vtk is not None:
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            from vtk_legacy import read_vtk
            from matplotlib.collections import PolyCollection
            v = read_vtk(mesh_vtk); pts = np.asarray(v["nodes"], dtype=float); cells = np.asarray(v["cells"], dtype=int)
            quads = pts[cells[:, :4], :2]
            files = sorted((d / "holdout").glob("*.csv"))
            fig, axes = plt.subplots(len(files), 3, figsize=(9.6, 1.75 * len(files) + 0.4), dpi=200, squeeze=False)
            for k, f in enumerate(files):
                rr = list(csv.DictReader(f.open()))
                yt = np.array([float(x["c0_true"]) for x in rr]); yp = np.array([float(x["c0_pred"]) for x in rr])
                err = float(np.linalg.norm(yp - yt) / max(np.linalg.norm(yt), 1e-12))
                for j, (arr, cm, lim, ttl) in enumerate(((yt, "inferno", (0, 1), f"FE, held-out sample {f.stem}: max d {yt.max():.2f}"),
                                                         (yp, "inferno", (0, 1), f"surrogate: max d {yp.max():.2f}, rel. L2 {100 * err:.1f}%"),
                                                         (yp - yt, "coolwarm", (-0.2, 0.2), f"difference (max |.| {np.abs(yp - yt).max():.2f})"))):
                    ax = axes[k, j]
                    ax.add_collection(PolyCollection(quads, array=arr[cells[:, :4]].mean(axis=1), cmap=cm, edgecolors="none", clim=lim))
                    ax.set_aspect("equal"); ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
                    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(ttl, fontsize=7)
            fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_{name}_holdout_fields.png"); plt.close(fig)
    if rows:
        with (out_dir / f"{prefix}_summary.csv").open("w", newline="") as fh:
            keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in rows[0], k))
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) + [k for k in keys if k not in rows[0]]); w.writeheader(); w.writerows(rows)
        # cost per query: routes (one run each) vs surrogate (amortised over N queries)
        fig, ax = plt.subplots(figsize=(6.0, 3.6), dpi=200)
        n = np.logspace(0, 4, 50)
        route_colors = ["#111111", "#1f4e9c", "#d1495b", "#00798c", "#edae49", "#6a4c93", "#2a9d8f"]
        for i, r in enumerate(routes):
            ax.axhline(float(r["total_solve_time_s"]), lw=1.0, ls=":", color=route_colors[i % len(route_colors)], label=f"{r['label']} (one run)")
        sur_styles = ["-", "--", "-."]
        for i, rec in enumerate(rows):
            if np.isfinite(rec["query_s"]):
                ax.plot(n, (rec["data_cost_s"] + rec["train_time_s"]) / n + rec["query_s"], lw=1.6, ls=sur_styles[i % 3], color="#d1495b" if i == 0 else "#edae49" if i == 1 else "#6a4c93", label=rec["label"])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("queries answered"); ax.set_ylabel("cost per query (s)"); ax.grid(False); ax.legend(frameon=False, fontsize=6.5)
        fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_cost_per_query.png"); plt.close(fig)
        print("wrote", out_dir / f"{prefix}_summary.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
