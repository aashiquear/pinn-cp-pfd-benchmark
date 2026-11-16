#!/usr/bin/env python3
"""Stress-strain curves of Ni CP+PFD runs against the experiment (docs/12 §6), driven by one deck.

  [Figures]
    out            = ../../../figures/ni/pfd/curves_L1_cal.png
    markers        = ../../../data/digitized/ni_fig3_experiment_markers.csv
    experiment     = ../../../data/digitized/ni_fig3_specimenL_calibration.csv  # full Fig. 3 trace
    cp_only        = ../../../results/ni_study/cal_h020/gn2/base/results.csv   # calibrated CP, no damage
    failure_target = 0.463          # experiment failure strain, Fig. 4(b): the experiment ends here
    damage         = true           # second row: max d against gauge strain
    d_markers      = true           # markers where a d > 0.55 .. 0.85 zone first cuts the section (thresholds.json)
    strain_max     = 60
    stress_max     = 300
    stress_tick    = 50
    legend         = inside         # inside (lower right) | below (under each panel)
    dpi            = 200
  []
  [Runs]
    [./l06]
      dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.6/gc30.95   # an ni_specimen.py CP+PFD run
      label = '$\\ell_0$ = 0.6 mm'
      panel = 0                      # optional: runs with the same panel share an axis (default 0)
    []
  []

Each run is drawn over every step with a converged mechanics solve (results.csv, eng_stress
finite), so the load drop of the failing step shows; its failure (summary.json, the first
step whose d > 0.85 zone cuts the grips apart) is marked with a cross. A run that stopped
without separating ends with an open square. Style as the CP calibration figures
(scripts/ni_calibration_figures.py): curves from (0, 0), closed axes box, no grid.

Usage: ni_pfd_curve_figures.py FIGURES.i
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ni_calibration_figures as ncf  # noqa: E402
from ni_calibration_figures import box, exp_marker, read_curve, read_markers  # noqa: E402
from report_figures import _resolve, _str  # noqa: E402

D_MARKERS = {0.55: "o", 0.65: "s", 0.75: "^", 0.85: "D"}
RUN_COLORS = ("#d1495b", "#edae49", "#00798c", "#1f4e9c", "#66a182", "#8d5a97", "#9a9a9a")


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "figures" not in cfg:
        raise SystemExit(f"{path}: missing [Figures] block")
    if "runs" not in cfg and "curves" not in cfg:
        raise SystemExit(f"{path}: missing [Runs] or [Curves] block")
    cfg.setdefault("runs", {})
    return cfg


def read_run(run: Path) -> dict:
    """Gauge strain (%), engineering stress and max d over the converged steps, with (0, 0, 0)."""
    rows = [r for r in csv.DictReader(open(run / "results.csv")) if np.isfinite(float(r["eng_stress"]))]
    e = np.array([0.0] + [100.0 * float(r["gauge_strain"]) for r in rows])
    s = np.array([0.0] + [float(r["eng_stress"]) for r in rows])
    d = np.array([0.0] + [float(r["max_damage"]) for r in rows])
    summ = json.loads((run / "summary.json").read_text())
    gc = json.loads((run / "case.json").read_text())["substitutions"]["GC"]
    fs = summ.get("failure_gauge_strain")
    marks = run / "thresholds.json"
    th = json.loads(marks.read_text())["markers"] if marks.exists() else []
    return {"e": e, "s": s, "d": d, "Gc": gc, "fail": None if fs is None else 100.0 * fs, "th": th}


def read_experiment(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [r for r in csv.DictReader(ln for ln in open(path) if not ln.startswith("#"))
            if r["series"] == "experiment_specimenL" and float(r["strain"]) >= 0.0]
    return (np.array([100.0 * float(r["strain"]) for r in rows]),
            np.array([float(r["stress_MPa"]) for r in rows]))


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    f, runs = cfg["figures"], cfg["runs"]
    # damage-free curves are cut at their first load cutback only on request (runs made before
    # PINNACLE ddc76c3 doubled the strain rate after a cutback; later runs are drawn whole)
    ncf.CUT_AT_CUTBACK = _str(f.get("cut_at_first_cutback", "false")).lower() in ("true", "1", "yes")
    panels = sorted({int(r.get("panel", 0)) for r in runs.values()}) or [0]
    with_d = _str(f.get("damage", "true")).lower() == "true"
    nrow = 2 if with_d else 1
    fig, axes = plt.subplots(nrow, len(panels), figsize=(4.6 * len(panels), 3.6 + 2.0 * with_d),
                             squeeze=False, height_ratios=[1.8, 1.0][:nrow],
                             dpi=int(f.get("dpi", 200)), sharex="col")
    ex, ey = read_markers(_resolve(_str(f["markers"]), deck))
    tx, ty = read_experiment(_resolve(_str(f["experiment"]), deck))
    # cp_only is optional; a [Curves] block draws runs without damage (results.csv, label, color)
    has_cp = "cp_only" in f
    ce, cs = (read_curve(_resolve(_str(f["cp_only"]), deck))[:2] if has_cp
              else (np.zeros(1), np.zeros(1)))
    curves = cfg.get("curves", {})
    target = 100.0 * float(f["failure_target"])
    # the experiment is shown only up to the adopted failure strain (author 2026-10-20: the
    # Fig. 3 trace runs on to 51% and its 50% marker confuses the comparison)
    ex, ey = ex[ex <= target + 1e-9], ey[ex <= target + 1e-9]
    ty, tx = ty[tx <= target], tx[tx <= target]
    ymax = max(ty.max(), cs.max())
    data = {k: read_run(_resolve(_str(r["dir"]), deck)) for k, r in runs.items()}
    plain = {k: read_curve(_resolve(_str(c["path"]), deck))[:2] for k, c in curves.items()}
    for j, p in enumerate(panels):
        ax = axes[0, j]
        if _str(f.get("hide_experiment", "false")).lower() != "true":   # other specimens: no digitized curve
            ax.plot(tx, ty, "-", color="#111111", lw=0.7, alpha=0.6, label="Experiment, Fig. 3 trace to the failure strain")
            exp_marker(ax, ex, ey)
            if has_cp:
                ax.plot(ce, cs, "--", color="#9a9a9a", lw=1.3, label="CP only (calibrated, no damage)")
        for i, (k, (pe, ps)) in enumerate(plain.items()):
            ax.plot(pe, ps, "-", color=_str(curves[k].get("color", RUN_COLORS[i % len(RUN_COLORS)])),
                    lw=1.5, label=_str(curves[k].get("label", k)))
            ymax = max(ymax, ps.max())
        ax.axvline(target, color="#111111", lw=0.9, ls=":",
                   label="experiment failure strain %.1f%% (Fig. 4(b))" % target)
        keys = [k for k, r in runs.items() if int(r.get("panel", 0)) == p]
        for i, k in enumerate(keys):
            r, c = data[k], _str(runs[k].get("color", RUN_COLORS[i % len(RUN_COLORS)]))
            end = ("failure at %.1f%%" % r["fail"]) if r["fail"] is not None else (
                "stopped unseparated at %.1f%%" % r["e"][-1])
            head = _str(runs[k].get("label", ""))
            ax.plot(r["e"], r["s"], "-", color=c, lw=1.5,
                    label=r"%s$G_c$ = %s MPa mm: %s" % (head + ", " if head else "", r["Gc"], end))
            ax.plot(r["e"][-1], r["s"][-1], "x" if r["fail"] is not None else "s", color=c,
                    ms=7 if r["fail"] is not None else 5, mew=1.6, mfc="none", zorder=6)
            ymax = max(ymax, r["s"].max())
            if _str(f.get("d_markers", "false")).lower() == "true":
                for m in r["th"]:
                    ax.plot(100.0 * m["gauge_strain"], m["eng_stress"], D_MARKERS.get(m["d"], "o"), color=c,
                            ms=5.5, mfc="white", mew=1.2, zorder=7, ls="none",
                            label=r"$d > %.2f$ across the section at %.1f%%" % (m["d"], 100.0 * m["gauge_strain"]))
            if with_d:
                axes[1, j].plot(r["e"], r["d"], "-", color=c, lw=1.5)
                axes[1, j].plot(r["e"][-1], r["d"][-1], "x" if r["fail"] is not None else "s",
                                color=c, ms=7 if r["fail"] is not None else 5, mew=1.6, mfc="none")
        box(ax, f, ymax, xlabel=not with_d, ylabel=(j == 0))
        if "panel_title" in f:
            ax.set_title(_str(f["panel_title"]).split("|")[j].strip(), fontsize=9)
        if _str(f.get("legend", "inside")) == "below":
            axes[-1, j].legend(*ax.get_legend_handles_labels(), loc="upper center", frameon=False,
                               fontsize=6.5, bbox_to_anchor=(0.5, -0.32))
        else:
            ax.legend(loc="lower right", frameon=False, fontsize=6.5)
        if with_d:
            a = axes[1, j]
            a.axvline(target, color="#111111", lw=0.9, ls=":")
            a.axhline(0.85, color="#777777", lw=0.8, ls="--")
            a.text(1.0, 0.86, "d = 0.85 (failure, across the section)", fontsize=6.5, va="bottom")
            a.set_xlim(0.0, float(f["strain_max"]))
            a.set_ylim(0.0, 1.05)
            a.set_xlabel(r"Engineering (gauge) strain $\varepsilon$ (%)")
            if j == 0:
                a.set_ylabel(r"max phase field $d$")
    out = _resolve(_str(f["out"]), deck)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
