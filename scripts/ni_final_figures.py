#!/usr/bin/env python3
"""Last-stage figures of the Ni specimens: strain-to-failure parity and fracture location per route.

Deck-driven (PINNACLE .i syntax):

  [Figures]
    out_dir    = ../../../figures/ni/final
    prefix     = fem
    parity     = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv   # experiment column only
    fracture   = ../../../data/digitized/ni_fig5_fracture_descriptors.csv     # experiment, Fig. 5 panels
    label      = 'FEM'
    dpi        = 200
  []
  [Routes]                       # one or more routes; each a study root with <specimen>/summary.json
    [./fem]    root = ../../../results/ni_study/final/fem    label = 'FEM'       color = '#1f4e9c' []
  []

Writes <out_dir>/<prefix>_parity.png (experimental against predicted strain to failure, one
marker per specimen, colored by orientation group and shaped by thickness variant as the
paper's Fig. 4(b); the x = y line and the +-10% band), <out_dir>/<prefix>_location.png
(experimental against predicted fracture location along the gage, with the +-0.08 reading
uncertainty of the experiment) and <out_dir>/<prefix>_table.csv (per specimen: experimental and
predicted failure strain, error, location, mode). Only the experiment's columns of the data
files are read: the paper's own simulation values are never used (author 2026-10-25).

Usage: ni_final_figures.py FIGURES.i
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
from report_figures import _resolve, _str  # noqa: E402

GROUP_COLORS = {"9": "#1f4e9c", "10": "#2a9d8f", "H": "#8d5a97", "L": "#d1495b", "P": "#111111"}
VARIANT_MARKERS = {"1": "o", "3": "x", "4": "*"}


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("figures", "routes"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))


def experiment(f: dict, deck: Path) -> dict[str, dict]:
    exp = {r["specimen"]: {"strain": float(r["exp_failure_strain"])} for r in _rows(_resolve(_str(f["parity"]), deck))}
    for r in _rows(_resolve(_str(f["fracture"]), deck)):
        exp.setdefault(r["specimen"], {}).update({"location": float(r["location_fraction"]),
                                                  "location_unc": float(r["location_uncertainty"]), "mode": r["mode"]})
    return exp


def predicted(root: Path) -> dict[str, dict]:
    out = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir() and (p / "summary.json").exists()):
        s = json.loads((d / "summary.json").read_text())
        fr = json.loads((d / "fracture.json").read_text()) if (d / "fracture.json").exists() else {}
        spec = json.loads((d / "case.json").read_text())["substitutions"].get("SPECIMEN", d.name.upper()) if (d / "case.json").exists() else d.name.upper()
        out[str(spec)] = {"strain": s.get("failure_gauge_strain"), "mode": s.get("failure_mode"), "final": s.get("final_gauge_strain"),
                          "location": fr.get("location_fraction"), "angle": fr.get("angle_to_axis_deg"), "peak": s.get("peak_eng_stress")}
    return out


def _marker(spec: str):
    g, v = spec.split(".")
    return GROUP_COLORS.get(g, "#777777"), VARIANT_MARKERS.get(v, "o")


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    f, routes = cfg["figures"], cfg["routes"]
    out_dir = _resolve(_str(f["out_dir"]), deck)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = _str(f.get("prefix", "final"))
    exp = experiment(f, deck)
    dpi = int(f.get("dpi", 200))
    table = []
    # parity
    fig, axes = plt.subplots(1, len(routes), figsize=(4.0 * len(routes), 3.9), squeeze=False, dpi=dpi)
    figl, axl = plt.subplots(1, len(routes), figsize=(4.0 * len(routes), 3.9), squeeze=False, dpi=dpi)
    for j, (name, r) in enumerate(routes.items()):
        pred = predicted(_resolve(_str(r["root"]), deck))
        ax, al = axes[0, j], axl[0, j]
        lim = (0.3, 0.7)
        ax.plot(lim, lim, "-", color="#111111", lw=0.9)
        ax.fill_between(lim, [0.9 * x for x in lim], [1.1 * x for x in lim], color="#dddddd", lw=0, label=r"$\pm$10%")
        al.plot([0, 1], [0, 1], "-", color="#111111", lw=0.9)
        al.fill_between([0, 1], [-0.08, 0.92], [0.08, 1.08], color="#dddddd", lw=0, label=r"$\pm$0.08 (reading of Fig. 5)")
        errs = []
        for spec, e in sorted(exp.items()):
            p = pred.get(spec)
            c, m = _marker(spec)
            kw = dict(color=c, marker=m, ms=8 if m == "*" else 6, mfc="none" if m == "o" else c, mew=1.4, ls="none")
            if p and p.get("strain") is not None:
                ax.plot(e["strain"], p["strain"], **kw)
                errs.append((p["strain"] - e["strain"]) / e["strain"])
            elif p and p.get("final") is not None:
                ax.plot(e["strain"], p["final"], marker="s", color=c, mfc="none", ms=6, ls="none")  # no failure: run end
            if p and p.get("location") is not None and "location" in e:
                al.plot(e["location"], p["location"], **kw)
            table.append({"route": name, "specimen": spec, "exp_failure_strain": e["strain"], "pred_failure_strain": p.get("strain") if p else None,
                          "pred_final_strain": p.get("final") if p else None, "pred_mode": p.get("mode") if p else None,
                          "rel_error": ((p["strain"] - e["strain"]) / e["strain"]) if p and p.get("strain") is not None else None,
                          "exp_location": e.get("location"), "pred_location": p.get("location") if p else None,
                          "exp_mode": e.get("mode"), "pred_angle_deg": p.get("angle") if p else None, "pred_peak_MPa": p.get("peak") if p else None})
        for g, c in GROUP_COLORS.items():
            ax.plot([], [], "s", color=c, label=f"group {g}")
        for v, m in VARIANT_MARKERS.items():
            ax.plot([], [], m, color="#555555", mfc="none" if m == "o" else "#555555", label=f"variant x.{v}")
        ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
        ax.set_xlabel("experimental strain to failure (Fig. 4(b))"); ax.set_ylabel("predicted strain to failure" if j == 0 else "")
        rmse = float(np.sqrt(np.mean(np.square(errs)))) if errs else float("nan")
        ax.set_title("%s: %d of %d specimens, RMS error %.1f%%" % (_str(r.get("label", name)), len(errs), len(exp), 100 * rmse), fontsize=9)
        ax.legend(loc="upper left", frameon=False, fontsize=6, ncol=2)
        al.set_xlim(0, 1); al.set_ylim(0, 1); al.set_aspect("equal")
        al.set_xlabel("experimental fracture location (of the gage, Fig. 5)"); al.set_ylabel("predicted fracture location" if j == 0 else "")
        al.set_title(_str(r.get("label", name)), fontsize=9)
        al.legend(loc="upper left", frameon=False, fontsize=6)
    for fg, nm in ((fig, "parity"), (figl, "location")):
        fg.tight_layout()
        fg.savefig(out_dir / f"{prefix}_{nm}.png", bbox_inches="tight")
        print("wrote", out_dir / f"{prefix}_{nm}.png")
    with (out_dir / f"{prefix}_table.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0])); w.writeheader(); w.writerows(table)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
