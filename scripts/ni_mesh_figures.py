#!/usr/bin/env python3
"""Figures of the Ni mesh study (Part III, docs/12 §3), driven by one deck.

Reads ``decks/ni/figures/mesh_figures.i`` and writes, in the style of the
calibration figures (scripts/ni_calibration_figures.py):

  mesh_curves.png       every mesh's curve against the Fig. 3 experiment markers
  mesh_convergence.png  engineering stress at fixed gauge strains against the gage
                        cell size (one-layer meshes joined; two-layer mesh marked)

Usage: ni_mesh_figures.py decks/ni/figures/mesh_figures.i
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ni_calibration_figures import box, exp_marker, read_curve, read_markers  # noqa: E402
from report_figures import _resolve, _str  # noqa: E402

MESH_COLORS = {"h040": "#1f4e9c", "h020": "#d1495b", "h010": "#00798c", "h020z2": "#edae49"}


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("figures", "mesh"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    f, meshes = cfg["figures"], cfg["mesh"]
    import ni_calibration_figures as ncf

    ncf.CUT_AT_CUTBACK = str(f.get("cut_at_first_cutback", "true")).lower() in ("true", "1", "yes")
    out = _resolve(_str(f["out_dir"]), deck)
    out.mkdir(parents=True, exist_ok=True)
    ex, ey = read_markers(_resolve(_str(f["markers"]), deck))
    curves = {n: read_curve(_resolve(_str(m["results"]), deck)) for n, m in meshes.items()}

    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=int(f["dpi"]))
    for n, m in meshes.items():
        e, s, _cut = curves[n]
        ls = "--" if int(m["layers"]) > 1 else "-"
        ax.plot(e, s, ls, color=MESH_COLORS.get(n, "#555"), lw=1.7,
                label="%s, gage cell %.1f mm" % (_str(m["label"]), float(m["h"])))
    exp_marker(ax, ex, ey)
    box(ax, f, max([ey.max()] + [c[1].max() for c in curves.values()]))
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(out / "mesh_curves.png")
    plt.close(fig)

    strains = [float(x) for x in _str(f["strains"]).split()]
    fig, ax = plt.subplots(figsize=(5.6, 4.0), dpi=int(f["dpi"]))
    markers = ("o", "s", "^", "D")
    for k, eps in enumerate(strains):
        hs, ss = [], []
        for n, m in meshes.items():
            e, s, _cut = curves[n]
            if eps > e[-1]:
                continue
            val = float(np.interp(eps, e, s))
            if int(m["layers"]) > 1:
                ax.plot(float(m["h"]), val, markers[k % 4], mfc="white", color="#555", ms=7)
            else:
                hs.append(float(m["h"]))
                ss.append(val)
        order = np.argsort(hs)
        ax.plot(np.array(hs)[order], np.array(ss)[order], "-" + markers[k % 4], lw=1.4, ms=5,
                label=r"$\varepsilon$ = %g%%" % eps)
    ax.plot([], [], "o", mfc="white", color="#555", label="h020, two layers (open)")
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.1, 0.2, 0.4])
    ax.set_xticklabels(["0.1", "0.2", "0.4"])
    ax.set_xlim(0.08, 0.5)
    ax.set_xlabel("gage cell size $h$ (mm)")
    ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)")
    ax.grid(False)
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out / "mesh_convergence.png")
    plt.close(fig)
    print("wrote", sorted(p.name for p in out.glob("*.png")))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
