#!/usr/bin/env python3
"""Fracture shape of Ni CP+PFD runs against the experiment's DIC panel (Fig. 5).

Deck-driven (PINNACLE .i syntax):

  [Figures]
    out         = ../../../figures/ni/pfd/shape_L1.png
    panel       = ../../../data/digitized/ni_fig5_panels/L.1.png   # Fig. 5 crop (Exp. above Sim.)
    panel_label = 'L.1, experiment (DIC) and the paper's CP-PFD (Rezwan 2026, Fig. 5)'
    exx_range   = '0.23 0.98'          # color limits of E_xx, the experiment panel's scale
    x_range     = '5.0 20.4'           # mm, the gage, fillets and grip ends (undeformed x of the cells)
    crack       = 0.85                 # the d = crack outline (black on E_xx, cyan on d): the fracture zone
    color       = exx                  # exx | damage (the phase field d, 0..1, with the d = crack contour)
    dpi         = 200
  []
  [Runs]
    [./l04]
      dir   = ../../../results/ni_study/pfd/cal_L1_h020/l0.4/gc15.56   # an ni_specimen.py CP+PFD run
      label = 'l0 = 0.4 mm'
    []
  []

Each run is drawn at its failure step (the snapshot that `ni_specimen.py fracture`
found separated; fracture.json): plan view of the deformed top face, cells colored
by the Green-Lagrange axial strain E_xx = H_xx + (H_xx^2 + H_yx^2 + H_zx^2) / 2 from
the displacement gradient H the run records (the measure a DIC code such as Vic-2D
reports as e_xx), with the fracture zone outlined in black at d = crack (or, with
color = damage, the phase field d itself with the outline in cyan).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(STUDY_ROOT / "scripts"))

from vtk_legacy import read_vtk  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.labelsize": 9, "mathtext.fontset": "dejavusans"})


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    return Translator().translate(InputParser().parse(path.read_text()))


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    """A deck value as text (the parser makes a list of a quoted string with spaces)."""
    if isinstance(v, (list, tuple)):
        v = " ".join(str(x) for x in v)
    return str(v).strip().strip("'\"")


def _pair(v) -> tuple[float, float]:
    a, b = (float(x) for x in (v if isinstance(v, (list, tuple)) else _str(v).split()))
    return a, b


def failure_frame(run: Path) -> tuple[dict, dict, dict]:
    frac = json.loads((run / "fracture.json").read_text())
    summ = json.loads((run / "summary.json").read_text())
    if not frac.get("separated"):
        raise SystemExit(f"{run}: no separated snapshot (fracture.json)")
    return read_vtk(run / "solver" / frac["snapshot"]), frac, summ


def draw_damage(ax, v: dict, f: dict):
    """The phase field d on the deformed top face (nodal, Gouraud-shaded) with its d = crack contour."""
    import matplotlib.tri as mtri

    X0 = v["nodes"]
    x = X0 + v["point"]["displacement"]
    cells = np.asarray(v["cells"], dtype=int)
    top = cells[:, 4:8]
    d = np.asarray(v["point"]["damage"], dtype=float)
    x0, x1 = _pair(f["x_range"])
    keep = np.where((X0[cells, 0].mean(axis=1) >= x0) & (X0[cells, 0].mean(axis=1) <= x1))[0]
    tris = np.concatenate([top[keep][:, [0, 1, 2]], top[keep][:, [0, 2, 3]]])
    used = np.unique(tris)
    remap = -np.ones(len(x), dtype=int)
    remap[used] = np.arange(len(used))
    tri = mtri.Triangulation(x[used, 0], x[used, 1], remap[tris])
    tpc = ax.tripcolor(tri, d[used], shading="gouraud", cmap="inferno", vmin=0.0, vmax=1.0)
    ax.tricontour(tri, d[used], levels=[float(f.get("crack", 0.85))], colors="cyan", linewidths=0.8)
    ax.set_xlim(x[used, 0].min(), x[used, 0].max())
    ax.set_ylim(x[used, 1].min() - 0.2, x[used, 1].max() + 0.2)
    ax.set_aspect("equal")
    ax.set_yticks([])
    return tpc


def draw_run(ax, v: dict, f: dict, lo: float, hi: float):
    if _str(f.get("color", "exx")) == "damage":
        return draw_damage(ax, v, f)
    X0 = v["nodes"]
    x = X0 + v["point"]["displacement"]
    cells = np.asarray(v["cells"], dtype=int)
    top = cells[:, 4:8]                                   # z+ face of each HEX8
    # snapshots carry H per element (cell data); the end-of-run file per node
    src = v["cell"].get("strain")
    H = np.asarray(src if src is not None else v["point"]["strain"], dtype=float).reshape(-1, 3, 3)
    Exx = H[:, 0, 0] + 0.5 * (H[:, 0, 0] ** 2 + H[:, 1, 0] ** 2 + H[:, 2, 0] ** 2)
    if src is None:
        Exx = Exx[cells].mean(axis=1)
    d = np.asarray(v["point"]["damage"], dtype=float)
    x0, x1 = _pair(f["x_range"])
    xc = X0[cells, 0].mean(axis=1)
    keep = (xc >= x0) & (xc <= x1)
    polys = [x[top[k]][:, :2] for k in np.where(keep)[0]]
    vals = Exx[keep]
    pc = PolyCollection(polys, array=vals, cmap="jet", edgecolors="none")
    pc.set_clim(lo, hi)
    ax.add_collection(pc)
    # the fracture zone as the outline of d = crack (author 2026-10-20: a phase field is
    # diffuse, so the damaged area's outline is shown rather than a single crack path)
    import matplotlib.tri as mtri

    kk = np.where(keep)[0]
    tris = np.concatenate([top[kk][:, [0, 1, 2]], top[kk][:, [0, 2, 3]]])
    used = np.unique(tris)
    remap = -np.ones(len(x), dtype=int)
    remap[used] = np.arange(len(used))
    tri = mtri.Triangulation(x[used, 0], x[used, 1], remap[tris])
    ax.tricontour(tri, d[used], levels=[float(f.get("crack", 0.85))], colors="black", linewidths=1.0)
    ax.set_xlim(min(p[:, 0].min() for p in polys), max(p[:, 0].max() for p in polys))
    ax.set_ylim(min(p[:, 1].min() for p in polys) - 0.2, max(p[:, 1].max() for p in polys) + 0.2)
    ax.set_aspect("equal")
    ax.set_yticks([])
    return pc


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    f, runs = cfg["figures"], cfg["runs"]
    lo, hi = _pair(f["exx_range"])
    n = len(runs)
    fig = plt.figure(figsize=(7.2, 1.55 * (n + 1.4)), dpi=int(f.get("dpi", 200)))
    gs = fig.add_gridspec(n + 1, 2, width_ratios=[1.0, 0.025], height_ratios=[1.4] + [1.0] * n, hspace=0.55, wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(plt.imread(_resolve(_str(f["panel"]), deck)))
    ax.set_axis_off()
    ax.set_title(_str(f.get("panel_label", "experiment (Fig. 5)")), fontsize=8)
    pc = None
    for i, (_, r) in enumerate(runs.items()):
        v, frac, summ = failure_frame(_resolve(_str(r["dir"]), deck))
        a = fig.add_subplot(gs[i + 1, 0])
        pc = draw_run(a, v, f, lo, hi)
        g = json.loads((_resolve(_str(r["dir"]), deck) / "case.json").read_text())["substitutions"]
        a.set_title(r"%s, $G_c$ = %s MPa mm: failure at %.1f%% gauge strain, at %.2f of the gage"
                    % (_str(r["label"]), g["GC"], 100 * summ["failure_gauge_strain"], frac["location_fraction"]),
                    fontsize=8)
        if i == n - 1:
            a.set_xlabel("x (mm), deformed")
    cax = fig.add_subplot(gs[1:, 1])
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r"phase field $d$ (cyan: $d$ = %g)" % float(f.get("crack", 0.85))
                 if _str(f.get("color", "exx")) == "damage" else r"$E_{xx}$ (Green-Lagrange)")
    out = _resolve(_str(f["out"]), deck)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
