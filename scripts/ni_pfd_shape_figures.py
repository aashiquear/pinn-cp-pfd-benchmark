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
    ridge       = true                 # crack line: at each node row across the width, the x of max d in the gage
    ridge_always = true                # optional: draw the ridge on a run that never met the failure definition
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
color = damage, the phase field d itself with the outline in cyan). With ridge = true
the crack is also drawn as one line across the gage width (author 2026-10-23): for
every node row across the width, the node of that row with the largest d within the
gage, joined row to row and colored by that d (so the line carries the maximum of d
along x at each y; it is the likely fracture path when the field stays below crack).
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
    """The separated snapshot; for a run that never separated, its last snapshot (so marked)."""
    fj = run / "fracture.json"
    frac = json.loads(fj.read_text()) if fj.exists() else {"separated": False}   # a killed run has none
    summ = json.loads((run / "summary.json").read_text())
    if not frac.get("separated"):
        frames = sorted((run / "solver").glob("fields_step*.vtk"))
        if not frames:
            raise SystemExit(f"{run}: no separated snapshot and no snapshots (fracture.json)")
        frac = {**frac, "snapshot": frames[-1].name, "last_frame": True}
        return read_vtk(frames[-1]), frac, summ
    return read_vtk(run / "solver" / frac["snapshot"]), frac, summ


def ridge(v: dict, f: dict, run: Path) -> tuple[np.ndarray, np.ndarray]:
    """Deformed (x, y) of the max-d node of each row across the width, and its d.

    Rows are the undeformed y of the top-face nodes within the gage (the dog-bone mesh is a
    mapped grid: a row keeps its y along the straight gage)."""
    X0 = v["nodes"]
    x = X0 + v["point"]["displacement"]
    d = np.asarray(v["point"]["damage"], dtype=float)
    cells = np.asarray(v["cells"], dtype=int)
    top = np.unique(cells[:, 4:8])
    g0, g1 = _pair(f.get("ridge_x", f.get("gage_x", "8.7 16.7")))
    sel = top[(X0[top, 0] >= g0) & (X0[top, 0] <= g1)]
    rows = np.round(X0[sel, 1], 4)
    pts, vals = [], []
    for y in np.unique(rows):
        ids = sel[rows == y]
        # the crest of d along the row, between nodes: centroid of the nodes within 2 % of the
        # row maximum, weighted by their excess (a bare argmax jumps node to node on a flat crest)
        dm = float(d[ids].max())
        w = np.clip(d[ids] - 0.98 * dm, 0.0, None)
        w = w / w.sum() if w.sum() > 0 else np.full(len(ids), 1.0 / len(ids))
        pts.append(w @ x[ids, :2]); vals.append(dm)
    return np.array(pts), np.array(vals)


def draw_ridge(ax, v: dict, f: dict, run: Path):
    from matplotlib.collections import LineCollection

    pts, vals = ridge(v, f, run)
    # an undamaged run (mechanics breakdown before initiation) has no crack line to draw
    if len(pts) < 2 or vals.max() < float(f.get("ridge_min", 0.05)):
        return None
    # author 2026-11-15: the ridge is ONE thin line through the layer of maximum damage (the
    # x of max d at each row across the width); the d = crack outline is dashed beside it
    ax.plot(pts[:, 0], pts[:, 1], "-", color="black", lw=float(f.get("ridge_width", 0.9)), zorder=6,
            solid_capstyle="round")
    return None


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
    if f.get("crack") is not None:
        ax.tricontour(tri, d[used], levels=[float(f["crack"])], colors="cyan", linewidths=0.8, linestyles="--")
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
    if f.get("crack") is not None:
        ax.tricontour(tri, d[used], levels=[float(f["crack"])], colors="black", linewidths=0.7, linestyles="--")
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
    lc = None
    for i, (_, r) in enumerate(runs.items()):
        v, frac, summ = failure_frame(_resolve(_str(r["dir"]), deck))
        a = fig.add_subplot(gs[i + 1, 0])
        # the crack outline and the ridge line are drawn once the run has met the failure
        # definition (author 2026-10-28); an unfailed run shows its field only
        fr = {"crack": 0.85, **f}
        if summ.get("failure_gauge_strain") is None:
            fr["crack"] = None
            # the W_c = 0 baseline keeps the ridge on an unfailed run (author 2026-11-02)
            if _str(fr.get("ridge_always", "false")).lower() != "true":
                fr["ridge"] = "false"
        pc = draw_run(a, v, fr, lo, hi)
        if _str(fr.get("ridge", "false")).lower() == "true":
            lc = draw_ridge(a, v, fr, _resolve(_str(r["dir"]), deck)) or lc
        g = json.loads((_resolve(_str(r["dir"]), deck) / "case.json").read_text())["substitutions"]
        if frac.get("last_frame") and summ.get("failure_mode") == "load_drop_cut":
            title = r"%s, $G_c$ = %s MPa mm: failure (load drop) at %.1f%% gauge strain; last snapshot %.1f%%" % (
                _str(r["label"]), g["GC"], 100 * summ["failure_gauge_strain"], 100 * frac.get("snapshot_strain", summ["final_gauge_strain"]))
        elif frac.get("last_frame"):
            title = r"%s, $G_c$ = %s MPa mm: no separation; last snapshot, %.1f%% gauge strain" % (
                _str(r["label"]), g["GC"], 100 * summ["final_gauge_strain"])
        else:
            title = r"%s, $G_c$ = %s MPa mm: failure at %.1f%% gauge strain, at %.2f of the gage" % (
                _str(r["label"]), g["GC"], 100 * summ["failure_gauge_strain"], frac["location_fraction"])
        a.set_title(title, fontsize=8)
        if i == n - 1:
            a.set_xlabel("x (mm), deformed")
    cax = fig.add_subplot(gs[1:, 1])
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r"phase field $d$ (dashed: $d$ = %g; line: ridge of max $d$)" % float(f.get("crack", 0.85))
                 if _str(f.get("color", "exx")) == "damage" else r"$E_{xx}$ (dashed: $d$ = 0.85; line: ridge of max $d$)")
    if lc is not None and _str(f.get("color", "exx")) == "exx":
        cax2 = fig.add_axes([cax.get_position().x1 + 0.09, cax.get_position().y0,
                             0.018, cax.get_position().height])
        cb2 = fig.colorbar(lc, cax=cax2)
        cb2.set_label(r"crack line: max $d$ across the gage at each row")
    out = _resolve(_str(f["out"]), deck)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
