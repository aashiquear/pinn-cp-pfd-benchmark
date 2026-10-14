#!/usr/bin/env python3
"""Build and render a plate mesh with a stress concentrator (hexaform ``plate``).

Deck-driven (PINNACLE ``.i`` syntax):

  [Mesh]
    out_dir    = ../../../results/pfd_study/mesh/h40_hole5
    width      = 100
    height     = 50
    thickness  = 1
    h          = 2.5                       # in-plane background cell (h40 = 100/2.5)
    hole       = circle                    # circle | ellipse | none
    diameter   = 5                         # circle: diameter; ellipse: 'major minor' (full axes)
    angle      = 0                         # ellipse: rotation of the major axis from +x, degrees
    centre     = '50 25'                   # optional, default plate centre
    ring_cells = 2                         # background cells between the hole box and the O-grid rectangle
    n_r        = auto                      # radial layers (auto: inner = hole arc, outer = h)
    grading    = auto
    points     = ../../../results/cp_study/cal/gn2/base/orientations_points.csv   # optional: grain ids for the figure
    figure     = ../../../figures/pfd/mesh_h40_hole5.png

Writes ``<out_dir>/mesh.e`` (Exodus, read by PINNACLE ``[Mesh] type = file``),
``mesh.vtk``, ``mesh.json`` (hexaform's plate summary and validation) and the
figure: the whole plate coloured by grain id (nearest point of the
orientation point file at each element centroid, which is exactly how
PINNACLE's ``point_file`` source assigns orientations) and a zoom on the
ring around the hole. Needs hexaform on the path
(``/mnt/extended/projects/hexaform/src``) and matplotlib for the figure.

Usage: plate_hole_mesh.py decks/pfd_study/mesh/h40_hole5.i
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, "/mnt/extended/projects/hexaform/src")


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "mesh" not in cfg:
        raise SystemExit(f"{path}: missing [Mesh] block")
    return cfg["mesh"]


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _words(v) -> list[str]:
    return [str(x) for x in v] if isinstance(v, (list, tuple)) else str(v).split()


def _opt(v, cast):
    return None if v is None or str(v).lower() == "auto" else cast(v)


def build(m: dict, deck: Path):
    from hexaform import HoleSpec, PlateOptions, plate_mesh, plate_summary, validate
    kind = str(m.get("hole", "circle")).lower()
    centre = tuple(float(x) for x in _words(m["centre"])) if m.get("centre") else None
    hole = None
    if kind == "circle":
        hole = HoleSpec("circle", 0.5 * float(m["diameter"]), centre=centre)
    elif kind == "ellipse":
        a, b = (float(x) for x in _words(m["diameter"]))
        hole = HoleSpec("ellipse", 0.5 * a, 0.5 * b, centre=centre, angle=float(m.get("angle", 0.0)))
    elif kind != "none":
        raise SystemExit(f"hole = {kind}: use circle, ellipse or none")
    opts = PlateOptions(width=float(m.get("width", 100.0)), height=float(m.get("height", 50.0)),
                        thickness=float(m.get("thickness", 1.0)), h=float(m["h"]), hole=hole,
                        ring_cells=int(m.get("ring_cells", 2)), n_r=_opt(m.get("n_r"), int),
                        grading=_opt(m.get("grading"), float))
    mesh = plate_mesh(opts)
    report = validate(mesh, cell_size=opts.h)
    summary = plate_summary(mesh)
    summary["validation"] = report.as_dict()
    return mesh, summary, report


def grain_ids(mesh, points_csv: Path) -> np.ndarray:
    """Nearest point-file grain id at each element centroid (PINNACLE's nearest interpolation)."""
    from scipy.spatial import cKDTree
    rows = list(csv.DictReader((ln for ln in points_csv.open() if not ln.startswith("#"))))
    xy = np.array([[float(r["x"]), float(r["y"])] for r in rows]); gid = np.array([int(r["grain_id"]) for r in rows])
    cent = mesh.points[mesh.cells[:, :4]].mean(axis=1)[:, :2]
    return gid[cKDTree(xy).query(cent)[1]]


def render(mesh, gid: np.ndarray | None, fig_path: Path, summary: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "mathtext.fontset": "dejavusans"})
    quads = mesh.points[mesh.cells[:, :4]][:, :, :2]
    hole = summary.get("hole")
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.4), dpi=220, gridspec_kw={"height_ratios": [1.0, 1.15]})
    for ax, zoom in zip(axes, (False, True)):
        if gid is not None:
            rng = np.random.default_rng(3)
            cols = plt.get_cmap("tab20")(rng.permutation(np.unique(gid).size % 20 + 20)[gid % 20] / 20.0) if False else plt.get_cmap("tab20")((gid * 7) % 20 / 20.0)
            pc = PolyCollection(quads, facecolors=cols, edgecolors="k", linewidths=0.25 if zoom else 0.12)
        else:
            pc = PolyCollection(quads, facecolors="#dfe6ee", edgecolors="k", linewidths=0.25 if zoom else 0.12)
        ax.add_collection(pc)
        ax.set_aspect("equal")
        if zoom and hole:
            cx, cy = hole["centre"]; kx, ky = hole["transition_cells"]; h = summary["h"]
            ax.set_xlim(cx - (kx + 1.5) * h, cx + (kx + 1.5) * h); ax.set_ylim(cy - (ky + 1.5) * h, cy + (ky + 1.5) * h)
            ax.set_title(f"O-grid ring: {hole['n_r']} layers, {hole['n_theta']} segments, first layer {hole['first_layer']:.2f} µm, hole arc {hole['arc_at_hole']:.2f} µm", fontsize=8.5)
        else:
            ax.set_xlim(0, summary["width"]); ax.set_ylim(0, summary["height"])
            ax.set_title(f"{summary['n_elements']} HEX8, {summary['n_nodes']} nodes, h = {summary['h']} µm, min scaled Jacobian {summary['min_scaled_jacobian']:.2f}", fontsize=8.5)
        ax.set_xlabel(r"$x$ (µm)"); ax.set_ylabel(r"$y$ (µm)")
        for s in ax.spines.values():
            s.set_visible(True)
    fig.tight_layout(); fig_path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(fig_path); plt.close(fig)


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); m = read_deck(deck)
    out = _resolve(str(m["out_dir"]), deck); out.mkdir(parents=True, exist_ok=True)
    mesh, summary, report = build(m, deck)
    from hexaform.io import exodus, vtk
    exodus.write(mesh, out / "mesh.e", title="plate with stress concentrator (hexaform plate)")
    vtk.write(mesh, out / "mesh.vtk")
    (out / "mesh.json").write_text(json.dumps(summary, indent=1, default=str))
    gid = grain_ids(mesh, _resolve(str(m["points"]), deck)) if m.get("points") else None
    if gid is not None:
        np.savetxt(out / "elem_grain.csv", np.column_stack([np.arange(len(gid)), gid]), fmt="%d", delimiter=",", header="element_id,grain_id", comments="")
    if m.get("figure"):
        render(mesh, gid, _resolve(str(m["figure"]), deck), summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "validation"}, default=str))
    print("validation passed:", report.passed, report.failures)
    return 0 if report.passed else 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
