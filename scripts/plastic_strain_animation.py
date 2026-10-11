#!/usr/bin/env python3
"""Equivalent plastic strain on the deformed plate, animated over the loading.

Reads the ``fields_step{NNNN}.vtk`` series PINNACLE writes when the deck sets
``snapshots_every`` (decks/cp_study/fields_h*.i), one series per mesh, and
draws every element as a quad in the deformed configuration (nodal
displacement from the file, magnified by ``scale``), colored by the
element-averaged equivalent plastic strain. One panel per mesh, one frame per
snapshot, a common color scale, the nominal strain of the frame in the title.

Also writes the still frames the report uses (``--stills``: nominal strains
in percent; the nearest recorded snapshot is taken).

Driven by a deck (PINNACLE ``.i`` syntax):

  [Animation]
    out_gif    = ../../figures/report/plastic_strain_mesh.gif
    out_stills = ../../figures/report/plastic_strain_mesh_{strain}.png
    stills     = '2 8 15 20'        # nominal strain (%) of each still
    strain_max = 0.20               # nominal strain at load factor 1
    vmax       = 0.6                # top of the color scale (equivalent plastic strain)
    scale      = 1.0                # displacement magnification
    fps        = 8
  []
  [Series]
    [./h20]  dir = ../../results/cp_study/fields/h20/solver   label = h20 []
    ...
  []

Usage: plastic_strain_animation.py decks/report/plastic_strain_animation.i
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

STUDY = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- deck
def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("animation", "series"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _list(v) -> list[float]:
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(x) for x in str(v).split()]


# --------------------------------------------------------------------------- vtk
def read_vtk(path: Path) -> dict:
    """Legacy ASCII unstructured grid as PINNACLE writes it.

    Returns nodes (n, 3), cells (m, k), point vectors and cell scalars by name.
    Only the records the exporter emits are parsed (POINTS, CELLS, CELL_TYPES,
    POINT_DATA/CELL_DATA with SCALARS, VECTORS, TENSORS).
    """
    tok = path.read_text().split("\n")
    i, n = 0, len(tok)
    out = {"point": {}, "cell": {}}
    section = None
    n_points = n_cells = 0
    while i < n:
        parts = tok[i].split()
        if not parts:
            i += 1; continue
        key = parts[0]
        if key == "POINTS":
            n_points = int(parts[1]); vals = []
            i += 1
            while len(vals) < 3 * n_points:
                vals.extend(float(v) for v in tok[i].split()); i += 1
            out["nodes"] = np.asarray(vals).reshape(n_points, 3); continue
        if key == "CELLS":
            n_cells = int(parts[1]); rows = []
            i += 1
            for _ in range(n_cells):
                r = [int(v) for v in tok[i].split()]; rows.append(r[1:1 + r[0]]); i += 1
            out["cells"] = np.asarray(rows); continue
        if key == "CELL_TYPES":
            i += 1 + n_cells; continue
        if key == "POINT_DATA":
            section = "point"; i += 1; continue
        if key == "CELL_DATA":
            section = "cell"; i += 1; continue
        if key in ("SCALARS", "VECTORS", "TENSORS"):
            name = parts[1]
            count = n_points if section == "point" else n_cells
            width = {"SCALARS": 1, "VECTORS": 3, "TENSORS": 9}[key]
            i += 1
            if key == "SCALARS" and tok[i].startswith("LOOKUP_TABLE"):
                i += 1
            vals = []
            while len(vals) < width * count:
                vals.extend(float(v) for v in tok[i].split()); i += 1
            arr = np.asarray(vals).reshape(count, width) if width > 1 else np.asarray(vals)
            out[section][name] = arr; continue
        i += 1
    return out


def series_files(d: Path) -> list[tuple[int, Path]]:
    files = []
    for p in sorted(d.glob("fields_step*.vtk")):
        m = re.search(r"fields_step(\d+)\.vtk", p.name)
        if m:
            files.append((int(m.group(1)), p))
    if not files:
        raise SystemExit(f"{d}: no fields_step*.vtk (did the deck set snapshots_every?)")
    return files


def pvd_times(d: Path) -> dict[int, float]:
    """step index (file order) -> time from fields.pvd, if present."""
    p = d / "fields.pvd"
    if not p.exists():
        return {}
    out = {}
    for m in re.finditer(r'timestep="([^"]+)"[^>]*file="fields_step(\d+)\.vtk"', p.read_text()):
        out[int(m.group(2))] = float(m.group(1))
    return out


# --------------------------------------------------------------------------- draw
def top_quads(nodes: np.ndarray, cells: np.ndarray, disp: np.ndarray, scale: float) -> np.ndarray:
    """In-plane quads of the HEX8 plate: the four nodes with the largest z, in cycle order."""
    xy = nodes[:, :2] + scale * disp[:, :2]
    quads = []
    for c in cells:
        z = nodes[c, 2]
        top = c[np.argsort(z)[-4:]] if len(c) == 8 else c
        pts = xy[top]
        ang = np.arctan2(pts[:, 1] - pts[:, 1].mean(), pts[:, 0] - pts[:, 0].mean())
        quads.append(pts[np.argsort(ang)])
    return np.asarray(quads)


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    an = cfg["animation"]; series = cfg["series"]
    strain_max = float(an.get("strain_max", 0.20))
    vmax = float(an.get("vmax", 0.6)); scale = float(an.get("scale", 1.0))
    stills = _list(an.get("stills", "2 8 15 20"))

    # Load every series: frames keyed by step number.
    loaded = {}
    for name, s in series.items():
        d = _resolve(str(s["dir"]), deck)
        files = series_files(d)
        frames = []
        for k, p in files:
            v = read_vtk(p)
            frames.append({"index": k, "nodes": v["nodes"], "cells": v["cells"],
                           "disp": v["point"]["displacement"],
                           "peeq": v["cell"]["equivalent_plastic_strain"]})
        loaded[name] = {"label": str(s.get("label", name)), "frames": frames}
    n_frames = min(len(v["frames"]) for v in loaded.values())
    # Nominal strain of frame j: the load factor is written in the history, but
    # the series is uniform in load factor, so eps_j = strain_max * (j+1)/n_frames
    # when the run reached the end (the deck's n_steps / snapshots_every frames).
    strains = [strain_max * (j + 1) / n_frames for j in range(n_frames)]

    names = list(loaded)
    ncol = 2; nrow = int(np.ceil(len(names) / ncol))
    cmap = plt.get_cmap("viridis")

    def draw(j: int, fig_path: Path, dpi: int):
        fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.9 * nrow), dpi=dpi)
        for ax, name in zip(np.ravel(axes), names):
            fr = loaded[name]["frames"][j]
            quads = top_quads(fr["nodes"], fr["cells"], fr["disp"], scale)
            pc = PolyCollection(quads, array=fr["peeq"], cmap=cmap, edgecolors="none")
            pc.set_clim(0.0, vmax)
            ax.add_collection(pc)
            ax.set_xlim(-2, 102 + 100 * strain_max * scale + 2); ax.set_ylim(-2, 52)
            ax.set_aspect("equal"); ax.set_xticks([0, 50, 100]); ax.set_yticks([0, 25, 50])
            ax.set_title(loaded[name]["label"], fontsize=9, loc="left")
        for ax in np.ravel(axes)[len(names):]:
            ax.axis("off")
        fig.suptitle(r"Equivalent plastic strain, $\varepsilon$ = %.1f%%" % (100 * strains[j]), fontsize=10)
        cb = fig.colorbar(pc, ax=np.ravel(axes).tolist(), shrink=0.8, pad=0.02)
        cb.set_label(r"$\bar{\varepsilon}^p$")
        fig.savefig(fig_path); plt.close(fig)

    out_gif = _resolve(str(an["out_gif"]), deck); out_gif.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_gif.parent / "_frames"; tmp.mkdir(exist_ok=True)
    frame_paths = []
    for j in range(n_frames):
        p = tmp / f"frame_{j:03d}.png"; draw(j, p, 90); frame_paths.append(p)
    try:
        from PIL import Image
        imgs = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in frame_paths]
        imgs[0].save(out_gif, save_all=True, append_images=imgs[1:],
                     duration=int(1000 / float(an.get("fps", 8))), loop=0)
        print("wrote", out_gif, f"({n_frames} frames)")
    except ImportError:
        raise SystemExit("Pillow is needed to assemble the GIF (pip install pillow)")
    for p in frame_paths:
        p.unlink()
    tmp.rmdir()
    if "out_stills" in an:
        for eps in stills:
            j = int(np.argmin([abs(100 * s - eps) for s in strains]))
            p = _resolve(str(an["out_stills"]).replace("{strain}", f"{eps:g}"), deck)
            draw(j, p, 200); print("wrote", p, f"(frame {j}, {100 * strains[j]:.2f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
