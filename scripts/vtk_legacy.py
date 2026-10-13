#!/usr/bin/env python3
"""Legacy-ASCII VTK reader for PINNACLE field series (no plotting dependencies).

``read_vtk(path)`` returns nodes, cells, point and cell data; ``series_files(dir)``
lists the ``fields_step*.vtk`` frames in order. Used by the animation, the field
comparison and the surrogate data collection.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np


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

