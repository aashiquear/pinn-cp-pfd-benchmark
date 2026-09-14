#!/usr/bin/env python3
"""Voronoi against relaxed morphology: the maps, and the numbers behind them.

Six panels. Top row is the grain map before and after relaxation plus the
junction-angle distributions; bottom row is the grain-count and free-energy
trajectories and the metric table. The distributions carry the claim -- the
maps only illustrate it -- because "looks natural" is settled by the spread
of triple-junction angles about 120 degrees, not by appearance.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grain_morphology import summarize, triple_junction_angles  # noqa: E402
from map_grains import label_from_eta, voronoi_map  # noqa: E402

STUDY = Path(__file__).resolve().parent.parent
PLATE = (100.0, 50.0)


def load_case(n: int):
    """Grain maps before and after, on a common grid."""
    z = np.load(STUDY / f"results/graingrowth/n{n}/final_field.npz")
    field = np.asarray(z["field"], dtype=float)
    hist = json.loads((STUDY / f"results/graingrowth/n{n}/growth_history.json").read_text())
    div_from_nodes = None
    # recover the grid from the deck so the reshape is not guessed
    deck = (STUDY / f"results/graingrowth/n{n}/deck.i").read_text()
    for line in deck.splitlines():
        if "divisions" in line and "=" in line:
            a, b = line.split("=", 1)[1].strip().strip("'").split()
            div_from_nodes = (int(a), int(b))
            break
    if div_from_nodes is None:
        raise ValueError("divisions not found in deck")
    gid_new, n_new = label_from_eta(field, div_from_nodes)
    seeds = np.asarray(np.load(STUDY / f"meshes/plate/poly_d3_n{n}.npz")["seeds"], float)
    gid_old = voronoi_map(seeds, div_from_nodes)
    hx = PLATE[0] / div_from_nodes[0]
    hy = PLATE[1] / div_from_nodes[1]
    return dict(n=n, div=div_from_nodes, hx=hx, hy=hy, gid_old=gid_old,
                gid_new=gid_new, n_new=n_new, hist=hist, field=field)
