#!/usr/bin/env python3
"""Map a relaxed grain structure onto the mechanics mesh, and compare morphology.

The growth solve runs on its own fine grid because the diffuse boundary needs
h <= w_GB/3 while the mechanics mesh is set by the mechanics. Transfer is by
nearest-cell lookup at each mechanics element centroid, which is exact for a
labelled (piecewise-constant) grain field -- no interpolation is meaningful
across a grain boundary.

Orientations follow the grain, not the cell: each relaxed grain inherits the
orientation of whichever original Voronoi grain contributed most of its area,
so the relaxed polycrystal is the SAME polycrystal evolved and the IPF
coverage measured for the Voronoi set still describes it. Where relaxation has
removed a grain its orientation simply disappears; that is reported rather
than back-filled.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grain_morphology import summarize  # noqa: E402

STUDY = Path(__file__).resolve().parent.parent
PLATE = (100.0, 50.0)


def label_from_eta(eta: np.ndarray, div: tuple[int, int], threshold: float = 0.5):
    """Grain ids by connected-component labelling within each colour.

    Colouring puts several NON-ADJACENT grains on one order parameter, so an
    argmax over order parameters returns the colour, not the grain. Grains are
    the connected components of each colour's own indicator.
    """
    nx, ny = div[0] + 1, div[1] + 1
    if eta.ndim != 2 or eta.shape[0] != nx * ny:
        raise ValueError(f"eta shape {eta.shape} does not match grid {(nx, ny)}")
    col = np.argmax(eta, axis=1).reshape(nx, ny, order="C")
    val = np.max(eta, axis=1).reshape(nx, ny, order="C")
    gid = np.full((nx, ny), -1, dtype=int)
    nxt = 0
    for k in range(eta.shape[1]):
        m = (col == k) & (val >= threshold)
        if not m.any():
            continue
        lab, n = ndimage.label(m)
        for i in range(1, n + 1):
            sel = lab == i
            # drop single-cell specks: they are labelling noise on the diffuse
            # boundary, not grains
            if sel.sum() < 4:
                continue
            gid[sel] = nxt
            nxt += 1
    # fill any unassigned cell (inside a diffuse boundary) from its nearest
    # labelled neighbour, so the map is a partition
    if (gid < 0).any():
        known = np.argwhere(gid >= 0)
        if len(known):
            tree = cKDTree(known)
            miss = np.argwhere(gid < 0)
            _, j = tree.query(miss)
            gid[miss[:, 0], miss[:, 1]] = gid[known[j, 0], known[j, 1]]
    return gid, int(gid.max()) + 1


def voronoi_map(seeds: np.ndarray, div: tuple[int, int]):
    """The unrelaxed structure on the same grid, for a like-for-like comparison."""
    nx, ny = div[0] + 1, div[1] + 1
    xs = np.linspace(0.0, PLATE[0], nx)
    ys = np.linspace(0.0, PLATE[1], ny)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return cKDTree(seeds).query(np.column_stack([X.ravel(), Y.ravel()]))[1].reshape(nx, ny)


def inherit_orientations(gid_new: np.ndarray, gid_old: np.ndarray,
                         euler_old: np.ndarray):
    """Each relaxed grain takes the orientation of its majority ancestor."""
    out, prov = {}, []
    for g in range(int(gid_new.max()) + 1):
        m = gid_new == g
        if not m.any():
            continue
        anc, cnt = np.unique(gid_old[m], return_counts=True)
        win = int(anc[np.argmax(cnt)])
        frac = float(cnt.max() / cnt.sum())
        out[g] = euler_old[win]
        prov.append({"new_gid": g, "ancestor": win, "majority_fraction": frac,
                     "n_ancestors": int(len(anc))})
    return out, prov


def to_mechanics(gid: np.ndarray, div: tuple[int, int], nodes: np.ndarray,
                 cells: np.ndarray):
    """Grain id per mechanics ELEMENT, by nearest grid cell to its centroid."""
    cen = nodes[cells].mean(axis=1)
    nx, ny = div[0] + 1, div[1] + 1
    ix = np.clip(np.round(cen[:, 0] / PLATE[0] * (nx - 1)).astype(int), 0, nx - 1)
    iy = np.clip(np.round(cen[:, 1] / PLATE[1] * (ny - 1)).astype(int), 0, ny - 1)
    return gid[ix, iy]
