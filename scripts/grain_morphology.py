#!/usr/bin/env python3
"""Measure whether a relaxed grain structure is actually natural.

"Looks natural" has testable content, and three measures carry it. All are
computed on the labelled grain map, before and after relaxation, so the claim
rests on numbers rather than on the picture.

1. TRIPLE-JUNCTION ANGLES. In an isotropic polycrystal at local equilibrium
   the three boundaries meeting at a junction make 120 degrees (Young's law
   with equal boundary energies). A Voronoi tessellation has no such
   constraint -- its junctions are set by seed positions -- so the spread
   about 120 degrees is the sharpest single discriminator.

2. NEIGHBOUR COUNT. A space-filling 2-D network averages exactly 6 sides per
   grain (Euler's relation), both before and after, so the MEAN is not
   informative -- the DISTRIBUTION is: growth removes many-sided and few-sided
   grains preferentially, narrowing it.

3. BOUNDARY CURVATURE / PERIMETER. Curvature-driven growth reduces total
   boundary length monotonically and replaces straight Voronoi edges with
   curved ones. Circularity 4*pi*A/P^2 rises toward the equiaxed limit.

A caveat carried into the output: at these grain counts the sample is small
(8-32 grains), so the distributions are indicative, not statistically
converged. Junction angles are the exception -- each grain contributes several
junctions, so there are enough of them to compare meaningfully.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import ndimage

STUDY = Path(__file__).resolve().parent.parent



def perimeter_crofton(m: np.ndarray, hx: float, hy: float) -> float:
    """Boundary length of a binary region by the 4-direction Crofton estimate.

    ``P = (pi/2) * mean_theta (N_theta * s_theta)``, with ``N`` the number of
    boundary intersections along a family of parallel test lines and ``s`` that
    family's perpendicular spacing. Counting face crossings alone (the
    city-block perimeter) overestimates by about 24% -- it measures the
    staircase, not the boundary -- which drove the circularity of the hexagonal
    control to 0.5885 against its analytic 0.9069.

    Validated against shapes with known perimeter: disc +0.21%, square -4.52%,
    diamond -4.52%, and the regular-hexagon circularity +1.79%. The residual
    few percent on straight edges is the known bias of a 4-direction estimate.
    It assumes near-isotropic cells; at hx/hy = 2 the error grows to +13%, so
    the growth grids are built isotropic (h equal on both axes to 3 figures).
    """
    b = m.astype(np.int8)
    d = float(np.hypot(hx, hy))
    s_d = hx * hy / d
    n_x = int(np.abs(np.diff(b, axis=0)).sum())
    n_y = int(np.abs(np.diff(b, axis=1)).sum())
    n_d1 = int(np.abs(b[1:, 1:] - b[:-1, :-1]).sum())
    n_d2 = int(np.abs(b[1:, :-1] - b[:-1, 1:]).sum())
    return float((np.pi / 2.0) * (n_x * hy + n_y * hx + (n_d1 + n_d2) * s_d) / 4.0)


def grain_props(gid: np.ndarray, h: tuple[float, float],
                interior_only: bool = True) -> dict:
    """Per-grain area, perimeter and circularity from a labelled map.

    Restricted to interior grains by default. A grain cut by the domain edge
    has the cut portion of its boundary uncounted, so its perimeter is
    truncated while its area is not -- which drove circularity to 1.0585 for
    the 8-grain structure, above the geometric maximum of 1, and is a
    measurement artifact rather than a shape.
    """
    ids = np.unique(gid[gid >= 0])
    if interior_only:
        edge = set()
        for band in (gid[0, :], gid[-1, :], gid[:, 0], gid[:, -1]):
            edge.update(int(v) for v in np.unique(band) if v >= 0)
        keep = [i for i in ids if int(i) not in edge]
        ids = np.array(keep, dtype=int) if len(keep) >= 1 else np.array([], dtype=int)
    cell = h[0] * h[1]
    out = []
    for i in ids:
        m = gid == i
        area = float(m.sum()) * cell
        per = perimeter_crofton(m, h[0], h[1])
        out.append(dict(gid=int(i), area=area, perimeter=per,
                        circularity=(4 * np.pi * area / per**2) if per > 0 else np.nan,
                        equiv_diam=2.0 * np.sqrt(area / np.pi)))
    return {"grains": out, "n": len(out)}


def neighbour_counts(gid: np.ndarray, interior_only: bool = True) -> dict:
    """How many distinct grains each grain touches.

    Grains touching the domain edge are excluded by default: they are cut off
    by the boundary and have artificially few neighbours, which pulls the mean
    below the value Euler's relation requires (6 for a space-filling 2-D
    network) and would show up as a spurious difference between structures
    with different numbers of edge grains.
    """
    ids = np.unique(gid[gid >= 0])
    edge = set()
    if interior_only:
        for band in (gid[0, :], gid[-1, :], gid[:, 0], gid[:, -1]):
            edge.update(int(v) for v in np.unique(band) if v >= 0)
    nb = {int(i): set() for i in ids}
    for a, b in ((gid[1:, :], gid[:-1, :]), (gid[:, 1:], gid[:, :-1])):
        d = a != b
        for u, v in zip(a[d].ravel(), b[d].ravel()):
            if u >= 0 and v >= 0:
                nb[int(u)].add(int(v)); nb[int(v)].add(int(u))
    keep = [int(i) for i in ids if int(i) not in edge]
    # No silent fallback. When every grain touches the domain edge there is no
    # interior population and the neighbour statistic is NOT measurable: an
    # earlier version quietly fell back to all grains, which reported a mean of
    # 3.50 against the 6 Euler's relation requires and looked like a result.
    if len(keep) < 2:
        return {"counts": [], "n_interior": 0, "mean": float("nan"),
                "std": float("nan"), "min": -1, "max": -1,
                "note": "no interior grains: every grain touches the domain edge"}
    c = np.array([len(nb[i]) for i in keep], dtype=float)
    return {"counts": c.tolist(), "n_interior": len(keep), "mean": float(c.mean()),
            "std": float(c.std()), "min": int(c.min()), "max": int(c.max())}


def triple_junction_angles(gid: np.ndarray, hx: float, hy: float,
                           radius: float = 8.0) -> dict:
    """Angles between the three BOUNDARIES meeting at a triple junction.

    The angle must be measured between boundary TANGENTS, not between
    directions to grain centroids: a centroid direction is set by the whole
    grain's shape, and on a regular hexagonal tiling -- where the answer is
    exactly 120 degrees -- that estimator returns 16.6 degrees of RMS scatter,
    which is larger than the effect being looked for. Here each of the three
    boundaries is located directly, as the midpoints of adjacent pixel pairs
    carrying that grain pair, and its direction is the mean offset of those
    midpoints from the junction.

    Only the SPREAD about 120 degrees is informative. The mean is 120 by
    construction, since the three angles at a junction sum to 360.
    """
    nx, ny = gid.shape
    # candidate junctions: 2x2 windows with exactly three distinct grains
    cand = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            u = set(int(v) for v in gid[i:i+2, j:j+2].ravel() if v >= 0)
            if len(u) == 3:
                cand.append((i + 0.5, j + 0.5, tuple(sorted(u))))
    # cluster candidates of the same grain triple that are close together
    used, junc = [], []
    for cx, cy, key in cand:
        if any(k == key and (cx - x) ** 2 + (cy - y) ** 2 < radius ** 2
               for x, y, k in used):
            continue
        used.append((cx, cy, key)); junc.append((cx, cy, key))

    # precompute boundary midpoints per grain pair
    angs = []
    xy = []
    for cx, cy, key in junc:
        i0, i1 = max(0, int(cx - radius)), min(nx, int(cx + radius) + 2)
        j0, j1 = max(0, int(cy - radius)), min(ny, int(cy + radius) + 2)
        dirs = {}
        for pair in ((key[0], key[1]), (key[1], key[2]), (key[0], key[2])):
            pts = []
            for ii in range(i0, i1 - 1):
                for jj in range(j0, j1):
                    a_, b_ = int(gid[ii, jj]), int(gid[ii + 1, jj])
                    if {a_, b_} == set(pair):
                        pts.append((ii + 0.5, jj))
            for ii in range(i0, i1):
                for jj in range(j0, j1 - 1):
                    a_, b_ = int(gid[ii, jj]), int(gid[ii, jj + 1])
                    if {a_, b_} == set(pair):
                        pts.append((ii, jj + 0.5))
            if len(pts) < 3:
                dirs = {}; break
            p = np.asarray(pts, dtype=float)
            r2 = (p[:, 0] - cx) ** 2 + (p[:, 1] - cy) ** 2
            p = p[r2 <= radius ** 2]
            if len(p) < 3:
                dirs = {}; break
            dx = (p[:, 0].mean() - cx) * hx
            dy = (p[:, 1].mean() - cy) * hy
            n = np.hypot(dx, dy)
            if n < 1e-12:
                dirs = {}; break
            dirs[pair] = (dx / n, dy / n)
        if len(dirs) != 3:
            continue
        th = sorted(np.degrees(np.arctan2(d[1], d[0])) % 360.0 for d in dirs.values())
        angs.extend([th[1] - th[0], th[2] - th[1], 360.0 - th[2] + th[0]])
        # Position in MICROMETRES, one entry per angle (three per junction), so
        # a caller can select a spatial subset -- e.g. junctions in a band along
        # the domain edge, where an edge condition on the diffusion distorts
        # them and a domain-wide average would hide it.
        xy.extend([(cx * hx, cy * hy)] * 3)
    if not angs:
        return {"n_junctions": len(junc), "angles": [], "rms_dev_from_120": np.nan,
                "frac_within_15deg": np.nan, "std": np.nan, "angle_xy": []}
    a = np.asarray(angs, dtype=float)
    return {"n_junctions": len(junc), "angles": a.tolist(),
            "mean": float(a.mean()), "std": float(a.std()),
            "rms_dev_from_120": float(np.sqrt(np.mean((a - 120.0) ** 2))),
            "frac_within_15deg": float(np.mean(np.abs(a - 120.0) <= 15.0)),
            "angle_xy": [[float(u), float(v)] for u, v in xy]}


def summarize(gid: np.ndarray, hx: float, hy: float, label: str) -> dict:
    gp = grain_props(gid, (hx, hy))
    nb = neighbour_counts(gid)
    tj = triple_junction_angles(gid, hx, hy)
    n_all = int(len(np.unique(gid[gid >= 0])))
    if gp["n"] == 0:
        return {"label": label, "n_grains": n_all, "n_interior": 0,
                "shape_metrics": "not measurable: no interior grains",
                "total_boundary_length_um": float(
                    sum(g["perimeter"] for g in grain_props(
                        gid, (hx, hy), interior_only=False)["grains"]) / 2.0),
                "n_junctions": tj["n_junctions"],
                "junction_rms_dev_from_120": tj["rms_dev_from_120"],
                "junction_frac_within_15deg": tj.get("frac_within_15deg")}
    area = np.array([g["area"] for g in gp["grains"]])
    circ = np.array([g["circularity"] for g in gp["grains"]])
    per = np.array([g["perimeter"] for g in gp["grains"]])
    return {
        "label": label, "n_grains": n_all, "n_shape_grains": gp["n"],
        "total_boundary_length_um": float(
            sum(g["perimeter"] for g in grain_props(
                gid, (hx, hy), interior_only=False)["grains"]) / 2.0),
        "area_mean_um2": float(area.mean()), "area_cv": float(area.std() / area.mean()),
        "circularity_mean": float(np.nanmean(circ)),
        "circularity_std": float(np.nanstd(circ)),
        "neighbours_mean": nb["mean"], "neighbours_std": nb["std"],
        "neighbours_min": nb["min"], "neighbours_max": nb["max"], "n_interior": nb["n_interior"],
        "n_junctions": tj["n_junctions"],
        "junction_angle_mean": tj["mean"],
        "junction_rms_dev_from_120": tj["rms_dev_from_120"],
        "junction_frac_within_15deg": tj.get("frac_within_15deg"),
    }
