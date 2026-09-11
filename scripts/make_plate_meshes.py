#!/usr/bin/env python3
"""Pseudo-3D perforated-plate meshes: 100 x 50 x 10 um, one element through z.

GEOMETRY. A rectangular plate with a central through-hole of radius 5 um, the
hole axis along z. One linear hex through the thickness, so the out-of-plane
kinematics are a constant through-thickness strain -- a generalized
plane-stress representation, with both z faces left traction-free. Mesh
density is refined only in x and y; the z direction is always one element, by
construction rather than by a parameter, so no density level can accidentally
add through-thickness DOF.

WHY A BLOCK-STRUCTURED O-GRID RATHER THAN A CARVED GRID. Deleting elements
whose centroid falls inside the hole leaves a staircase boundary, and the
whole point of the hole is that damage initiates at its edge: a staircase
puts spurious stress concentrations at every step and the initiation site
becomes a meshing artifact. The O-grid below is conforming -- every node on
the hole lies at exactly r = 5 -- and lets the radial layers be graded
independently of the background, which is what resolves l0 in the process
zone. The earlier spherical-void family taught this the hard way: refining
the background while the first element layer at the surface barely moved
converged the far field and left the damage band under-resolved.

BLOCK LAYOUT (all conforming, all structured, no distorted transition):

    +---------+---------------+---------+   y = 50
    | corner  |  top strip    | corner  |
    +---------+---------------+---------+   y = 40
    |         |               |         |
    |  left   |    O-grid     |  right  |
    |  strip  |  (hole here)  |  strip  |
    |         |               |         |
    +---------+---------------+---------+   y = 10
    | corner  | bottom strip  | corner  |
    +---------+---------------+---------+   y = 0
  x=0      x=35            x=65     x=100

The O-grid maps the circle to the square [35,65] x [10,40] by a transfinite
blend, so the four side strips and four corner blocks are plain rectangles.
Interfaces match because the square's side division count is shared with the
strips it touches.

Usage:
    make_plate_meshes.py --densities 1 2 3 4 5 [--out meshes/plate]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

WIDTH, HEIGHT, THICK = 100.0, 50.0, 10.0
HOLE_R = 5.0
CX, CY = WIDTH / 2.0, HEIGHT / 2.0
# Transition square around the hole; its half-extents set where the O-grid
# hands over to the rectangular blocks.
SQ_X0, SQ_X1 = 35.0, 65.0
SQ_Y0, SQ_Y1 = 10.0, 40.0


def _graded(n: int, ratio: float) -> np.ndarray:
    """``n+1`` points on [0,1] with spacing growing by ``ratio`` end to end."""
    if abs(ratio - 1.0) < 1e-12:
        return np.linspace(0.0, 1.0, n + 1)
    r = ratio ** (1.0 / max(1, n - 1))
    w = r ** np.arange(n)
    s = np.concatenate([[0.0], np.cumsum(w)])
    return s / s[-1]


def _square_perimeter(n_side: int) -> np.ndarray:
    """Points on the transition square, counter-clockwise from its +x mid-edge.

    Returns ``4*n_side`` points, matching the O-grid's angular resolution and
    the side strips' division counts, which is what makes the interfaces
    conform.
    """
    pts = []
    xs = np.linspace(SQ_X0, SQ_X1, n_side + 1)
    ys = np.linspace(SQ_Y0, SQ_Y1, n_side + 1)
    # right edge, bottom->top (excluding last), then top R->L, left T->B, bottom L->R
    pts += [(SQ_X1, y) for y in ys[:-1]]
    pts += [(x, SQ_Y1) for x in xs[::-1][:-1]]
    pts += [(SQ_X0, y) for y in ys[::-1][:-1]]
    pts += [(x, SQ_Y0) for x in xs[:-1]]
    return np.asarray(pts, dtype=float)


def _ogrid_nodes(n_side: int, n_r: int, grading: float) -> np.ndarray:
    """O-grid node grid of shape ``(n_r+1, 4*n_side, 2)``, hole ring first.

    Radial station 0 is the hole (exactly r = HOLE_R by construction, so the
    boundary is conforming rather than approximated); station ``n_r`` is the
    transition square.
    """
    sq = _square_perimeter(n_side)
    n_theta = sq.shape[0]
    # Angle of each square point about the hole center: the circle ring uses
    # the SAME angles, so radial lines are straight and the blend is smooth.
    ang = np.arctan2(sq[:, 1] - CY, sq[:, 0] - CX)
    circ = np.stack([CX + HOLE_R * np.cos(ang), CY + HOLE_R * np.sin(ang)], axis=1)
    t = _graded(n_r, grading)[:, None, None]
    return (1.0 - t) * circ[None, :, :] + t * sq[None, :, :]


def _rect_nodes(x0, x1, y0, y1, nx, ny):
    """Structured node grid on a rectangle, shape ``(nx+1, ny+1, 2)``."""
    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    return np.stack([gx, gy], axis=-1)


class _NodePool:
    """Deduplicating 2-D node pool, so blocks share interface nodes exactly."""

    def __init__(self, tol: float = 1.0e-9) -> None:
        self.pts: list[tuple[float, float]] = []
        self._index: dict[tuple[int, int], int] = {}
        self.tol = tol

    def add(self, x: float, y: float) -> int:
        key = (int(round(x / self.tol)), int(round(y / self.tol)))
        got = self._index.get(key)
        if got is not None:
            return got
        self._index[key] = len(self.pts)
        self.pts.append((float(x), float(y)))
        return self._index[key]


def build_plate(density: int, *, n_r: int, grading: float, hole: bool = True) -> dict:
    """Build one pseudo-3D plate mesh at the given in-plane density level.

    ``hole=False`` returns the same plate with no perforation: a plain
    structured grid. CP parameter calibration runs on that one, because a
    stress concentration would fold the geometry into the fitted hardening
    parameters -- the calibration needs the material's response, not the
    notch's.
    """
    n_side = 4 * density          # per side of the transition square
    nx_out = 3 * density          # x divisions in each side strip
    ny_out = max(1, density)      # y divisions in each corner/strip band

    pool = _NodePool()
    quads: list[tuple[int, int, int, int]] = []

    if not hole:
        # Plain plate: one structured block, same in-plane resolution as the
        # background of the perforated family so the two are comparable.
        nx = int(round(WIDTH / (HEIGHT / (2 * n_side))))
        ny = 2 * n_side
        g = _rect_nodes(0.0, WIDTH, 0.0, HEIGHT, nx, ny)
        ids = np.empty((nx + 1, ny + 1), dtype=int)
        for i in range(nx + 1):
            for j in range(ny + 1):
                ids[i, j] = pool.add(*g[i, j])
        for i in range(nx):
            for j in range(ny):
                quads.append((ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]))
        og = np.zeros((2, 1, 2))   # keeps the first-layer diagnostic well defined
    else:
        # --- O-grid around the hole ---
        og = _ogrid_nodes(n_side, n_r, grading)
        n_theta = og.shape[1]
        og_id = np.empty((n_r + 1, n_theta), dtype=int)
        for i in range(n_r + 1):
            for j in range(n_theta):
                og_id[i, j] = pool.add(*og[i, j])
        for i in range(n_r):
          for j in range(n_theta):
            jn = (j + 1) % n_theta
            # Radial-then-angular ordering, NOT angular-then-radial. The
            # square perimeter runs counter-clockwise, so theta increases with
            # j; taking the angular step first traverses +theta at the inner
            # radius, outward, then -theta, which encloses a NEGATIVE area
            # (theta_hat x r_hat = -z_hat). The rectangular blocks below wind
            # counter-clockwise, so mixing the two gives inverted hexes on
            # extrusion -- caught here by a triple-product sign check, not by
            # anything the solver would have reported.
            quads.append((og_id[i, j], og_id[i + 1, j], og_id[i + 1, jn], og_id[i, jn]))

    # --- eight rectangular blocks around the transition square ---
    blocks = [] if not hole else [
        (0.0, SQ_X0, SQ_Y0, SQ_Y1, nx_out, n_side),      # left strip
        (SQ_X1, WIDTH, SQ_Y0, SQ_Y1, nx_out, n_side),    # right strip
        (SQ_X0, SQ_X1, 0.0, SQ_Y0, n_side, ny_out),      # bottom strip
        (SQ_X0, SQ_X1, SQ_Y1, HEIGHT, n_side, ny_out),   # top strip
        (0.0, SQ_X0, 0.0, SQ_Y0, nx_out, ny_out),        # corners
        (SQ_X1, WIDTH, 0.0, SQ_Y0, nx_out, ny_out),
        (0.0, SQ_X0, SQ_Y1, HEIGHT, nx_out, ny_out),
        (SQ_X1, WIDTH, SQ_Y1, HEIGHT, nx_out, ny_out),
    ]
    for x0, x1, y0, y1, nx, ny in blocks:
        g = _rect_nodes(x0, x1, y0, y1, nx, ny)
        ids = np.empty((nx + 1, ny + 1), dtype=int)
        for i in range(nx + 1):
            for j in range(ny + 1):
                ids[i, j] = pool.add(*g[i, j])
        for i in range(nx):
            for j in range(ny):
                quads.append((ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]))

    pts2 = np.asarray(pool.pts, dtype=float)
    n2 = pts2.shape[0]

    # --- extrude one element through the thickness ---
    nodes = np.zeros((2 * n2, 3), dtype=float)
    nodes[:n2, :2] = pts2
    nodes[:n2, 2] = 0.0
    nodes[n2:, :2] = pts2
    nodes[n2:, 2] = THICK
    cells = np.asarray(
        [[a, b, c, d, a + n2, b + n2, c + n2, d + n2] for (a, b, c, d) in quads],
        dtype=int,
    )

    # --- node sets ---
    tol = 1.0e-7
    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    r = np.hypot(x - CX, y - CY)
    sets = {
        "xmin": np.flatnonzero(np.abs(x - 0.0) < tol),
        "xmax": np.flatnonzero(np.abs(x - WIDTH) < tol),
        "ymin": np.flatnonzero(np.abs(y - 0.0) < tol),
        "ymax": np.flatnonzero(np.abs(y - HEIGHT) < tol),
        "zmin": np.flatnonzero(np.abs(z - 0.0) < tol),
        "zmax": np.flatnonzero(np.abs(z - THICK) < tol),
        "hole": np.flatnonzero(np.abs(r - HOLE_R) < 1.0e-6),
        # The outer boundary excluding the hole: this is what a damage
        # Dirichlet condition should use if one is imposed at all. Pinning
        # damage on a set that included the hole would forbid initiation at
        # the exact place the hole exists to create it.
        # A SINGLE node, not a face. Plane stress leaves u_y free everywhere
        # (only u_x is prescribed on the loaded faces and u_z on one z face),
        # so a y-translation mode remains. Constraining a whole face in y
        # would suppress Poisson contraction there and quietly raise the
        # transverse stress -- the failure mode already seen on the cube
        # study -- so exactly one node is pinned. Rotation about z is already
        # removed by holding the whole xmin face at u_x = 0.
        "pin_y": np.flatnonzero(
            (np.abs(x) < tol) & (np.abs(y - HEIGHT / 2.0) < tol) & (np.abs(z) < tol)
        ),
        # Damage Dirichlet set: the LOADED ENDS ONLY. The crack initiates at
        # the hole and runs transversely to the free edges at y = 0 and
        # y = HEIGHT, so pinning d = 0 on the whole outer boundary would
        # forbid the crack from ever reaching the edge and would cap the
        # simulation short of final fracture. The grips are far from the
        # crack path and do not fracture in the experiment either.
        "grips": np.flatnonzero((np.abs(x) < tol) | (np.abs(x - WIDTH) < tol)),
        "outer": np.flatnonzero(
            (np.abs(x) < tol) | (np.abs(x - WIDTH) < tol)
            | (np.abs(y) < tol) | (np.abs(y - HEIGHT) < tol)
        ),
    }

    # --- quality: signed area of each quad, and the first radial layer ---
    q = pts2[np.asarray(quads, dtype=int)]
    # SIGNED shoelace area, not its absolute value: the sign is the check.
    # An inverted element has a perfectly respectable |area| and would pass a
    # magnitude test while giving the solver a negative Jacobian.
    qn = np.roll(q, -1, axis=1)
    area = 0.5 * np.sum(q[:, :, 0] * qn[:, :, 1] - qn[:, :, 0] * q[:, :, 1], axis=1)
    first_layer = float(np.min(np.linalg.norm(og[1] - og[0], axis=1)))
    edge = np.linalg.norm(q - np.roll(q, -1, axis=1), axis=2)

    return {
        "density": density,
        "n_elements": int(cells.shape[0]),
        "n_nodes": int(nodes.shape[0]),
        "n_dofs_mech": int(3 * nodes.shape[0]),
        "n_side": n_side,
        "n_r": n_r,
        "grading": grading,
        "first_layer_um": first_layer,
        "background_h_um": float(np.max(edge)),
        "min_signed_area": float(area.min()),
        "n_inverted": int(np.sum(area <= 0.0)),
        "max_aspect": float(np.max(edge.max(axis=1) / np.maximum(edge.min(axis=1), 1e-30))),
        "hole_nodes": int(sets["hole"].size),
        "nodes": nodes,
        "cells": cells,
        "sets": sets,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--densities", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--n-r", type=int, nargs="+", default=None,
                    help="Radial layers in the O-grid, one per density level.")
    ap.add_argument("--grading", type=float, default=3.0,
                    help="Radial spacing growth from the hole outward; >1 clusters "
                         "elements at the hole, where the damage band forms.")
    ap.add_argument("--out", default="meshes/plate")
    ap.add_argument("--no-hole", action="store_true",
                    help="Plain plate with no perforation, for CP calibration.")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    n_r_list = args.n_r or [4 * d for d in args.densities]
    if len(n_r_list) != len(args.densities):
        raise SystemExit("--n-r must give one value per density level")

    records = []
    for d, n_r in zip(args.densities, n_r_list):
        m = build_plate(d, n_r=n_r, grading=args.grading, hole=not args.no_hole)
        if m["n_inverted"]:
            raise SystemExit(
                f"density {d}: {m['n_inverted']} inverted elements "
                f"(min signed area {m['min_signed_area']:.4g}). Refusing to "
                "write a mesh the solver would silently accept."
            )
        tag = "plain" if args.no_hole else "hole"
        np.savez_compressed(
            out / f"plate_{tag}_d{d}.npz",
            nodes=m["nodes"], cells=m["cells"],
            **{f"set_{k}": v for k, v in m["sets"].items()},
        )
        rec = {k: v for k, v in m.items() if k not in ("nodes", "cells", "sets")}
        rec["file"] = f"plate_{tag}_d{d}.npz"
        rec["hole"] = not args.no_hole
        records.append(rec)
        print(f"d={d:2d}  elems {rec['n_elements']:6d}  nodes {rec['n_nodes']:6d}  "
              f"dofs {rec['n_dofs_mech']:7d}  first_layer {rec['first_layer_um']:.4f} um  "
              f"bg_h {rec['background_h_um']:.3f} um  aspect {rec['max_aspect']:.2f}  "
              f"hole_nodes {rec['hole_nodes']}")

    # The first radial layer must shrink with refinement, or the process zone
    # at the hole is not being refined even though the element count grows.
    fl = [r["first_layer_um"] for r in records]
    if not args.no_hole and not all(a > b for a, b in zip(fl, fl[1:])):
        print(f"WARNING: first layer not monotonically decreasing: {fl}")
    (out / ("plate_plain_family.json" if args.no_hole else "plate_family.json")).write_text(json.dumps(records, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
