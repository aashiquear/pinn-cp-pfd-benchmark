#!/usr/bin/env python3
"""Polycrystal grain structures with IPF-covering random orientations.

Three sets are generated: 8, 16 and 32 grains on the cuboid plain plate.

TWO DELIBERATE CHOICES, both recorded because they affect interpretation.

1. Orientations are drawn HAAR-UNIFORM (uniform on SO(3) via quaternion
   sampling), not uniform in Euler angles -- a naive uniform draw in
   (phi1, Phi, phi2) over-samples the poles and is not a random texture.

2. From a large Haar-uniform pool a MAXIMIN-SPREAD subset is selected, so the
   standard triangle is covered rather than clumped. At 8 grains a pure random
   draw routinely leaves half the triangle empty, and the purpose of these runs
   is to see how orientation changes slip, so coverage is the requirement. The
   cost is a mild anti-clustering bias: these sets are NOT a statistically
   representative texture sample and must not be used to estimate texture
   statistics. The pure-random coverage is measured alongside for comparison.

GRAIN SHAPE. The plate is one element thick and grains are Voronoi cells of the
in-plane element centroids, so every grain spans the full thickness as a prism.
With cuboid (near-cubic) elements the in-plane grain size is a multiple of the
element size while the through-thickness size is exactly one element, so the
grains are NOT equiaxed -- they are flattened in z by roughly the in-plane
grain width in elements. The measured aspect ratio is reported per set.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

STUDY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(STUDY / "scripts"))


def haar_uniform_euler(n: int, seed: int) -> np.ndarray:
    """Haar-uniform orientations as Bunge ZXZ Euler angles (degrees)."""
    from pinnacle.materials.orientations import sample_orientations
    return sample_orientations({"type": "mackenzie"}, n, seed=seed)


def euler_to_matrix(ang: np.ndarray) -> np.ndarray:
    phi1, Phi, phi2 = np.radians(ang)
    c1, s1 = np.cos(phi1), np.sin(phi1)
    cP, sP = np.cos(Phi), np.sin(Phi)
    c2, s2 = np.cos(phi2), np.sin(phi2)
    Rz1 = np.array([[c1, -s1, 0.0], [s1, c1, 0.0], [0.0, 0.0, 1.0]])
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, cP, -sP], [0.0, sP, cP]])
    Rz2 = np.array([[c2, -s2, 0.0], [s2, c2, 0.0], [0.0, 0.0, 1.0]])
    return Rz2 @ Rx @ Rz1


def loading_direction(ang: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Crystal direction lying along the sample loading axis."""
    return euler_to_matrix(ang).T @ axis


def triangle_xy(d: np.ndarray) -> np.ndarray:
    """Standard-triangle stereographic coordinates of a direction."""
    d = np.asarray(d, float)
    d = d / np.linalg.norm(d)
    z, x, y = np.sort(np.abs(d))[::-1]
    return np.array([x / (1.0 + z), y / (1.0 + z)])


def triangle_mask(nx: int = 70):
    """Boolean grid over the bounding box, True inside the standard triangle."""
    xs = np.linspace(0.0, 0.4143, nx)
    ys = np.linspace(0.0, 0.3661, nx)
    GX, GY = np.meshgrid(xs, ys, indexing="ij")
    inside = np.zeros_like(GX, dtype=bool)
    for i in range(nx):
        for j in range(nx):
            X, Y = GX[i, j], GY[i, j]
            s = X * X + Y * Y
            v = np.array([2 * X / (1 + s), 2 * Y / (1 + s), (1 - s) / (1 + s)])
            v = v / np.linalg.norm(v)
            # inside the standard triangle iff z >= x >= y >= 0 for the
            # representative built from the sorted components
            zz, xx, yy = np.sort(np.abs(v))[::-1]
            inside[i, j] = np.allclose([xx / (1 + zz), yy / (1 + zz)], [X, Y], atol=2e-3)
    return GX, GY, inside


def coverage(pts_xy: np.ndarray, GX, GY, inside, tol: float = 0.045) -> float:
    """Fraction of triangle area within ``tol`` (projection units) of a point."""
    grid = np.column_stack([GX[inside], GY[inside]])
    if len(pts_xy) == 0:
        return 0.0
    d = cKDTree(pts_xy).query(grid)[0]
    return float((d <= tol).mean())


def maximin_subset(pool_xy: np.ndarray, n: int, seed: int) -> np.ndarray:
    """Greedy maximin selection: each pick is farthest from those already in."""
    rng = np.random.default_rng(seed)
    idx = [int(rng.integers(len(pool_xy)))]
    while len(idx) < n:
        d = cKDTree(pool_xy[idx]).query(pool_xy)[0]
        idx.append(int(np.argmax(d)))
    return np.array(idx)


def voronoi_grains(cent2: np.ndarray, n: int, seed: int):
    rng = np.random.default_rng(seed)
    lo, hi = cent2.min(axis=0), cent2.max(axis=0)
    seeds = rng.uniform(lo, hi, size=(n, 2))
    gid = cKDTree(seeds).query(cent2)[1]
    return gid.astype(int), seeds


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--density", type=int, default=3)
    ap.add_argument("--grain-sets", type=int, nargs="+", default=[8, 16, 32])
    ap.add_argument("--grain-seed", type=int, default=20260914)
    ap.add_argument("--orientation-seed", type=int, default=20260914)
    ap.add_argument("--pool", type=int, default=4000,
                    help="Haar-uniform pool size the maximin subset is drawn from.")
    ap.add_argument("--out", default="meshes/plate")
    args = ap.parse_args()

    mesh = np.load(STUDY / f"meshes/plate/plate_plain_d{args.density}.npz")
    nodes, cells = mesh["nodes"], mesh["cells"]
    cent = nodes[cells].mean(axis=1)
    thick = float(nodes[:, 2].max() - nodes[:, 2].min())
    # representative in-plane element edge
    e0 = np.linalg.norm(nodes[cells[:, 1]] - nodes[cells[:, 0]], axis=1)
    h_el = float(np.median(e0))

    GX, GY, inside = triangle_mask()
    # one Haar-uniform pool, reused for every set so the sets are nested in
    # construction and differences between them are grain count alone
    pool_eul = haar_uniform_euler(args.pool, args.orientation_seed)
    X_AXIS = np.array([1.0, 0.0, 0.0])   # plate decks pull along x
    pool_xy = np.array([triangle_xy(loading_direction(a, X_AXIS)) for a in pool_eul])

    out = Path(STUDY / args.out)
    records = []
    for n in args.grain_sets:
        gid, seeds = voronoi_grains(cent[:, :2], n, args.grain_seed)
        counts = np.bincount(gid, minlength=n)
        sel = maximin_subset(pool_xy, n, args.orientation_seed)
        eul = pool_eul[sel]
        xy = pool_xy[sel]
        # pure-random comparison at the same n
        rnd = haar_uniform_euler(n, args.orientation_seed + 1)
        rnd_xy = np.array([triangle_xy(loading_direction(a, X_AXIS)) for a in rnd])

        # grain geometry: in-plane equivalent width from element count
        width = np.sqrt(counts * h_el * h_el)
        aspect = width / thick
        rec = {
            "n_grains": int(n),
            "density": int(args.density),
            "n_elements": int(len(cells)),
            "element_edge_um": h_el,
            "thickness_um": thick,
            "elements_per_grain_min": int(counts.min()),
            "elements_per_grain_median": float(np.median(counts)),
            "elements_per_grain_max": int(counts.max()),
            "empty_grains": int((counts == 0).sum()),
            "grain_width_um_median": float(np.median(width)),
            "aspect_inplane_over_thickness_median": float(np.median(aspect)),
            "ipf_coverage_maximin": coverage(xy, GX, GY, inside),
            "ipf_coverage_pure_random": coverage(rnd_xy, GX, GY, inside),
            "grain_seed": args.grain_seed,
            "orientation_seed": args.orientation_seed,
            "euler_angles": eul.tolist(),
            "loading_directions": [loading_direction(a, X_AXIS).tolist() for a in eul],
            "triangle_xy": xy.tolist(),
        }
        records.append(rec)
        np.savez_compressed(out / f"poly_d{args.density}_n{n}.npz",
                            grain_id=gid, euler_angles=eul, seeds=seeds,
                            triangle_xy=xy)
        print(f"n={n:3d}: elements/grain {counts.min()}-{counts.max()} "
              f"(median {np.median(counts):.0f}), empty {int((counts==0).sum())}, "
              f"width {np.median(width):.2f} um, aspect {np.median(aspect):.2f}:1, "
              f"IPF coverage maximin {rec['ipf_coverage_maximin']*100:.1f}% vs "
              f"random {rec['ipf_coverage_pure_random']*100:.1f}%")

    (out / f"polycrystal_sets_d{args.density}.json").write_text(json.dumps(records, indent=1))
    print(f"\nelement edge {h_el:.3f} um, thickness {thick:.3f} um "
          f"-> elements are {'near-cubic' if abs(h_el/thick-1)<0.25 else 'NOT cubic'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
