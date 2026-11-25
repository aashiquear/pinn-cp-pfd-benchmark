#!/usr/bin/env python3
"""Part V: Voronoi grains with uniform random orientations as a PINNACLE orientation point file.

usage: scripts/lpbf_orientations.py <out.csv> <xmax> <ymax> <zmax> <grain size> <seed>

Writes one row per grain seed, `x,y,z,phi1,Phi,phi2,grain_id` (Bunge angles in degrees).
With `[OrientationField] type = point_file, interpolation = nearest` every quadrature
point takes the orientation of the nearest seed: a Voronoi polycrystal. The number of
grains is the box volume over the cube of the grain size. Orientations are uniform over
the rotation group (random unit quaternions), the choice of decision D4.
"""

from __future__ import annotations

import sys

import numpy as np


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    out = sys.argv[1]
    box = np.array([float(v) for v in sys.argv[2:5]])
    size, seed = float(sys.argv[5]), int(sys.argv[6])
    count = max(int(round(np.prod(box) / size**3)), 1)
    rng = np.random.default_rng(seed)
    points = rng.random((count, 3)) * box
    q = rng.normal(size=(count, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q.T
    # rotation matrix of the quaternion, then its Bunge angles (R = Rz(phi1) Rx(Phi) Rz(phi2))
    r22 = 1 - 2 * (x * x + y * y)
    r20, r21 = 2 * (x * z - y * w), 2 * (y * z + x * w)
    r02, r12 = 2 * (x * z + y * w), 2 * (y * z - x * w)
    big = np.degrees(np.arccos(np.clip(r22, -1.0, 1.0)))
    phi1 = np.degrees(np.arctan2(r02, -r12)) % 360.0
    phi2 = np.degrees(np.arctan2(r20, r21)) % 360.0
    with open(out, "w") as fh:
        fh.write(f"# {count} Voronoi grain seeds, uniform random orientations, seed {seed}, "
                 f"box {box.tolist()}, grain size {size}\n")
        fh.write("x,y,z,phi1,Phi,phi2,grain_id\n")
        for i in range(count):
            fh.write(f"{points[i,0]:.9e},{points[i,1]:.9e},{points[i,2]:.9e},"
                     f"{phi1[i]:.6f},{big[i]:.6f},{phi2[i]:.6f},{i}\n")
    print(f"wrote {out}: {count} grains")


if __name__ == "__main__":
    main()
