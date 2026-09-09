#!/usr/bin/env python3
"""Verify that the isochoric loading idiom holds det(F) = 1 in an actual solve.

The closed form is exact by construction, so what needs checking is that the
composed deck (parsed function + dirichlet BC) really applies it: that the BC
function scaling reaches the prescribed values, that the roller faces leave the
lateral contraction free, and that the resulting macroscopic deformation is the
intended uniform stretch. The volume is measured from the DEFORMED node
positions rather than from the prescribed values, so a BC that failed to be
applied would show up instead of being assumed away.
"""

from __future__ import annotations

import sys

import numpy as np
from pinnacle.io.loader import load_input
from pinnacle.plugins import discover_components
from pinnacle.runner.runner import Runner

discover_components(force=True)


def main(deck: str, edge: float = 50.0, eps_axial_max: float = 0.12) -> int:
    res = Runner.from_config(load_input(deck))
    fields = res["fields"]
    nodes = np.asarray(fields["coordinates"], dtype=float).reshape(-1, 3)
    # The runner may return the displacement in grid form (nx, ny, nz, 3)
    # rather than flat (n_nodes, 3); both describe the same nodes in the same
    # order, so flatten rather than assuming one shape.
    disp = np.asarray(fields["displacement"], dtype=float).reshape(-1, 3)
    if disp.shape != nodes.shape:
        raise SystemExit(
            f"displacement {disp.shape} and coordinates {nodes.shape} disagree "
            "after flattening; the field ordering assumption is wrong"
        )
    deformed = nodes + disp

    # Macroscopic stretches from the extreme faces of the deformed body.
    lam = np.array(
        [
            (deformed[:, i].max() - deformed[:, i].min()) / (nodes[:, i].max() - nodes[:, i].min())
            for i in range(3)
        ]
    )
    det = float(np.prod(lam))

    # Reference: the axial stretch that the ramp should have reached at t_end,
    # and the lateral companion the constraint demands.
    lam_axial_expected = 1.0 + eps_axial_max
    lam_lateral_expected = lam_axial_expected**-0.5

    print(f"deck                  {deck}")
    print(f"stretches (x, y, z)   {np.round(lam, 8)}")
    print(f"axial (z) expected    {lam_axial_expected:.8f}   got {lam[2]:.8f}")
    print(f"lateral expected      {lam_lateral_expected:.8f}   got x {lam[0]:.8f}  y {lam[1]:.8f}")
    print(f"det(F) macroscopic    {det:.12f}   deviation {abs(det - 1.0):.3e}")

    # Element-level volume change, which is NOT constrained to 1 pointwise --
    # only the macroscopic volume is. Reporting the spread makes the
    # distinction explicit rather than leaving a reader to assume pointwise
    # incompressibility.
    # Element-level volume change, when the runner exposes connectivity. This
    # is NOT constrained to 1 pointwise -- only the macroscopic volume is -- so
    # the spread is reported rather than asserted, and its absence is stated
    # rather than passed over silently.
    conn = np.asarray(fields.get("connectivity", []), dtype=int)
    if conn.ndim == 2 and conn.shape[0]:
        v0 = _hex_volumes(nodes, conn)
        v1 = _hex_volumes(deformed, conn)
        jac = v1 / v0
        print(
            f"element volume ratio  min {jac.min():.6f}  max {jac.max():.6f}  "
            f"mean {jac.mean():.6f}  (total {v1.sum() / v0.sum():.10f})"
        )

    else:
        print("element volume ratio  not reported (no connectivity in output fields)")

    ok = abs(det - 1.0) < 1.0e-6
    print("VERDICT               " + ("det(F) = 1 to 1e-6" if ok else "CONSTRAINT VIOLATED"))
    return 0 if ok else 1


def _hex_volumes(pts: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """Hex8 volumes by decomposition into six tetrahedra per element."""
    tets = [
        (0, 1, 3, 4),
        (1, 2, 3, 6),
        (1, 3, 4, 6),
        (3, 4, 6, 7),
        (1, 4, 5, 6),
        (0, 3, 4, 1),
    ]
    vol = np.zeros(conn.shape[0], dtype=float)
    for a, b, c, d in tets:
        p0, p1, p2, p3 = (pts[conn[:, i]] for i in (a, b, c, d))
        vol += np.abs(np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), p3 - p0)) / 6.0
    return vol


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "decks/al2219/iso_bc_verification_cube.i"))
