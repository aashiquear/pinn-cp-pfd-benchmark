#!/usr/bin/env python3
"""Write the mesh of case E5: a liquid layer on a rough (sawtooth) wall, Gmsh 2.2 ASCII quads.

The wall is the bottom boundary y = a * tri(x / P) (a triangular wave of amplitude a and period
P, grooves and ridges with straight flanks); the top is flat. The grid is a mapped structured
grid of nx x ny bilinear quadrilaterals: the columns stay vertical and the rows follow the wall
near the bottom and are horizontal at the top. All boundaries are natural (no flux), which for
the phase field is a wall with a contact angle of 90 degrees.

Usage: python3 scripts/cryst_rough_wall_mesh.py [nx=400] [ny=200] [dx=6.25e-3] [amplitude_cells=30]
           [period_cells=80] [out=decks/crystallization/walls_rough/e5_rough_wall.msh]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    opts = dict(a.split("=") for a in sys.argv[1:])
    nx, ny = int(opts.get("nx", 400)), int(opts.get("ny", 200))
    dx = float(opts.get("dx", 6.25e-3))
    amp, period = float(opts.get("amplitude_cells", 30)) * dx, float(opts.get("period_cells", 80)) * dx
    out = ROOT / opts.get("out", "decks/crystallization/walls_rough/e5_rough_wall.msh")
    x = np.arange(nx + 1) * dx
    phase = (x / period) % 1.0
    wall = amp * (1.0 - 2.0 * np.abs(2.0 * phase - 1.0)) * 0.5 + 0.5 * amp  # 0 (groove) .. amp (ridge)
    top = ny * dx + 0.5 * amp
    eta = np.linspace(0.0, 1.0, ny + 1)
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$Nodes", str((nx + 1) * (ny + 1))]
    for j in range(ny + 1):
        for i in range(nx + 1):
            y = wall[i] + eta[j] * (top - wall[i])
            lines.append(f"{j * (nx + 1) + i + 1} {x[i]:.10g} {y:.10g} 0")
    lines += ["$EndNodes", "$Elements", str(nx * ny)]
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i + 1
            lines.append(f"{j * nx + i + 1} 3 2 1 1 {n0} {n0 + 1} {n0 + nx + 2} {n0 + nx + 1}")
    lines.append("$EndElements")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print(out.relative_to(ROOT), f"{(nx + 1) * (ny + 1)} nodes, {nx * ny} quads; wall amplitude {amp:.4g}, period {period:.4g}")


if __name__ == "__main__":
    main()
