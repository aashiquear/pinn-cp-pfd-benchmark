#!/usr/bin/env python3
"""Convert a PINNACLE binary-model run (stage C) to the frames.npz / history.csv layout.

Usage: python3 scripts/cryst_binary_collect.py <run directory under results/crystallization/binary>
       python3 scripts/cryst_binary_collect.py orientation/<run>      (also stores theta)
Reads the fields_step*.vtk series (phi, c on a structured box mesh) and fields.pvd; writes
frames.npz (phi, c, steps, dx, dt, c_inf, c_l, c_s) and history.csv with the same columns as
scripts/cryst_binary_reference.py, so scripts/cryst_binary_figures.py serves both.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from vtk_legacy import read_vtk, series_files  # noqa: E402

C_L, C_S = 0.466219, 0.399112


def main() -> None:
    arg = sys.argv[1]
    run = ROOT / "results/crystallization" / (arg if "/" in arg else f"binary/{arg}")
    times = [float(t) for t in re.findall(r'timestep="([^"]+)"', (run / "fields.pvd").read_text())]
    files = series_files(run)
    first = read_vtk(files[0][1])
    x, y = first["nodes"][:, 0], first["nodes"][:, 1]
    xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
    ix, iy = np.searchsorted(xs, np.round(x, 9)), np.searchsorted(ys, np.round(y, 9))
    dx = float(xs[1] - xs[0])
    with (run / "solver_history.csv").open() as fh:
        dt = float(next(csv.DictReader(fh))["dt"])
    phis, cs, thetas, rows = [], [], [], []
    solute0 = None
    for (step, path), t in zip(files, times):
        data = read_vtk(path)["point"]
        phi = np.zeros((len(ys), len(xs)))
        c = np.zeros_like(phi)
        phi[iy, ix] = data["phi"]
        c[iy, ix] = data["c"]
        # drop the duplicated periodic row and column: cells = nodes - 1
        phi, c = phi[:-1, :-1], c[:-1, :-1]
        if "theta" in data:
            theta = np.zeros((len(ys), len(xs)))
            theta[iy, ix] = data["theta"]
            thetas.append(theta[:-1, :-1].astype(np.float32))
        phis.append(phi.astype(np.float32))
        cs.append(c.astype(np.float32))
        solute0 = c.sum() if solute0 is None else solute0
        solid = float(np.mean(1.0 - phi))
        step = int(round(t / dt))  # the series is numbered by frame, not by step
        rows.append({"step": step, "t": t, "solid_fraction": solid,
                     "equivalent_radius_dx": float(np.sqrt(solid * phi.size / np.pi)),
                     "solute_drift": float(c.sum() / solute0 - 1.0),
                     "phi_min": float(phi.min()), "phi_max": float(phi.max()),
                     "c_min": float(c.min()), "c_max": float(c.max())})
    c_inf = float(np.median(cs[0]))
    extra = {"theta": np.array(thetas)} if thetas else {}
    np.savez_compressed(run / "frames.npz", phi=np.array(phis), c=np.array(cs), **extra,
                        steps=np.array([r["step"] for r in rows]), dx=dx, dt=dt, c_inf=c_inf, c_l=C_L, c_s=C_S)
    with (run / "history.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows[1:])
    print(run.name, "frames", len(files), "last:", rows[-1])


if __name__ == "__main__":
    main()
