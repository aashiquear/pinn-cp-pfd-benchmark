#!/usr/bin/env python3
"""Collect the PINNACLE critical-nucleus runs (cases A3, A4) against the reference table.

Usage: python3 scripts/cryst_nucleus_collect.py [route]     (route: fem | pinn | hybrid)
Reads results/crystallization/single_field/a{3,4}_ni_nucleus_*_<route>/ and writes
results/crystallization/single_field/nucleus_<route>_table.csv.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from vtk_legacy import read_vtk  # noqa: E402

K_B = 1.380649e-23
route = sys.argv[1] if len(sys.argv) > 1 else "fem"
base = ROOT / "results/crystallization/single_field"
ref = {(r["well"], r["interpolation"], round(float(r["reduced_undercooling"]), 3)): r
       for r in csv.DictReader((ROOT / "results/crystallization/reference/ni_nucleus_3d/nucleus_table.csv").open())
       if r["temperature_scaling"] == "proportional"}
rows = []
for run in sorted(base.glob(f"a[34]_ni_nucleus_*_{route}")):
    m = re.match(rf"a[34]_ni_nucleus_(?:(\w+?)_(standard|gl_bcc|gl_fcc)_)?x([0-9.]+)_{route}", run.name)
    well, interp = (m.group(1) or "quartic"), (m.group(2) or "standard")
    x = round(float(m.group(3)), 3)
    r = ref[(well, interp, x)]
    if not (run / "postprocessors.csv").exists():
        print("missing or failed:", run.name)
        continue
    hist_rows = list(csv.DictReader((run / "solver_history.csv").open()))
    # FE and hybrid: the corrector's Newton row; PINN: the training history (no such row)
    newton = [h for h in hist_rows if h.get("newton_iters") not in (None, "")]
    hist = newton[-1] if newton else {}
    mesh = read_vtk(run / "fields.vtk")
    radius = mesh["nodes"][:, 0]
    phi = mesh["point"]["phi"]
    order = np.argsort(radius)
    radius, phi = radius[order], phi[order]
    r_e = (3.0 * np.trapezoid(radius**2 * phi, radius) / phi[0]) ** (1.0 / 3.0)
    kt = K_B * float(r["T_K"]) * 1e18
    w_star = float((run / "postprocessors.csv").read_text().splitlines()[-1])
    w_ref = float(r["W_star_J"]) * 1e18
    rows.append({
        "well": well, "interpolation": interp, "reduced_undercooling": x,
        "W_star_aJ": w_star, "W_star_over_kT": w_star / kt, "W_ref_over_kT": w_ref / kt,
        "W_rel_error": w_star / w_ref - 1.0,
        "phi_center": float(phi[0]), "phi_center_ref": float(r["phi_center"]),
        "R_equimolar_nm": r_e, "R_equimolar_ref_nm": float(r["R_equimolar_m"]) * 1e9,
        "newton_iters": int(float(hist.get("newton_iters", 0) or 0)),
        "converged": int(float(hist.get("converged", 1) or 1)),
    })
out = base / f"nucleus_{route}_table.csv"
with out.open("w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    wr.writeheader()
    wr.writerows(rows)
for r in rows:
    print(f"{r['well']:8s} {r['interpolation']:9s} x={r['reduced_undercooling']:.3f} "
          f"W*/kT {r['W_star_over_kT']:8.3f} ref {r['W_ref_over_kT']:8.3f} ({r['W_rel_error']:+.2e}) "
          f"phi0 {r['phi_center']:.5f}/{r['phi_center_ref']:.5f} R_e {r['R_equimolar_nm']:.4f}/"
          f"{r['R_equimolar_ref_nm']:.4f} it {r['newton_iters']}")
