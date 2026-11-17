#!/usr/bin/env python3
"""Collect case B1 (heterogeneous nucleus on a wall) against the digitized Fig. 6 curves.

W_I = 2 (free_energy + wall_energy); W_II = 2 (F_nucleus - F_layer). Both are divided by the
sharp-interface circular cap W_cap = (gamma^2/|Delta_f|)(Psi - sin Psi cos Psi). The papers
do not say which gamma the cap uses: `ratio_gammaT` takes gamma(T) = gamma T/T_m (the model's
own interface energy), `ratio_gamma0` the melting-point value (cap larger by (T_m/T)^2).

Usage: python3 scripts/cryst_wall_collect.py
Output: results/crystallization/walls/wall_table.csv
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/walls"
PSI = math.radians(61.2)
T_M = 1728.0
L_OVER_W = 2.666e9 / (12.0 * 0.364 / (2.0e-9 / (2.0 * math.log(9.0))))
dig = [ln.strip().split(",") for ln in (ROOT / "data/digitized/cryst_review_fig6_barriers.csv").open()
       if ln.startswith("model")]
pub = {c: (np.array([float(r[1]) for r in dig if r[0] == c]), np.array([float(r[2]) for r in dig if r[0] == c]))
       for c in ("model_I", "model_II")}


def post(name: str) -> dict | None:
    f = RES / name / "postprocessors.csv"
    if not f.exists():
        return None
    head, vals = f.read_text().split("\n")[:2]
    return dict(zip(head.split(","), map(float, vals.split(","))))


def saddle(name: str) -> str:
    log = RES / name / "pinnacle.log"
    m = re.search(r"pf_fem_saddle: (\d+) volume iterations \((\d+) constrained", log.read_text()) if log.exists() else None
    return f"{m.group(1)}/{m.group(2)}" if m else ""


rows = []
for d_t in range(10, 101, 10):
    big_d = -L_OVER_W * d_t / (T_M - d_t)
    w_cap = (1.0 / 12.0) ** 2 / abs(big_d) * (PSI - math.sin(PSI) * math.cos(PSI))
    scale0 = ((T_M - d_t) / T_M) ** 2
    for model in ("model_I", "model_II"):
        if model == "model_I":
            a = post(f"b1_wall_model1_dT{d_t:03d}_fem")
            w = 2.0 * (a["free_energy"] + a["wall_energy"]) if a else None
            its = saddle(f"b1_wall_model1_dT{d_t:03d}_fem")
            angle = a["angle"] if a else float("nan")
        else:
            n, lay = post(f"b1_wall_model2_dT{d_t:03d}_nucleus_fem"), post(f"b1_wall_model2_dT{d_t:03d}_layer_fem")
            w = 2.0 * (n["free_energy"] - lay["free_energy"]) if n and lay else None
            its = saddle(f"b1_wall_model2_dT{d_t:03d}_nucleus_fem")
            angle = float("nan")
        if w is None:
            continue
        xs, ys = pub[model]
        p = float(np.interp(d_t, xs, ys)) if xs.min() <= d_t <= xs.max() else float("nan")
        rows.append({"model": model, "undercooling_K": d_t, "D": big_d, "W": w, "W_cap": w_cap,
                     "ratio_gammaT": w / w_cap, "ratio_gamma0": w / w_cap * scale0, "published": p,
                     "wall_angle_deg": angle, "volume_iters/newton_iters": its})
with (RES / "wall_table.csv").open("w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    wr.writeheader()
    wr.writerows(rows)
for r in rows:
    print(f"{r['model']:8s} {r['undercooling_K']:3d} K  W/W_cap {r['ratio_gammaT']:.4f} (gamma(T))  "
          f"{r['ratio_gamma0']:.4f} (gamma(T_m))  published {r['published']:.3f}  angle {r['wall_angle_deg']:.1f}  {r['volume_iters/newton_iters']}")
