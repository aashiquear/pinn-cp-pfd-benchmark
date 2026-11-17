#!/usr/bin/env python3
"""Digitize the two barrier curves of the bottom panel of Fig. 6 of (granasy 2019 nucleation).

Input: the raster of the review's Fig. 6, `pdfimages -png -f 11 -l 11
data/literature/crystallization_review.pdf results/crystallization/digitize/p11`
(869 x 777 px, 300 ppi; not shipped). The bottom panel shows W / W_cap against the undercooling
for Model I (solid) and Model II (dashed) at a contact angle of 61.2 deg.
Method: the axes box is the pair of longest dark rows and columns of the panel; inside it,
every pixel column is scanned for dark runs. The upper run is the solid line; a second,
lower run is a dash of the dashed line. Axis ranges from the tick labels: x 0..100 K,
y 0..1. Output: data/digitized/cryst_review_fig6_barriers.csv (curve, undercooling, ratio);
the half thickness of the stroke is the extraction uncertainty.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/crystallization/digitize/p11-000.png"
img = np.asarray(Image.open(src).convert("L")).astype(float)
top = 430
panel = img[top:, :]
dark = panel < 120
rows = np.where(dark.sum(axis=1) > 0.8 * dark.sum(axis=1).max())[0]
cols = np.where(dark.sum(axis=0) > 0.6 * dark.sum(axis=0).max())[0]
# Long dark rows of the crop: the lower border of the panel above (first rows), then the
# top and the bottom frame of this panel. The frame is the last two groups.
groups, cur = [], [rows[0]]
for r in rows[1:]:
    if r - cur[-1] <= 2:
        cur.append(r)
    else:
        groups.append(float(np.mean(cur)))
        cur = [r]
groups.append(float(np.mean(cur)))
y_top, y_bottom = int(round(groups[-2])), int(round(groups[-1]))
x_left, x_right = cols.min(), cols.max()


def to_x(px):
    return (px - x_left) / (x_right - x_left) * 100.0


def to_y(py):
    return (y_bottom - py) / (y_bottom - y_top) * 1.0


out = []
MARGIN = 8  # px: tick marks reach this far into the box (bottom)
TOP = 3  # px below the top frame; a run touching it is a tick mark or the frame
for px in range(x_left + 3, x_right - 2):
    if to_x(px) < 10.0:
        continue  # both curves start at 10 K
    col = dark[y_top + TOP:y_bottom - MARGIN, px]
    runs, start = [], None
    for i, v in enumerate(col):
        if v and start is None:
            start = i
        if (not v) and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(col) - 1))
    runs = [r for r in runs if r[1] - r[0] + 1 <= 8 and r[0] > 0]
    if not runs:
        continue
    # tick columns (every 20 K) carry a tick mark that merges with the curves: skip them
    if min(abs(to_x(px) - t) for t in (0, 20, 40, 60, 80, 100)) < 0.5:
        continue
    for a, b in (runs[0], runs[-1]) if len(runs) > 1 else (runs[0],):
        ratio = to_y(y_top + TOP + 0.5 * (a + b))
        # the two strokes are separated by the line 0.96 - 0.004 (dT - 10 K): the solid one
        # (Model I) runs from 1.00 to 0.80, the dashed one (Model II) from 0.92 to 0
        curve = "model_I" if ratio > 0.96 - 0.004 * (to_x(px) - 10.0) else "model_II"
        out.append((curve, to_x(px), ratio, 0.5 * (b - a + 1) / (y_bottom - y_top)))
dest = ROOT / "data/digitized/cryst_review_fig6_barriers.csv"
with dest.open("w") as fh:
    fh.write("# W / W_cap of Fig. 6 (bottom), Granasy et al. 2019; provenance: digitized\n")
    fh.write(f"# axes box px: x {x_left}..{x_right} = 0..100 K; y {top + y_bottom}..{top + y_top} = 0..1; "
             f"{100.0 / (x_right - x_left):.3f} K and {1.0 / (y_bottom - y_top):.4f} per px\n")
    fh.write("curve,undercooling_K,ratio,half_stroke\n")
    for c, x, y, e in out:
        fh.write(f"{c},{x:.3f},{y:.4f},{e:.4f}\n")
print("box", x_left, x_right, y_top, y_bottom, "points", len(out))
for c in ("model_I", "model_II"):
    xs = np.array([o[1] for o in out if o[0] == c]); ys = np.array([o[2] for o in out if o[0] == c])
    print(c, len(xs), "range", xs.min().round(1), xs.max().round(1),
          [(t, round(float(np.interp(t, xs, ys)), 3)) for t in (10, 20, 40, 60, 80, 100) if xs.min() <= t <= xs.max()])
