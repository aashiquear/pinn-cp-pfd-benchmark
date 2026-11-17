#!/usr/bin/env python3
"""Digitize the 'standard PF' curve of Fig. 3(c) of (granasy 2019 nucleation).

Input: the raster of the review's Fig. 3 extracted with
`pdfimages -png -f 9 -l 9 data/literature/crystallization_review.pdf <prefix>`
(1371 x 1256 px, 300 ppi; not shipped, the paper is not redistributed).
Method: the axes box of panel (c) is located from its longest dark rows and columns; the
'standard PF' line is the only stroke thicker than 6 px, so in every pixel column the vertical
dark run of 7 to 30 px is taken and its center mapped linearly to (Delta T / T_f, W / kT).
Output: data/digitized/cryst_review_fig3c_standard_pf.csv with the per-column value and the
half-thickness of the stroke as the extraction uncertainty.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/crystallization/digitize/rev-000.png"
img = np.asarray(Image.open(src).convert("L")).astype(float)
panel = img[620:, :700]  # panel (c): lower-left quadrant
dark = panel < 110

# Long dark columns / rows of the panel, grouped: columns = [left frame, legend left edge,
# right frame, legend right edge]; rows = [legend top, top frame, legend bottom, bottom frame]
# (the legend overlaps the upper-right corner of the axes box and sticks out of it).
def groups(lengths, minimum):
    idx = np.where(lengths > minimum)[0]
    out, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 2:
            cur.append(i)
        else:
            out.append(float(np.mean(cur)))
            cur = [i]
    out.append(float(np.mean(cur)))
    return out


cols = groups(dark.sum(axis=0), 200)
rws = groups(dark.sum(axis=1), 250)
assert len(cols) == 4 and len(rws) == 4, (cols, rws)
x_left, x_right = cols[0], cols[2]
y_top, y_bottom = rws[1], rws[3]
# axis ranges read from the tick labels: x 0.1 .. 0.3, y 0 .. 120
X0, X1, Y0, Y1 = 0.1, 0.3, 0.0, 120.0


def to_x(px):
    return X0 + (px - x_left) / (x_right - x_left) * (X1 - X0)


def to_y(py):
    return Y0 + (y_bottom - py) / (y_bottom - y_top) * (Y1 - Y0)


rows = []
for px in range(int(x_left) + 4, int(x_right) - 3):
    y_a = int(y_top) + 4
    col = dark[y_a:int(y_bottom) - 3, px]
    runs = []
    start = None
    for i, v in enumerate(col):
        if v and start is None:
            start = i
        if (not v) and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(col) - 1))
    # strokes inside the legend rectangle (columns cols[1]..cols[3], rows rws[0]..rws[2]) are not data
    in_legend = cols[1] - 2 <= px <= cols[3] + 2
    thick = [(a, b) for a, b in runs if 7 <= b - a + 1 <= 30
             and not (in_legend and y_a + a <= rws[2] + 3)]
    if len(thick) != 1:
        continue  # legend, markers or a merged stroke: skip the column
    a, b = thick[0]
    yc = y_a + 0.5 * (a + b)
    rows.append((to_x(px), to_y(yc), 0.5 * (b - a + 1) / (y_bottom - y_top) * (Y1 - Y0)))

out = ROOT / "data/digitized/cryst_review_fig3c_standard_pf.csv"
with out.open("w", newline="") as fh:
    fh.write("# standard PF line of Fig. 3(c), Granasy et al. 2019; provenance: digitized\n")
    fh.write(f"# axes box px: x {x_left}..{x_right} = {X0}..{X1}; y {y_bottom}..{y_top} = {Y0}..{Y1}; "
             f"resolution {(X1 - X0) / (x_right - x_left):.2e} per px in x, "
             f"{(Y1 - Y0) / (y_bottom - y_top):.3f} kT per px in y\n")
    wr = csv.writer(fh)
    wr.writerow(["reduced_undercooling", "W_over_kT", "half_stroke_kT"])
    for r in rows:
        wr.writerow([f"{r[0]:.5f}", f"{r[1]:.3f}", f"{r[2]:.3f}"])
print("box", x_left, x_right, y_top, y_bottom, "points", len(rows))
xs = np.array([r[0] for r in rows]); ys = np.array([r[1] for r in rows])
for x in (0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3):
    i = np.argmin(abs(xs - x))
    print(f"{x:.3f} -> x {xs[i]:.4f}  W/kT {ys[i]:.2f}")
