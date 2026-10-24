#!/usr/bin/env python3
"""Experiment-only crops of the Fig. 5 panels (data/digitized/ni_fig5_panels/<spec>.png).

Each panel stacks the DIC experiment ("Exp.") above the paper's own simulation ("Sim."),
with one e_xx color bar between the two strips at the right. The author (2026-10-23)
wants only the experiment next to this study's fields, so for every panel the Exp. strip
is kept (its rows found from the colored pixels), the specimen label below it is cleared,
and the color-bar block (the bar with its labels, found by its saturated columns) is moved
up beside the strip, so that the scale stays complete. Written as ``<spec>_exp.png``.

Usage: ni_fig5_exp_crops.py
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent / "data/digitized/ni_fig5_panels"
SAT = 120  # a pixel with max-min RGB above this is colored (the DIC field, the color bar)


def _saturation(a: np.ndarray) -> np.ndarray:
    return a.max(axis=2).astype(int) - a.min(axis=2).astype(int)


def layout(a: np.ndarray) -> dict:
    """Rows of the Exp. strip and the box of the color-bar block, per panel.

    The strip is the saturated rows of the central columns in the upper half; the bar is
    the rightmost narrow run of saturated columns in the rows below the strip (the crack's
    red band would otherwise pass for it). The block takes the bar with its labels (14 px
    above, 18 px below, text to the right edge).
    """
    sat = _saturation(a)
    h, w = sat.shape
    mid = sat[: h // 2, w // 8 : w - 80]
    srows = np.where((mid > 80).sum(axis=1) > 20)[0]
    s0, s1 = int(srows.min()), int(srows.max())
    below = sat[s1 + 5 :, w // 2 :]
    cols = np.where((below > SAT).sum(axis=0) > 40)[0] + w // 2
    bx1 = int(cols.max())
    bx0 = bx1
    while bx0 - 1 in set(cols.tolist()):
        bx0 -= 1
    rows = np.where((sat[s1 + 5 :, bx0 : bx1 + 1] > SAT).sum(axis=1) > 3)[0] + s1 + 5
    by0, by1 = int(rows.min()), int(rows.max())
    return {"strip": (s0, s1), "bar": (bx0 - 10, by0 - 14, w, by1 + 18)}


def main(argv):
    for p in sorted(ROOT.glob("*.png")):
        if p.stem.endswith("_exp"):
            continue
        im = Image.open(p).convert("RGB")
        a = np.asarray(im)
        L = layout(a)
        s0, s1 = L["strip"]
        x0, y0, x1, y1 = L["bar"]
        block = im.crop((x0, y0, x1, y1))
        center = (s0 + s1) // 2
        top = max(s1 + 14, center + block.height // 2 + 2)
        out = im.crop((0, 0, im.width, top))
        # below the strip only the grips' speckle may stay: clear the label zone and the bar column
        out.paste((0, 0, 0), (im.width // 5, s1 + 8, x0, top))
        out.paste((0, 0, 0), (x0, 0, im.width, top))
        out.paste(block, (x0, center - block.height // 2))
        out.save(p.with_name(p.stem + "_exp.png"))
        print("wrote", p.with_name(p.stem + "_exp.png").name, "strip rows", L["strip"], "bar", L["bar"], "height", top)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
