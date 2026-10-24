#!/usr/bin/env python3
"""Experiment-only crops of the Fig. 5 panels (data/digitized/ni_fig5_panels/<spec>.png).

Each panel stacks the DIC experiment ("Exp.") above the paper's own simulation ("Sim.").
The author (2026-10-23) wants only the experiment next to this study's fields, so the
top ``FRACTION`` of every panel (the Exp. strip, the specimen label and the color scale,
which sits between the two strips) is written as ``<spec>_exp.png``.

Usage: ni_fig5_exp_crops.py [FRACTION]   (default 0.62)
"""
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent / "data/digitized/ni_fig5_panels"


def main(argv):
    frac = float(argv[1]) if len(argv) > 1 else 0.62
    for p in sorted(ROOT.glob("*.png")):
        if p.stem.endswith("_exp"):
            continue
        im = Image.open(p)
        im.crop((0, 0, im.width, int(round(frac * im.height)))).save(p.with_name(p.stem + "_exp.png"))
        print("wrote", p.with_name(p.stem + "_exp.png").name)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
