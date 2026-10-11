#!/usr/bin/env python3
"""Write a grain-growth deck for a NEW random Voronoi seeding of the plate.

The study's grain structure is a Voronoi seeding evolved by phase-field grain
growth (`decks/graingrowth/pilot/n32_bdf1_h40.i`). A second or third
realization for the grain-instance check keeps every numerical setting of
that deck and changes only the seed points: ``n`` points drawn uniformly in
the plate (a Poisson seeding, the same process the original 32 seeds came
from) with ``numpy.random.default_rng(seed)``. The drawn coordinates are
written INTO the new deck, so the structure is reproducible from the deck
alone and the seed is recorded in its header.

Usage: make_growth_instance.py TEMPLATE.i --seed INT --out NEW.i --output-dir DIR [--n 32]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("template")
    ap.add_argument("--seed", type=int, required=True, help="rng seed of the Voronoi points")
    ap.add_argument("--n", type=int, default=32, help="number of Voronoi points")
    ap.add_argument("--out", required=True, help="deck to write")
    ap.add_argument("--output-dir", required=True, help="[Global] output_dir of the new deck")
    a = ap.parse_args()
    text = Path(a.template).read_text()
    m = re.search(r"domain\s*=\s*'([^']+)'", text)
    if not m:
        raise SystemExit("template has no [Mesh] domain")
    x0, y0, x1, y1 = (float(v) for v in m.group(1).split())
    rng = np.random.default_rng(a.seed)
    pts = np.column_stack([rng.uniform(x0, x1, a.n), rng.uniform(y0, y1, a.n)])
    flat = " ".join(f"{v:.6f}" for v in pts.reshape(-1))
    text, n1 = re.subn(r"seeds\s*=\s*'[^']*'", f"seeds          = '{flat}'", text)
    text, n2 = re.subn(r"n_grains\s*=\s*\d+", f"n_grains       = {a.n}", text)
    text, n3 = re.subn(r"output_dir\s*=\s*\S+", f"output_dir       = {a.output_dir}", text)
    if not (n1 == n2 == n3 == 1):
        raise SystemExit(f"template substitution failed (seeds {n1}, n_grains {n2}, output_dir {n3})")
    header = (f"# GRAIN INSTANCE: {a.n} Voronoi points drawn uniformly in the plate with\n"
              f"# numpy.random.default_rng({a.seed}) by scripts/make_growth_instance.py from\n"
              f"# {Path(a.template).as_posix()}; every other setting is that deck's.\n")
    Path(a.out).write_text(header + text)
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
