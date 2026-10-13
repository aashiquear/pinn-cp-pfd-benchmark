#!/usr/bin/env python3
"""Write N random parameter rows of a design deck's [Parameters] space (for `evaluate`).

Usage: surrogate_queries.py decks/cp_study/surrogate/design_h40.i 1000 results/cp_study/surrogate/queries_1000.csv [seed]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from surrogate_design import read_deck  # noqa: E402


def main(argv):
    if len(argv) < 4:
        print(__doc__); return 2
    from pinnacle.pinn.parametric import ParameterSpace
    cfg = read_deck(Path(argv[1]).resolve()); n = int(argv[2]); out = Path(argv[3]); seed = int(argv[4]) if len(argv) > 4 else 1
    space = ParameterSpace(dict(cfg["parameters"])); rows = space.sample(n, "random", seed)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh); w.writerow(space.names); w.writerows([[f"{v:.12g}" for v in r] for r in rows])
    print("wrote", out, n, "rows"); return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
