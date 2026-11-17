#!/usr/bin/env python3
"""Guess scan of case A3: from which initial droplets does each route reach the nucleus?

For every route (fem, pinn, hybrid) the standalone deck
`decks/crystallization/single_field/a3_ni_nucleus_x0.200_<route>.i` is copied with its
initial droplet replaced by a full-amplitude tanh droplet (planar width) whose radius is the
classical radius times a factor, run, and the outcome recorded: W* against the reference, or
the failure (Newton did not converge, collapse onto the liquid, or another state).

Usage: python3 scripts/cryst_guess_scan.py [factors...]
Output: results/crystallization/single_field/guess_scan/guess_scan.csv
"""
from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECKS = ROOT / "decks/crystallization/single_field"
OUT = ROOT / "results/crystallization/single_field/guess_scan"
R_CLASSICAL = 0.9814  # nm, 2 gamma(T) / |Delta_f| at Delta T / T_f = 0.2
W_REF = 0.9137052  # aJ, shooting solution
factors = [float(a) for a in sys.argv[1:]] or [0.4, 0.55, 0.7, 0.85, 1.0, 1.15, 1.3, 1.6, 2.0]

rows = []
OUT.mkdir(parents=True, exist_ok=True)
for route in ("fem", "pinn", "hybrid"):
    base = (DECKS / f"a3_ni_nucleus_x0.200_{route}.i").read_text()
    for f in factors:
        name = f"{route}_r{f:.2f}"
        text = re.sub(r"expression = [^\n]*",
                      f"expression = 0.5*(1-tanh((x-{R_CLASSICAL * f:.4f})/(2*0.22756)))", base, count=1)
        text = text.replace(f"results/crystallization/single_field/a3_ni_nucleus_x0.200_{route}",
                            f"results/crystallization/single_field/guess_scan/{name}")
        deck = OUT / f"{name}.i"
        deck.write_text(text)
        t0 = time.time()
        env = dict(os.environ, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4")
        proc = subprocess.run(["python3.11", "-m", "pinnacle", "run", str(deck)], cwd=ROOT, env=env,
                              capture_output=True, text=True)
        wall = time.time() - t0
        post = OUT / name / "postprocessors.csv"
        w_star = float("nan")
        status = "ok"
        if proc.returncode != 0:
            last = (proc.stdout + proc.stderr).strip().splitlines()[-1]
            status = ("newton failed" if "Newton did not converge" in last
                      else "collapsed" if "trivial stationary state" in last else "failed: " + last[:80])
        elif post.exists():
            w_star = float(post.read_text().splitlines()[-1])
            if abs(w_star / W_REF - 1.0) > 1e-2:
                status = "other state"
        rows.append({"route": route, "radius_factor": f, "guess_radius_nm": round(R_CLASSICAL * f, 4),
                     "W_star_aJ": w_star, "rel_error": w_star / W_REF - 1.0, "status": status,
                     "wall_s": round(wall, 1)})
        print(rows[-1], flush=True)
with (OUT / "guess_scan.csv").open("w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    wr.writeheader()
    wr.writerows(rows)
