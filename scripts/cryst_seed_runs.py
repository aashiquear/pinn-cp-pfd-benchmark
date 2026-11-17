#!/usr/bin/env python3
"""Seed statistics of the network routes of stage A (PROJECT_GUIDELINES section 7).

Each listed deck is run with `[Global] seed` = 0..n-1 (the deck itself carries seed 0); the
free energy of every run is collected, and mean, standard deviation and worst case against
the reference are written to results/crystallization/single_field/seeds/seed_table.csv.

Usage: python3 scripts/cryst_seed_runs.py [n_seeds]
"""
from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DECKS = ROOT / "decks/crystallization/single_field"
OUT = ROOT / "results/crystallization/single_field/seeds"
CASES = {  # deck stem: reference free energy
    "a1_planar_interface_pinn": 0.302,
    "a1_planar_interface_hybrid": 0.302,
    "a3_ni_nucleus_x0.200_pinn": 0.9137052,
    "a3_ni_nucleus_x0.200_hybrid": 0.9137052,
}
n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
OUT.mkdir(parents=True, exist_ok=True)
rows, summary = [], []
for stem, ref in CASES.items():
    base = (DECKS / f"{stem}.i").read_text()
    values = []
    for seed in range(n_seeds):
        name = f"{stem}_seed{seed}"
        text = re.sub(r"seed       = 0", f"seed       = {seed}", base, count=1)
        text = text.replace(f"results/crystallization/single_field/{stem}",
                            f"results/crystallization/single_field/seeds/{name}")
        deck = OUT / f"{name}.i"
        deck.write_text(text)
        t0 = time.time()
        env = dict(os.environ, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4")
        proc = subprocess.run(["python3.11", "-m", "pinnacle", "run", str(deck)], cwd=ROOT, env=env,
                              capture_output=True, text=True)
        wall = time.time() - t0
        post = OUT / name / "postprocessors.csv"
        if proc.returncode != 0 or not post.exists():
            rows.append({"case": stem, "seed": seed, "value": float("nan"), "rel_error": float("nan"),
                         "wall_s": round(wall, 1), "status": "failed"})
            continue
        value = float(post.read_text().splitlines()[-1])
        values.append(value)
        log = (OUT / name / "pinnacle.log").read_text()
        m = re.search(r"prediction within ([0-9.e+-]+) .*first step: (\d+)", log)
        rows.append({"case": stem, "seed": seed, "value": value, "rel_error": value / ref - 1.0,
                     "wall_s": round(wall, 1), "status": "ok",
                     "prediction_l2": m.group(1) if m else "", "corrector_iters": m.group(2) if m else ""})
        print(rows[-1], flush=True)
    err = np.array(values) / ref - 1.0
    summary.append({"case": stem, "n_ok": len(values), "n": n_seeds, "mean_rel_error": err.mean(),
                    "std_rel_error": err.std(ddof=1) if len(err) > 1 else 0.0,
                    "worst_abs_rel_error": np.abs(err).max()})
fields = ["case", "seed", "value", "rel_error", "wall_s", "status", "prediction_l2", "corrector_iters"]
with (OUT / "seed_runs.csv").open("w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=fields)
    wr.writeheader()
    wr.writerows(rows)
with (OUT / "seed_table.csv").open("w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
    wr.writeheader()
    wr.writerows(summary)
for s in summary:
    print(s)
