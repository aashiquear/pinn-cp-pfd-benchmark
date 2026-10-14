#!/usr/bin/env python3
"""Calibrate Gc for each length scale so fracture initiates at the target strain.

Deck-driven (PINNACLE .i syntax):

  [Calibrate]
    base          = ../j2/base_h40.i           # pipeline deck (scripts/pfd_plate.py)
    out_root      = ../../../results/pfd_study/j2/cal_h40
    l0            = '2.5 5 10'                 # length scales (each gets its own bisection)
    Gc_lo         = 10                         # bracket of Gc (MPa um = J/m^2), log bisection
    Gc_hi         = 4000
    target_strain = 0.172                      # fracture strain: last Vaughan point
    criterion     = node                       # node | element (d >= damage_threshold)
    tolerance     = 0.0008                     # accept when |strain - target| <= one load step
    max_iter      = 9
    concurrent    = 7                          # bisections run in parallel over l0, one run each
    threads       = 2

Bisection in log Gc: a run whose fracture strain is below the target (or
that fails) has too small a Gc (lo moves up); a run that never initiates
up to the end of the path has too large a Gc (hi moves down). Every run is
kept under <out_root>/l<l0>/gc<Gc>/ with its pipeline deck, so the whole
history is reproducible; <out_root>/calibration.csv records every run and
<out_root>/calibrated.json the Gc(l0) pairs with their fracture strains.

Usage: pfd_calibrate.py CALIBRATE.i
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STUDY_ROOT = Path(__file__).resolve().parent.parent
PIPE = STUDY_ROOT / "scripts/pfd_plate.py"
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(STUDY_ROOT / "scripts"))

from pfd_study import _resolve, _set_keys  # noqa: E402


def read_cal(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "calibrate" not in cfg:
        raise SystemExit(f"{path}: missing [Calibrate]")
    return cfg["calibrate"]


def _levels(v) -> list[float]:
    return [float(x) for x in (v if isinstance(v, (list, tuple)) else str(v).split())]


def run_case(base: Path, cdir: Path, values: dict[str, str], threads: int) -> dict:
    cdir.mkdir(parents=True, exist_ok=True)
    text = base.read_text()
    for m in re.finditer(r"^(\s*(\w+)\s*=\s*)(\.\.?/[^\s#]+)", base.read_text(), re.M):
        text = text.replace(m.group(0), m.group(1) + os.path.relpath(_resolve(m.group(3), base), cdir), 1)
    text = _set_keys(text, {"out_dir": ".", "threads": str(threads), **values})
    pipe = cdir / "pipeline.i"
    pipe.write_text(text)
    if not (cdir / "fracture.json").exists():
        subprocess.run([sys.executable, str(PIPE), str(pipe), "all"], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if not (cdir / "fracture.json").exists():
        return {"failed": True}
    return json.loads((cdir / "fracture.json").read_text())


def bisect(l0: float, cal: dict, base: Path, out_root: Path, log: list[dict]) -> dict:
    lo, hi = math.log(float(cal["Gc_lo"])), math.log(float(cal["Gc_hi"]))
    target = float(cal.get("target_strain", 0.172))
    tol = float(cal.get("tolerance", 0.0008))
    key = "fracture_strain_node" if str(cal.get("criterion", "node")) == "node" else "fracture_strain_element"
    best = None
    for it in range(int(cal.get("max_iter", 9))):
        gc = math.exp(0.5 * (lo + hi))
        cdir = out_root / f"l{l0:g}" / f"gc{gc:.6g}"
        t0 = time.time()
        frac = run_case(base, cdir, {"l0": f"{l0:g}", "Gc": f"{gc:.6g}"}, int(cal.get("threads", 2)))
        fs = frac.get(key)
        rec = {"l0": l0, "iteration": it, "Gc": gc, "fracture_strain": fs, "failed": bool(frac.get("failed")),
               "steps_converged": frac.get("steps_converged"), "wall_s": time.time() - t0, "dir": str(cdir)}
        log.append(rec)
        print(time.strftime("%H:%M:%S"), {k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items() if k != "dir"}, flush=True)
        if frac.get("failed"):
            lo = math.log(gc)  # a failed run is treated as too brittle: raise Gc
            continue
        if fs is None:
            hi = math.log(gc)  # never initiated: Gc too high
            continue
        if best is None or abs(fs - target) < abs(best["fracture_strain"] - target):
            best = rec
        if abs(fs - target) <= tol:
            break
        if fs < target:
            lo = math.log(gc)
        else:
            hi = math.log(gc)
    return {"l0": l0, **({k: best[k] for k in ("Gc", "fracture_strain", "dir")} if best else {"Gc": None})}


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cal = read_cal(deck)
    base = _resolve(str(cal["base"]), deck)
    out_root = _resolve(str(cal["out_root"]), deck)
    out_root.mkdir(parents=True, exist_ok=True)
    log: list[dict] = []
    with ThreadPoolExecutor(max_workers=int(cal.get("concurrent", 4))) as ex:
        results = list(ex.map(lambda l0: bisect(l0, cal, base, out_root, log), _levels(cal["l0"])))
    with (out_root / "calibration.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log[0]))
        w.writeheader()
        w.writerows(sorted(log, key=lambda r: (r["l0"], r["iteration"])))
    (out_root / "calibrated.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
