#!/usr/bin/env python3
"""Calibrate G_c of the Ni CP+PFD specimen, at each length scale, to a failure strain.

Deck-driven (PINNACLE .i syntax):

  [Calibrate]
    base        = base_L1_h020.i                # specimen deck (scripts/ni_specimen.py) with Gc/l0 keys
    out_root    = ../../../results/ni_study/pfd/cal_L1_h020
    l0          = '0.4 0.6 0.8'                 # length scales (mm), one search each, run side by side
    Gc_guess    = '11 16 22'                    # first estimate per l0 (MPa mm)
    first       = '0.5 2'                       # the first round runs Gc_guess x each factor, concurrently
    target      = 0.463                         # failure (gauge) strain: Fig. 4(b), L.1
    tolerance   = 0.004                         # accept |failure strain - target| <= this (two load steps)
    max_runs    = 6                             # runs per length scale
    threads     = 2

The failure strain of a run is the gauge strain of its first separated step
(PINNACLE stop_separation: d > 0.85 cuts the grips apart; `ni_specimen.py
curve`, summary.json failure_gauge_strain). A run that fails below the target
has too small a G_c; a run that fails above it, or not at all by the end of
the path, too large a G_c. A run that dies before the target without failing
decides nothing and stops that length scale's search (noted, not guessed).

Search in ln G_c per length scale: the first round brackets the guess; then
regula falsi between the closest runs on either side of the target when both
failed (failure strain linear in ln G_c), the geometric midpoint when the
upper one never failed, and a factor 3 outward while one side is missing.

Every run is kept under <out_root>/l<l0>/gc<Gc>/ with its specimen deck;
<out_root>/calibration.csv records every run, <out_root>/calibrated.json the
G_c(l0) pairs with their failure strains and fracture descriptors.

Usage: ni_pfd_calibrate.py CALIBRATE.i
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STUDY_ROOT = Path(__file__).resolve().parent.parent
PIPE = STUDY_ROOT / "scripts/ni_specimen.py"
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(STUDY_ROOT / "scripts"))

from ni_study import _resolve, _set_keys  # noqa: E402

LOCK = threading.Lock()


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
    """Render the specimen deck into cdir, run deck/run/curve/fracture, return summary + fracture."""
    cdir.mkdir(parents=True, exist_ok=True)
    text = base.read_text()
    for m in re.finditer(r"^(\s*(\w+)\s*=\s*)(\.\.?/[^\s#]+)", base.read_text(), re.M):
        text = text.replace(m.group(0), m.group(1) + os.path.relpath(_resolve(m.group(3), base), cdir), 1)
    text = _set_keys(text, {"out_dir": ".", "threads": str(threads), **values})
    deck = cdir / "specimen.i"
    deck.write_text(text)
    if not (cdir / "fracture.json").exists():
        with (cdir / "pipeline.log").open("w") as fh:
            subprocess.run([sys.executable, str(PIPE), str(deck), "deck", "run", "curve", "fracture"],
                           check=False, stdout=fh, stderr=subprocess.STDOUT)
    if not (cdir / "summary.json").exists():
        return {"failed": True}
    out = json.loads((cdir / "summary.json").read_text())
    if (cdir / "fracture.json").exists():
        out["fracture"] = json.loads((cdir / "fracture.json").read_text())
    return out


def search(l0: float, guess: float, cal: dict, base: Path, out_root: Path, log: list[dict]) -> dict:
    target = float(cal.get("target", 0.463))
    tol = float(cal.get("tolerance", 0.004))
    threads = int(cal.get("threads", 2))
    pts: list[tuple[float, float | None]] = []      # (ln Gc, failure strain or None = no failure)

    def one(gc: float, it: int) -> dict:
        cdir = out_root / f"l{l0:g}" / f"gc{gc:.4g}"
        t0 = time.time()
        res = run_case(base, cdir, {"l0": f"{l0:g}", "Gc": f"{gc:.4g}"}, threads)
        fs = res.get("failure_gauge_strain")
        frac = res.get("fracture", {})
        rec = {"l0": l0, "run": it, "Gc": float(f"{gc:.4g}"), "failure_strain": fs,
               "final_gauge_strain": res.get("final_gauge_strain"), "run_died": bool(res.get("run_died")),
               "failed_pipeline": bool(res.get("failed")), "location_fraction": frac.get("location_fraction"),
               "angle_to_axis_deg": frac.get("angle_to_axis_deg"), "peak_eng_stress": res.get("peak_eng_stress"),
               "wall_h": (time.time() - t0) / 3600.0, "dir": str(cdir), "note": ""}
        with LOCK:
            log.append(rec)
            print(time.strftime("%H:%M:%S"), {k: (round(v, 4) if isinstance(v, float) else v)
                                              for k, v in rec.items() if k != "dir"}, flush=True)
        return rec

    def judge(rec: dict) -> str | None:
        """None when the run decides nothing; else 'low' (fails early), 'high' or 'hit'."""
        fs = rec["failure_strain"]
        if rec["failed_pipeline"]:
            rec["note"] = "pipeline failed"; return None
        if fs is None:
            stop = rec["final_gauge_strain"]
            if rec["run_died"] and (stop is None or float(stop) < target - tol):
                rec["note"] = "died at %.4f before the target without failing" % float(stop or 0.0)
                return None
            return "high"
        if abs(fs - target) <= tol:
            return "hit"
        return "low" if fs < target else "high"

    first = _levels(cal.get("first", "0.5 2"))
    with ThreadPoolExecutor(max_workers=len(first)) as ex:
        recs = list(ex.map(lambda a: one(guess * a[1], a[0]), enumerate(first)))
    n_run = len(recs)
    best = None
    while True:
        for rec in recs:
            verdict = judge(rec)
            if verdict is None:
                return {"l0": l0, "Gc": None, "note": rec["note"], "dir": rec["dir"]}
            pts.append((math.log(rec["Gc"]), rec["failure_strain"]))
            if rec["failure_strain"] is not None and (
                    best is None or abs(rec["failure_strain"] - target) < abs(best["failure_strain"] - target)):
                best = rec
            if verdict == "hit":
                return {"l0": l0, **{k: rec[k] for k in ("Gc", "failure_strain", "location_fraction", "angle_to_axis_deg", "dir")}}
        if n_run >= int(cal.get("max_runs", 6)):
            break
        lows = [p for p in pts if p[1] is not None and p[1] < target]
        highs = [p for p in pts if p[1] is None or p[1] > target]
        lo = max(lows) if lows else None
        hi = min(highs) if highs else None
        if lo and hi:
            if hi[1] is not None:
                x = lo[0] + (target - lo[1]) * (hi[0] - lo[0]) / (hi[1] - lo[1])
                w = hi[0] - lo[0]
                x = min(max(x, lo[0] + 0.1 * w), hi[0] - 0.1 * w)
            else:
                x = 0.5 * (lo[0] + hi[0])
        elif lo:
            x = lo[0] + math.log(3.0)
        else:
            x = hi[0] - math.log(3.0)
        recs = [one(math.exp(x), n_run)]
        n_run += 1
    out = {"l0": l0, "Gc": best["Gc"] if best else None, "note": "run budget spent; closest run kept"}
    if best:
        out.update({k: best[k] for k in ("failure_strain", "location_fraction", "angle_to_axis_deg", "dir")})
    return out


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cal = read_cal(deck)
    base = _resolve(str(cal["base"]), deck)
    out_root = _resolve(str(cal["out_root"]), deck)
    out_root.mkdir(parents=True, exist_ok=True)
    l0s, guesses = _levels(cal["l0"]), _levels(cal["Gc_guess"])
    if len(l0s) != len(guesses):
        raise SystemExit("Gc_guess needs one value per l0")
    log: list[dict] = []
    with ThreadPoolExecutor(max_workers=len(l0s)) as ex:
        results = list(ex.map(lambda a: search(a[0], a[1], cal, base, out_root, log), zip(l0s, guesses, strict=True)))
    with (out_root / "calibration.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log[0]))
        w.writeheader()
        w.writerows(sorted(log, key=lambda r: (r["l0"], r["run"])))
    (out_root / "calibrated.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
