#!/usr/bin/env python3
"""CP parameter study driver: one study deck -> per-case pipeline decks, runs, curves, errors.

The study deck (PINNACLE .i syntax) names a BASE pipeline deck and the
parameter levels to vary, one parameter at a time (OAT), so each curve isolates
what one parameter does to the stress-strain response:

  [Study]
    base      = ../cp_study/h40_base.i     # pipeline deck with the base [Mechanics]
    out_root  = ../../results/cp_study     # cases land in <out_root>/<name>
    markers   = ../../data/digitized/al2219_fig3_experiment_markers.csv
    concurrent = 4                          # runs at a time (each is single-threaded)
  [Sweep]
    rho0 = '55 110 440'                     # levels; the base value is the reference case
    k1   = '0.1716 0.6864'
  [Cases]                                   # optional named multi-parameter cases
    [./fit1]  rho0 = 110  k1 = 0.5  []

Stages:
  prepare  write <out_root>/<case>/pipeline.i (the base with one value replaced)
           and run label/orient/deck for it;
  run      launch `pinnacle run` for every prepared case that has no
           solver_history.csv yet, `concurrent` at a time, and wait;
  collect  build results.csv per case (pipeline `curve`), then
           <out_root>/<study>_summary.csv with the error of every case against the
           experiment markers: RMSE and max |dev| of the simulated stress interpolated
           at the marker strains, and the relative error at each marker.

Usage: cp_parameter_study.py STUDY.i prepare|run|collect|all
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(STUDY_ROOT / "scripts"))
from polycrystal_plate import read_pipeline_deck, _resolve  # noqa: E402

PIPE = STUDY_ROOT / "scripts/polycrystal_plate.py"
PINNACLE = os.environ.get("PINNACLE_BIN", "pinnacle")


def _read_study(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "study" not in cfg:
        raise SystemExit(f"{path}: missing [Study]")
    return cfg


def _levels(v) -> list[float]:
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(x) for x in str(v).split()]


def cases_of(cfg: dict) -> list[tuple[str, dict]]:
    out: list[tuple[str, dict]] = [("base", {})]
    for key, levels in (cfg.get("sweep") or {}).items():
        for v in _levels(levels):
            out.append((f"{key}_{v:g}", {key: v}))
    for name, params in (cfg.get("cases") or {}).items():
        out.append((name, {k: float(v) for k, v in params.items()}))
    return out


def stage_prepare(cfg: dict, study: Path) -> list[Path]:
    base = _resolve(str(cfg["study"]["base"]), study)
    out_root = _resolve(str(cfg["study"]["out_root"]), study)
    base_text = base.read_text()
    written = []
    for name, params in cases_of(cfg):
        case_dir = out_root / name
        case_dir.mkdir(parents=True, exist_ok=True)
        text = base_text
        for k, v in params.items():
            pat = re.compile(rf"^(\s*{re.escape(k)}\s*=\s*)\S+", re.M)
            if not pat.search(text):
                raise SystemExit(f"{base}: [Mechanics] has no key {k!r} to replace")
            text = pat.sub(rf"\g<1>{v:.12g}", text, count=1)
        text = re.sub(r"^(\s*out_dir\s*=\s*)\S+", rf"\g<1>{case_dir}", text, count=1, flags=re.M)
        # Paths in the base deck are relative to ITS directory; the case deck
        # lives elsewhere, so resolve them.
        for key in ("run_dir", "deck"):
            m = re.search(rf"^(\s*{key}\s*=\s*)(\S+)", text, re.M)
            if m:
                text = text[: m.start(2)] + str(_resolve(m.group(2), base)) + text[m.end(2):]
        hdr = f"# Case '{name}' of study {study.name}: {params or 'base values'}\n"
        pipe = case_dir / "pipeline.i"
        pipe.write_text(hdr + text)
        (case_dir / "case.json").write_text(json.dumps({"name": name, "params": params}, indent=1))
        subprocess.run([sys.executable, str(PIPE), str(pipe), "label", "orient", "deck"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        written.append(pipe)
        print("prepared", name, params)
    return written


def stage_run(cfg: dict, study: Path) -> None:
    out_root = _resolve(str(cfg["study"]["out_root"]), study)
    n_conc = int(cfg["study"].get("concurrent", 4))
    todo = []
    for name, _ in cases_of(cfg):
        d = out_root / name
        if not (d / "cpfe.i").exists():
            raise SystemExit(f"{d}: not prepared")
        if not (d / "solver" / "solver_history.csv").exists():
            todo.append(d)
    print(f"{len(todo)} cases to run, {n_conc} at a time")
    running: list[tuple[Path, subprocess.Popen]] = []
    while todo or running:
        while todo and len(running) < n_conc:
            d = todo.pop(0)
            log = (d / "run.log").open("w")
            p = subprocess.Popen([PINNACLE, "run", "cpfe.i"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
            running.append((d, p))
            print(time.strftime("%H:%M:%S"), "started", d.name)
        time.sleep(20)
        still = []
        for d, p in running:
            if p.poll() is None:
                still.append((d, p))
            else:
                print(time.strftime("%H:%M:%S"), "finished", d.name, "exit", p.returncode)
        running = still


def _markers(path: Path):
    lines = [ln for ln in path.open() if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    return np.array([float(r["strain_pct"]) / 100 for r in rows]), np.array([float(r["stress_MPa"]) for r in rows])


def stage_collect(cfg: dict, study: Path) -> None:
    out_root = _resolve(str(cfg["study"]["out_root"]), study)
    ex_e, ex_s = _markers(_resolve(str(cfg["study"]["markers"]), study))
    # The origin marker (0.0026% strain) sits below the first load step and
    # would dominate a relative error; the FE curve is compared from 0.1% on.
    keep = ex_e >= 1.0e-3
    ex_e, ex_s = ex_e[keep], ex_s[keep]
    rows = []
    for name, params in cases_of(cfg):
        d = out_root / name
        if not (d / "solver" / "solver_history.csv").exists():
            print("skip (no run)", name)
            continue
        subprocess.run([sys.executable, str(PIPE), str(d / "pipeline.i"), "curve"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        r = list(csv.DictReader((d / "results.csv").open()))
        e = np.array([float(x["nominal_strain"]) for x in r]); s = np.array([float(x["eng_stress"]) for x in r])
        reached = float(e.max())
        inside = ex_e <= reached + 1e-9
        sim = np.interp(ex_e[inside], e, s)
        dev = sim - ex_s[inside]
        rel = dev / ex_s[inside]
        rec = {"case": name, **{f"p_{k}": v for k, v in params.items()},
               "strain_reached": reached, "n_markers_covered": int(inside.sum()),
               "rmse_MPa": float(np.sqrt(np.mean(dev ** 2))) if inside.any() else np.nan,
               "max_abs_dev_MPa": float(np.abs(dev).max()) if inside.any() else np.nan,
               "mean_rel_err": float(np.mean(rel)) if inside.any() else np.nan,
               "stress_at_2pct": float(np.interp(0.02, e, s)) if reached >= 0.02 else np.nan,
               "stress_at_10pct": float(np.interp(0.10, e, s)) if reached >= 0.10 else np.nan,
               "stress_at_17pct": float(np.interp(0.17, e, s)) if reached >= 0.17 else np.nan,
               "cutbacks": int(sum(int(x["cutbacks"]) for x in r))}
        for k, (ee, ss, sv) in enumerate(zip(ex_e[inside], ex_s[inside], sim)):
            rec[f"rel_err_marker{k+1}_{ee*100:.2g}pct"] = float((sv - ss) / ss)
        rows.append(rec)
    keys = []
    for rec in rows:
        for k in rec:
            if k not in keys:
                keys.append(k)
    summ = out_root / f"{study.stem}_summary.csv"
    with summ.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print("wrote", summ)
    for rec in rows:
        print(rec["case"], "reached", round(rec["strain_reached"], 4), "rmse", round(rec["rmse_MPa"], 1) if rec["rmse_MPa"] == rec["rmse_MPa"] else "nan",
              "mean rel err", round(rec["mean_rel_err"], 3) if rec["mean_rel_err"] == rec["mean_rel_err"] else "nan")


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    study = Path(argv[0]).resolve(); cfg = _read_study(study)
    stages = ["prepare", "run", "collect"] if argv[1] == "all" else [argv[1]]
    for st in stages:
        {"prepare": stage_prepare, "run": stage_run, "collect": stage_collect}[st](cfg, study)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
