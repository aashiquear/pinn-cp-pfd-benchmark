#!/usr/bin/env python3
"""PFD study driver: one study deck -> per-case pipeline decks, concurrent runs, summary.

The study deck (PINNACLE .i syntax) names a BASE pipeline deck of
scripts/pfd_plate.py and the cases to run, one parameter at a time or named:

  [Study]
    base       = ../j2/base_h40.i            # pipeline deck with the base [Plate]
    out_root   = ../../../results/pfd_study/j2/sens_h40
    concurrent = 7                           # runs at a time (each uses [Plate] threads)
  [Sweep]
    Gc = '100 400'                           # levels of one key at a time (OAT); the
    l0 = '2.5 10'                            # base is the reference case
  [Cases]                                    # optional named multi-key cases
    [./fine]  l0 = 2.5  Gc = 50  []          # one key per line in a real deck

Any key of the base [Plate] block can be swept (mesh, Gc, l0, eta, W_c,
rho0, ...). Stages:
  prepare  write <out_root>/<case>/pipeline.i (the base with the values replaced,
           out_dir pointed at the case directory) and render its PINNACLE deck;
  run      launch `pfd_plate.py pipeline.i run` for every prepared case without a
           solver_history.csv, `concurrent` at a time, and wait;
  collect  curve every case and write <out_root>/<study>_summary.csv: the case's
           key values, fracture strain (node and element criteria), peak stress,
           fracture energy, final stress and damage, steps, Newton iterations,
           alternations and the wall time of the run (from the PINNACLE log).

Usage: pfd_study.py STUDY.i prepare|run|collect|all
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

STUDY_ROOT = Path(__file__).resolve().parent.parent
PIPE = STUDY_ROOT / "scripts/pfd_plate.py"
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")


def read_study(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "study" not in cfg:
        raise SystemExit(f"{path}: missing [Study]")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _levels(v) -> list[str]:
    return [str(x) for x in v] if isinstance(v, (list, tuple)) else str(v).split()


def cases_of(cfg: dict) -> list[tuple[str, dict[str, str]]]:
    out: list[tuple[str, dict[str, str]]] = [("base", {})]
    for key, levels in (cfg.get("sweep") or {}).items():
        for lv in _levels(levels):
            out.append((f"{key}_{lv}".replace(".", "p").replace("-", "m"), {key: lv}))
    for name, kv in (cfg.get("cases") or {}).items():
        out.append((str(name), {k: str(v) for k, v in kv.items()}))
    return out


def _set_keys(text: str, values: dict[str, str]) -> str:
    """Replace ``key = value`` lines inside the [Plate] block (added when absent)."""
    m = re.search(r"^\[Plate\]", text, re.M)
    if not m:
        raise SystemExit("base pipeline deck has no [Plate] block")
    i = m.start()
    j = text.index("\n[]", i)
    seg = text[i:j]
    for key, val in values.items():
        pat = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)[^#\n]*", re.M)
        if pat.search(seg):
            seg = pat.sub(rf"\g<1>{val}", seg, count=1)
        else:
            seg += f"\n  {key} = {val}"
    return text[:i] + seg + text[j:]


def stage_prepare(cfg: dict, study: Path) -> list[Path]:
    st = cfg["study"]
    base = _resolve(str(st["base"]), study)
    out_root = _resolve(str(st["out_root"]), study)
    out_root.mkdir(parents=True, exist_ok=True)
    text = base.read_text()
    pipes = []
    for name, kv in cases_of(cfg):
        cdir = out_root / name
        cdir.mkdir(parents=True, exist_ok=True)
        # paths in the base deck are relative to the base deck's directory;
        # the case deck lives elsewhere, so every relative path is re-based.
        case_text = text
        for m in re.finditer(r"^(\s*(\w+)\s*=\s*)(\.\.?/[^\s#]+)", text, re.M):
            case_text = case_text.replace(m.group(0), m.group(1) + os.path.relpath(_resolve(m.group(3), base), cdir), 1)
        case_text = _set_keys(case_text, {"out_dir": ".", **kv})
        pipe = cdir / "pipeline.i"
        pipe.write_text(case_text)
        (cdir / "case.json").write_text(json.dumps({"name": name, "params": kv, "base": str(base)}, indent=1))
        subprocess.run([sys.executable, str(PIPE), str(pipe), "deck"], check=True, stdout=subprocess.DEVNULL)
        pipes.append(pipe)
    print("prepared", len(pipes), "cases under", out_root)
    return pipes


def stage_run(cfg: dict, study: Path) -> None:
    st = cfg["study"]
    out_root = _resolve(str(st["out_root"]), study)
    n_conc = int(st.get("concurrent", 4))
    todo = [out_root / name for name, _ in cases_of(cfg) if not (out_root / name / "solver" / "solver_history.csv").exists()]
    running: list[tuple[Path, subprocess.Popen]] = []
    while todo or running:
        while todo and len(running) < n_conc:
            cdir = todo.pop(0)
            proc = subprocess.Popen([sys.executable, str(PIPE), str(cdir / "pipeline.i"), "run"], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            running.append((cdir, proc))
            print(time.strftime("%H:%M:%S"), "started", cdir.name)
        time.sleep(20)
        still = []
        for cdir, proc in running:
            if proc.poll() is None:
                still.append((cdir, proc))
            else:
                print(time.strftime("%H:%M:%S"), "finished", cdir.name, "exit", proc.returncode)
        running = still


def _wall_seconds(log: Path) -> float:
    if not log.exists():
        return float("nan")
    stamps = re.findall(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)", log.read_text(), re.M)
    if len(stamps) < 2:
        return float("nan")
    f = "%Y-%m-%d %H:%M:%S"
    return (datetime.strptime(stamps[-1], f) - datetime.strptime(stamps[0], f)).total_seconds()


def stage_collect(cfg: dict, study: Path) -> None:
    st = cfg["study"]
    out_root = _resolve(str(st["out_root"]), study)
    rows = []
    for name, kv in cases_of(cfg):
        cdir = out_root / name
        if not (cdir / "solver" / "solver_history.csv").exists():
            print("missing", cdir)
            continue
        subprocess.run([sys.executable, str(PIPE), str(cdir / "pipeline.i"), "curve"], check=True, stdout=subprocess.DEVNULL)
        frac = json.loads((cdir / "fracture.json").read_text())
        rows.append({"case": name, **{f"param_{k}": v for k, v in kv.items()}, **{k: v for k, v in frac.items() if not isinstance(v, (list, dict))},
                     "wall_s": _wall_seconds(cdir / "solver" / "pinnacle.log")})
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "case", not k.startswith("param_"), k))
    path = out_root / f"{study.stem}_summary.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k in ("case", "fracture_strain_node", "fracture_strain_element", "peak_eng_stress", "fracture_energy_MPa", "steps_converged", "newton_iterations", "wall_s")})
    print("wrote", path)


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    study = Path(argv[1]).resolve()
    cfg = read_study(study)
    stages = ["prepare", "run", "collect"] if argv[2] == "all" else [argv[2]]
    for s in stages:
        {"prepare": stage_prepare, "run": stage_run, "collect": stage_collect}[s](cfg, study)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
