#!/usr/bin/env python3
"""Design, run and collect the FE samples of a plate surrogate (ADR 0030 data source).

The PINNACLE ``parametric_surrogate`` executioner can generate its samples in
process, one after the other. For the plate (10 min per FE run) the samples
are better produced by the study pipeline, several at a time; this script
is the bridge. Deck-driven (PINNACLE ``.i`` syntax):

  [Design]
    base       = ../cal_center2.i                 # pipeline deck of the calibrated plate
    driver     = cp                               # cp: cp_parameter_study.py ([Mechanics] keys);
                                                  # pfd: pfd_study.py ([Plate] keys, Part II)
    out_root   = ../../../results/cp_study/surrogate/h40_data
    samples    = 24
    sampler    = lhs
    seed       = 0
    concurrent = 4
    columns    = 'eng_stress'                     # results.csv columns to store as y(t)
    abscissa   = nominal_strain
  [Parameters]                                    # PINNACLE syntax: path into the PIPELINE deck's [Mechanics]
    [./rho0]   path = mechanics/rho0   range = '23.99 95.97'  scale = log []
    ...

Stages:
  design   sample the rows (PINNACLE ParameterSpace, same sampler and seed
           semantics as the executioner), write <out_root>/design.csv and a
           study deck <out_root>/study.i with one [Cases] entry per row for
           scripts/cp_parameter_study.py;
  run      = cp_parameter_study.py study.i all (prepare, run concurrently, collect);
  collect  write <out_root>/samples/<i>/target.npz (t = abscissa, y = columns,
           run_time_s) and <out_root>/samples.csv in the layout the executioner
           reads with `samples_file = design.csv` and `data_dir = samples`;
           with `fields = displacement` also <i>/field.npz (the nodal field of
           the last VTK frame) for a field-target surrogate.

Usage: surrogate_design.py DESIGN.i design|run|collect|all
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(STUDY_ROOT / "scripts"))


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("design", "parameters"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _words(v) -> list[str]:
    return [str(x) for x in v] if isinstance(v, (list, tuple)) else str(v).split()


def stage_design(cfg: dict, deck: Path) -> Path:
    from pinnacle.pinn.parametric import ParameterSpace
    d = cfg["design"]
    out_root = _resolve(str(d["out_root"]), deck); out_root.mkdir(parents=True, exist_ok=True)
    space = ParameterSpace(dict(cfg["parameters"]))
    rows = space.sample(int(d.get("samples", 24)), str(d.get("sampler", "lhs")), int(d.get("seed", 0)))
    with (out_root / "design.csv").open("w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["sample", *space.names])
        for i, r in enumerate(rows):
            w.writerow([i, *[f"{v:.12g}" for v in r]])
    import os
    base = _resolve(str(d["base"]), deck)
    rel = lambda q: os.path.relpath(q, out_root)  # noqa: E731  # paths relative to study.i: host-independent
    driver = str(d.get("driver", "cp")).lower()
    lines = [f"# Surrogate design of {deck.name}: {len(rows)} {d.get('sampler', 'lhs')} samples, seed {d.get('seed', 0)}.",
             "[Study]", f"  base       = {rel(base)}", "  out_root   = cases",
             *([f"  markers    = {rel(STUDY_ROOT / 'data/digitized/al2219_fig3_experiment_markers.csv')}"] if driver == "cp" else []),
             f"  concurrent = {int(d.get('concurrent', 4))}", "[]", "[Cases]"]
    for i, r in enumerate(rows):
        lines.append(f"  [./s{i:04d}]")
        for name, path, v in zip(space.names, space.paths, r, strict=True):
            lines.append(f"    {path[-1]} = {v:.12g}")
        lines.append("  []")
    lines.append("[]")
    (out_root / "study.i").write_text("\n".join(lines) + "\n")
    print("design:", len(rows), "rows ->", out_root / "design.csv", out_root / "study.i")
    return out_root


def stage_run(cfg: dict, deck: Path) -> None:
    out_root = _resolve(str(cfg["design"]["out_root"]), deck)
    driver = "pfd_study.py" if str(cfg["design"].get("driver", "cp")).lower() == "pfd" else "cp_parameter_study.py"
    subprocess.run([sys.executable, str(STUDY_ROOT / "scripts" / driver), str(out_root / "study.i"), "all"], check=True)


def _converged(v) -> bool:
    """The history writes the flag as a float (``1.000e+00``) or a word."""
    t = str(v).strip().lower()
    return t == "true" or (t not in ("", "false") and float(t) == 1.0)


def _cell_to_nodes(values: np.ndarray, cells: np.ndarray, n_nodes: int) -> np.ndarray:
    """Average a cell field onto the nodes (each node takes the mean of its cells)."""
    acc = np.zeros((n_nodes, *values.shape[1:]), dtype=float); cnt = np.zeros(n_nodes, dtype=float)
    for c, val in zip(cells, values, strict=True):
        for n in c:
            acc[n] += val; cnt[n] += 1.0
    return acc / np.maximum(cnt, 1.0)[(slice(None),) + (None,) * (values.ndim - 1)]


def stage_collect(cfg: dict, deck: Path) -> None:
    d = cfg["design"]
    out_root = _resolve(str(d["out_root"]), deck)
    cols = _words(d.get("columns", "eng_stress")); absc = str(d.get("abscissa", "nominal_strain"))
    fields = _words(d.get("fields", ""))
    design = list(csv.DictReader((out_root / "design.csv").open()))
    sdir = out_root / "samples"; sdir.mkdir(exist_ok=True)
    rows_out = []
    for r in design:
        i = int(r["sample"]); case = out_root / "cases" / f"s{i:04d}"
        res = case / "results.csv"
        if not res.exists():
            print("missing", case); continue
        data = list(csv.DictReader(res.open()))
        t = np.array([float(x[absc]) for x in data]); y = np.array([[float(x[c]) for c in cols] for x in data])
        hist = [h for h in csv.DictReader((case / "solver" / "solver_history.csv").open()) if _converged(h.get("converged", "1"))]
        run_time = float(np.sum([float(h.get("newton_time_s", 0) or 0) for h in hist]))
        if run_time == 0.0 or (hist and "mech_iterations" in hist[0]):
            # coupled (staggered) runs: the Newton time is a fraction of the step (docs/10 §8);
            # the sample's cost is the wall time of the run from the PINNACLE log stamps
            from pfd_study import _wall_seconds
            run_time = float(_wall_seconds(case / "solver" / "pinnacle.log"))
        out = sdir / f"{i:04d}"; out.mkdir(exist_ok=True)
        np.savez(out / "target.npz", t=t, y=y, run_time_s=np.asarray(run_time))
        if fields:
            # One data set per field, in the executioner's layout (samples_<field>/<i>/target.npz,
            # y = the nodal field of the last frame), for field-target surrogate decks.
            from vtk_legacy import read_vtk, series_files
            files = series_files(case / "solver")
            if files:
                v = read_vtk(files[-1][1])
                for f in fields:
                    arr = v["point"].get(f)
                    if arr is None and f in v["cell"]:
                        arr = _cell_to_nodes(np.asarray(v["cell"][f]), np.asarray(v["cells"]), len(v["nodes"]))
                    if arr is None:
                        continue
                    fdir = out_root / f"samples_{f}" / f"{i:04d}"; fdir.mkdir(parents=True, exist_ok=True)
                    np.savez(fdir / "target.npz", y=np.asarray(arr, dtype=float).reshape(len(v["nodes"]), -1), run_time_s=np.asarray(run_time))
        rows_out.append({"sample": i, **{k: r[k] for k in r if k != "sample"}, "run_time_s": run_time, "steps": len(hist)})
    with (out_root / "samples.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0])); w.writeheader(); w.writerows(rows_out)
    print("collected", len(rows_out), "samples; FE time total %.0f s" % sum(r["run_time_s"] for r in rows_out))


def main(argv):
    if len(argv) < 3:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); cfg = read_deck(deck)
    stages = argv[2:]
    if stages == ["all"]:
        stages = ["design", "run", "collect"]
    for st in stages:
        {"design": stage_design, "run": stage_run, "collect": stage_collect}[st](cfg, deck)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
