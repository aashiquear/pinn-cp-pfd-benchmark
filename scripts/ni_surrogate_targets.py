#!/usr/bin/env python3
"""Training targets of the Ni specimen surrogate from finished specimen runs.

Deck-driven (PINNACLE .i syntax):

  [Targets]
    design     = ../../../results/ni_study/final/surrogate/data/design.csv   # ni_surrogate_design.py rows
    runs       = ../../../results/ni_study/final/surrogate/h040             # <specimen>/ run dirs (ni_study.py)
    runs_hf    = ../../../results/ni_study/final/surrogate/h020             # optional: the high-fidelity subset
    out        = ../../../results/ni_study/final/surrogate/targets          # <i>/target.npz, samples.csv
    abscissa   = gauge_strain
    columns    = eng_stress
    scalars    = 'failure_gauge_strain location_fraction'
    hf_scalars = location_fraction
    up_to      = failure        # failure: the curve up to the failure step (summary.json); end: the whole run
  []

Writes what PINNACLE's ``parametric_surrogate`` reads with ``samples_file`` and ``data_dir``
(ADR 0030 addendum): ``<out>/<i:04d>/target.npz`` with ``t`` (the abscissa), ``y`` (the
curve columns), ``e`` (the engineering strain, the abscissa itself here, for the
material-point baseline), ``s`` (the scalars, NaN where a run has none) and ``s_hf`` (the
high-fidelity scalars from ``runs_hf``, NaN where that run does not exist); and
``<out>/samples.csv`` with one row per design sample and the columns the surrogate deck's
[Parameters] use (``phi1 phi phi2 t1 t2 t3 t4 t5``; the parser lowercases ``Phi``).

A design row i is specimen S{i+1:02d} of ``ni_surrogate_design.py``; its run directory is
``<runs>/<specimen>`` or ``<runs>/<specimen lower-case>`` (ni_study.py's case naming). The
curve comes from the run's ``results.csv`` (ni_specimen.py: gauge_strain, eng_stress per
converged step), the scalars from ``summary.json`` (failure_gauge_strain, failure_mode, ...)
and ``fracture.json`` (location_fraction, angle_to_axis_deg). Runs that have not finished
are skipped with a note; the surrogate deck can regenerate or ignore them. No solver run.

Usage: ni_surrogate_targets.py TARGETS.i
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_figures import _resolve, _str  # noqa: E402

PARAM_COLUMNS = {"phi1": "phi1_deg", "phi": "Phi_deg", "phi2": "phi2_deg",
                 "t1": "t1_mm", "t2": "t2_mm", "t3": "t3_mm", "t4": "t4_mm", "t5": "t5_mm"}


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "targets" not in cfg:
        raise SystemExit(f"{path}: missing [Targets] block")
    return cfg["targets"]


def _words(v) -> list[str]:
    if isinstance(v, (list, tuple)):
        return [str(x) for x in v]
    return [w for w in str(v or "").replace("[", " ").replace("]", " ").replace("'", " ").replace(",", " ").split()]


def run_dir(root: Path, specimen: str) -> Path | None:
    for name in (specimen, specimen.lower(), specimen.upper()):
        if (root / name / "summary.json").exists():
            return root / name
    return None


def scalars_of(d: Path, names: list[str]) -> np.ndarray:
    s = json.loads((d / "summary.json").read_text())
    f = json.loads((d / "fracture.json").read_text()) if (d / "fracture.json").exists() else {}
    out = np.full(len(names), np.nan)
    for j, n in enumerate(names):
        v = s.get(n, f.get(n))
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            out[j] = float(v)
    return out


def curve_of(d: Path, abscissa: str, columns: list[str], up_to: str) -> tuple[np.ndarray, np.ndarray]:
    rows = list(csv.DictReader((d / "results.csv").open()))
    if up_to == "failure":
        step = json.loads((d / "summary.json").read_text()).get("failure_step")
        if step is not None:
            rows = [r for r in rows if int(float(r["step"])) <= int(step)]
    t = np.array([float(r[abscissa]) for r in rows])
    y = np.array([[float(r[c]) for c in columns] for r in rows])
    return t, y


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    design = list(csv.DictReader(_resolve(_str(cfg["design"]), deck).open()))
    runs = _resolve(_str(cfg["runs"]), deck)
    runs_hf = _resolve(_str(cfg["runs_hf"]), deck) if cfg.get("runs_hf") else None
    out = _resolve(_str(cfg["out"]), deck)
    out.mkdir(parents=True, exist_ok=True)
    abscissa = _str(cfg.get("abscissa", "gauge_strain"))
    columns = _words(cfg.get("columns", "eng_stress"))
    scalars = _words(cfg.get("scalars", ""))
    hf_scalars = _words(cfg.get("hf_scalars", ""))
    up_to = _str(cfg.get("up_to", "failure"))
    missing = [s for s in hf_scalars if s not in scalars]
    if missing:
        raise SystemExit(f"hf_scalars {missing} are not among scalars {scalars}")
    hf_cols = [scalars.index(s) for s in hf_scalars]
    written, skipped, with_hf = 0, [], 0
    with (out / "samples.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "specimen", *PARAM_COLUMNS, "run", "hf_run"])
        w.writeheader()
        for i, row in enumerate(design):
            spec = row.get("specimen", f"S{i + 1:02d}")
            params = {k: float(row[v]) for k, v in PARAM_COLUMNS.items() if v in row}
            d = run_dir(runs, spec)
            d_hf = run_dir(runs_hf, spec) if runs_hf is not None else None
            w.writerow({"sample": i, "specimen": spec, **params, "run": str(d or ""), "hf_run": str(d_hf or "")})
            if d is None:
                skipped.append(spec)
                continue
            t, y = curve_of(d, abscissa, columns, up_to)
            arrays = {"t": t, "y": y, "e": t if abscissa.endswith("strain") else t, "run_time_s": np.asarray(0.0)}
            if scalars:
                arrays["s"] = scalars_of(d, scalars)
            if hf_scalars:
                arrays["s_hf"] = scalars_of(d_hf, scalars)[hf_cols] if d_hf is not None else np.full(len(hf_cols), np.nan)
                with_hf += int(d_hf is not None)
            sdir = out / f"{i:04d}"
            sdir.mkdir(parents=True, exist_ok=True)
            np.savez(sdir / "target.npz", **arrays)
            written += 1
    print(f"wrote {written} of {len(design)} targets under {out} ({with_hf} with a high-fidelity run); "
          f"samples.csv has the [Parameters] columns {list(PARAM_COLUMNS)}")
    if skipped:
        print("no finished run for:", " ".join(skipped))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
