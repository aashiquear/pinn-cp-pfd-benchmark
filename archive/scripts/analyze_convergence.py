#!/usr/bin/env python3
"""Pairwise mesh convergence of a calibrated parameter set.

Convergence here is a property of the PARAMETER SET, not of a single mesh:
a set is converged at density d when densities d-1 and d produce
substantially the same result. Two measures, one per parameter group:

* CP set -- relative L2 difference between the stress-strain curves of
  consecutive densities, compared on a shared strain grid (the densities do
  not share load-step strains exactly once cutbacks occur).
* Gc -- absolute difference in damage-initiation strain between consecutive
  densities.

The reported answer is the LEAST DENSE mesh at which the pair meets the
tolerance, since that is the mesh the solver comparison should then use.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent


def _curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Stress-strain from a run, VALID STEPS ONLY.

    A run whose directory carries INVALID.txt has steps the harness refused
    to certify -- the boundary strain stopped tracking the prescribed value,
    so later rows are the solver's abandoned state, not a solution. Reading
    the CSV past that point produced negative stresses at small strain and
    6,784 MPa at 13%, and quoting those as a result is exactly the failure
    the marker exists to prevent. The run.json records the prescribed strain
    program, so the same invariant can be re-applied here.
    """
    rows = list(csv.DictReader(open(path)))
    # Infer the prescribed program from the data rather than from run.json:
    # a run that failed validation writes only a PARTIAL record, without the
    # strain program, so metadata is exactly what is missing when it matters.
    # Under displacement control boundary_strain / load_factor is a constant
    # (the target strain) on every step the solver actually solved, so the
    # modal ratio identifies the program and the outliers are the abandoned
    # steps.
    # Prefer the RECORDED target. The modal-ratio inference below assumes
    # most steps are valid, and silently inverts when they are not: a run
    # whose steps mostly failed has a modal ratio of 0, so the zeros are
    # kept as "valid" and the good steps discarded. That produced "20 valid
    # steps, strain to 0.00%" on a run that had genuinely solved 5.
    target = None
    meta_path = path.parent / "run.json"
    if meta_path.exists():
        meta = json.load(open(meta_path))
        target = meta.get("target_strain") or (meta.get("diagnostics") or {}).get(
            "target_strain"
        )
    ratios = [
        float(r["boundary_strain"]) / float(r["load_factor"])
        for r in rows if float(r["load_factor"]) > 0
    ]
    if target is None and ratios:
        target = float(np.median(ratios))
    if target:
        rows = [
            r for r in rows
            if float(r["load_factor"]) > 0
            and abs(float(r["boundary_strain"]) / float(r["load_factor"]) - target)
            <= 1e-3 * target
        ]
    eps = np.asarray([float(r["boundary_strain"]) for r in rows])
    sig = np.asarray([float(r["avg_stress_axial"]) for r in rows])
    order = np.argsort(eps)
    return eps[order], sig[order]


def curve_difference(a: Path, b: Path, n_grid: int = 200) -> dict:
    """Relative L2 difference between two stress-strain curves.

    Compared on the OVERLAPPING strain range only. Extrapolating the shorter
    curve would manufacture agreement or disagreement out of nothing, so a
    pair whose overlap is small is reported as such rather than scored.
    """
    ea, sa = _curve(a)
    eb, sb = _curve(b)
    lo, hi = max(ea.min(), eb.min()), min(ea.max(), eb.max())
    if not (hi > lo):
        return {"overlap": 0.0, "rel_l2": None, "note": "no overlapping strain range"}
    grid = np.linspace(lo, hi, n_grid)
    ia, ib = np.interp(grid, ea, sa), np.interp(grid, eb, sb)
    rel = float(np.sqrt(np.mean((ia - ib) ** 2)) / max(1e-30, np.sqrt(np.mean(ib**2))))
    return {
        "overlap_strain": [float(lo), float(hi)],
        "overlap_fraction_of_coarser": float((hi - lo) / max(1e-30, ea.max() - ea.min())),
        "rel_l2": rel,
        "max_abs_diff_MPa": float(np.max(np.abs(ia - ib))),
        "mean_stress_MPa": float(np.mean(ib)),
    }


def experiment_error(path: Path, markers: Path) -> dict:
    """RMSE of a run against the digitized experimental markers."""
    eps, sig = _curve(path)
    xe, ye = [], []
    with open(markers) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            # The digitized markers store strain in PERCENT ("strain_pct"),
            # while every simulation column is dimensionless strain. Reading
            # the column without converting would place the experiment at
            # 100x the strain and report a meaningless error.
            xe.append(float(row["strain_pct"]) / 100.0)
            ye.append(float(row["stress_MPa"]))
    xe, ye = np.asarray(xe), np.asarray(ye)
    keep = (xe >= eps.min()) & (xe <= eps.max())
    if not keep.any():
        return {"n_markers_in_range": 0, "rmse_MPa": None}
    pred = np.interp(xe[keep], eps, sig)
    return {
        "n_markers_in_range": int(keep.sum()),
        "n_markers_total": int(len(xe)),
        "rmse_MPa": float(np.sqrt(np.mean((pred - ye[keep]) ** 2))),
        "max_abs_MPa": float(np.max(np.abs(pred - ye[keep]))),
        "rel_rmse": float(np.sqrt(np.mean((pred - ye[keep]) ** 2)) / np.mean(np.abs(ye[keep]))),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True,
                    help="density=path/to/results.csv, in increasing density order")
    ap.add_argument("--tol-cp", type=float, default=0.02,
                    help="relative L2 tolerance on the stress-strain curve")
    ap.add_argument("--tol-init", type=float, default=0.005,
                    help="absolute tolerance on initiation strain")
    ap.add_argument("--markers", default="data/digitized/al2219_fig3_experiment_markers.csv")
    ap.add_argument("--out", default="results/convergence.json")
    args = ap.parse_args()

    runs = {}
    for spec in args.runs:
        d, p = spec.split("=", 1)
        runs[int(d)] = STUDY / p
    order = sorted(runs)

    rec: dict = {"densities": order, "pairs": [], "vs_experiment": {}}
    for d in order:
        if runs[d].exists():
            rec["vs_experiment"][d] = experiment_error(runs[d], STUDY / args.markers)

    converged_at = None
    for lo, hi in zip(order, order[1:]):
        if not (runs[lo].exists() and runs[hi].exists()):
            continue
        pair = {"coarse": lo, "fine": hi, **curve_difference(runs[lo], runs[hi])}
        init = {}
        for d in (lo, hi):
            f = runs[d].parent / "initiation.json"
            if f.exists():
                init[d] = json.load(open(f))["boundary_strain"]
        if len(init) == 2:
            pair["initiation_strain"] = {str(k): v for k, v in init.items()}
            pair["initiation_delta"] = abs(init[hi] - init[lo])
            pair["initiation_converged"] = bool(pair["initiation_delta"] <= args.tol_init)
        pair["cp_converged"] = bool(
            pair.get("rel_l2") is not None and pair["rel_l2"] <= args.tol_cp
        )
        rec["pairs"].append(pair)
        if converged_at is None and pair["cp_converged"]:
            converged_at = lo  # the COARSER of the first agreeing pair is usable

    rec["tolerances"] = {"cp_rel_l2": args.tol_cp, "initiation_strain_abs": args.tol_init}
    rec["least_dense_usable"] = converged_at
    rec["note"] = (
        "least_dense_usable is the COARSER member of the first consecutive pair that "
        "agrees within tolerance: agreement with the next finer mesh is what licenses "
        "using it. A null means no pair agreed, so the family does not yet demonstrate "
        "convergence and a finer mesh is required."
    )
    out = STUDY / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1)[:2000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
