#!/usr/bin/env python3
"""FE Gauss-Newton step of the Ni CP calibration, from a finished iteration deck.

Reads a GN study deck (decks/ni/cp/cal_gnN.i: [Study] base = the centre
specimen deck, markers = the Fig. 3 experiment markers, [Cases] = one
perturbed parameter per case; run by scripts/ni_study.py), builds the
finite-difference Jacobian of the engineering stress (against the GAUGE
strain, the DIC analogue) at the markers from the finished FE runs, takes a
damped (Levenberg-Marquardt) step in log-parameter space clipped to the
bounds below, and writes the next centre deck and the next GN deck. Markers
below 0.3% strain (the elastic rise) are excluded.

Bounds (log space, x0.1 .. x10 of the paper's set) keep the step inside the
range where the law's reading is the paper's. Writes
<out_root>/gauss_newton.json with the residuals, Jacobian, step, predicted
RMSE and the per-marker errors of the centre.

Usage: ni_gauss_newton.py decks/ni/cp/cal_gn1.i [--factor 1.2]
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(STUDY / "scripts"))
from ni_study import _resolve, cases_of, read_study  # noqa: E402

PAPER = {"rho0": 2.35e6, "k1": 0.0295, "y_c": 0.79e-6, "tau_c0": 0.26}
BOUNDS = {k: (0.1 * v, 10.0 * v) for k, v in PAPER.items()}


def _markers(path: Path):
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    return np.array([float(r["strain"]) for r in rows]), np.array([float(r["stress_MPa"]) for r in rows])


def _curve(case_dir: Path):
    rows = list(csv.DictReader((case_dir / "results.csv").open()))
    e = np.array([float(r["gauge_strain"]) for r in rows]); s = np.array([float(r["eng_stress"]) for r in rows])
    return e, s


def _deck_value(text: str, key: str) -> float:
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", text, re.M)
    if not m:
        raise SystemExit(f"key {key} not in deck")
    return float(m.group(1))


def main(argv):
    deck = Path(argv[0]).resolve()
    factor = float(argv[argv.index("--factor") + 1]) if "--factor" in argv else 1.2
    cfg = read_study(deck)
    base = _resolve(str(cfg["study"]["base"]), deck)
    out_root = _resolve(str(cfg["study"]["out_root"]), deck)
    ex_e, ex_s = _markers(_resolve(str(cfg["study"]["markers"]), deck))
    keep = ex_e >= 3.0e-3
    ex_e, ex_s = ex_e[keep], ex_s[keep]
    base_text = base.read_text()
    names = ["rho0", "k1", "y_c", "tau_c0"]
    p0 = np.array([_deck_value(base_text, k) for k in names])
    lo = np.log([BOUNDS[k][0] for k in names]); hi = np.log([BOUNDS[k][1] for k in names])
    e0, s0 = _curve(out_root / "base")
    inside = ex_e <= e0.max()
    if not inside.all():
        print("markers beyond the FE path dropped:", ex_e[~inside].tolist())
    ex_e, ex_s = ex_e[inside], ex_s[inside]
    f0 = np.interp(ex_e, e0, s0); r0 = f0 - ex_s
    J = np.zeros((len(ex_e), len(names)))
    pert = {}
    for name, params in cases_of(cfg):
        if name == "base":
            continue
        (k, v), = params.items()
        j = names.index(k)
        e, s = _curve(out_root / name)
        f = np.interp(ex_e, e, s)
        dz = np.log(float(v) / p0[j])
        J[:, j] = (f - f0) / dz
        pert[k] = {"value": float(v), "curve": f.tolist()}
    rmse0 = float(np.sqrt(np.mean(r0 ** 2)))
    z = np.log(p0)
    lam = 1e-3 * np.trace(J.T @ J) / len(z)
    step = -np.linalg.solve(J.T @ J + lam * np.eye(len(z)), J.T @ r0)
    step = np.clip(step, -0.7, 0.7)
    z1 = np.clip(z + step, lo, hi)
    p1 = np.exp(z1)
    pred = f0 + J @ (z1 - z)
    rmse_pred = float(np.sqrt(np.mean((pred - ex_s) ** 2)))
    n = int(re.search(r"gn(\d+)", deck.stem).group(1))
    rec = {"iteration": n, "center": dict(zip(names, p0.tolist())), "rmse_center_MPa": rmse0,
           "per_marker_center": [{"strain": float(a), "experiment": float(b), "fe": float(c), "rel_err": float((c - b) / b)}
                                 for a, b, c in zip(ex_e, ex_s, f0)],
           "jacobian_dS_dlnp": J.tolist(), "perturbations": pert, "step_ln": (z1 - z).tolist(),
           "next": dict(zip(names, p1.tolist())), "rmse_predicted_linear_MPa": rmse_pred}
    (out_root / "gauss_newton.json").write_text(json.dumps(rec, indent=1))
    print(f"iteration {n}: centre rmse {rmse0:.1f} MPa; predicted after step {rmse_pred:.1f} MPa")
    print("  centre:", {k: float(f"{v:.5g}") for k, v in zip(names, p0)})
    print("  next:  ", {k: float(f"{v:.5g}") for k, v in zip(names, p1)})
    print("  sensitivity dS/dln p at 2% / 10% / 30% (MPa):")
    for j, k in enumerate(names):
        idx = [int(np.argmin(np.abs(ex_e - t))) for t in (0.02, 0.10, 0.30)]
        print(f"    {k:>7}: " + "  ".join(f"{J[i, j]:7.1f}" for i in idx))
    center_text = base_text
    for k, v in zip(names, p1):
        center_text = re.sub(rf"^(\s*{re.escape(k)}\s*=\s*)\S+", rf"\g<1>{v:.6g}", center_text, count=1, flags=re.M)
    center_text = re.sub(r"^(\s*out_dir\s*=\s*)\S+", rf"\g<1>../../../results/ni_study/cal/center{n+1}", center_text, count=1, flags=re.M)
    center_text = re.sub(r"^# Ni single crystal, specimen L.1.*$", f"# Ni single crystal, specimen L.1: calibration CENTRE (iteration {n+1}, from gn{n}'s step).", center_text, count=1, flags=re.M)
    (STUDY / f"decks/ni/cp/cal_center{n+1}.i").write_text(center_text)
    gn = [f"# Ni CP calibration, Gauss-Newton iteration {n+1}: centre from iteration {n}'s step and one",
          "# multiplicative perturbation per free parameter (scripts/ni_gauss_newton.py).", "[Study]",
          f"  base       = cal_center{n+1}.i", f"  out_root   = ../../../results/ni_study/cal/gn{n+1}",
          "  markers    = ../../../data/digitized/ni_fig3_experiment_markers.csv", "  concurrent = 5", "[]", "[Cases]"]
    for k, v in zip(names, p1):
        f = 1.5 if k == "tau_c0" else factor
        gn += [f"  [./{k}_p]", f"    {k} = {v * f:.6g}", "  []"]
    gn += ["[]", ""]
    (STUDY / f"decks/ni/cp/cal_gn{n+1}.i").write_text("\n".join(gn))
    print(f"wrote decks/ni/cp/cal_center{n+1}.i and cal_gn{n+1}.i")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
