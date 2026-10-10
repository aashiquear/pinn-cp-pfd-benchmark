#!/usr/bin/env python3
"""FE Gauss-Newton step for the CP calibration, from a finished iteration deck.

Reads a GN study deck (decks/cp_study/cal_gnN.i: [Study] base = the center
pipeline deck, [Cases] = one perturbed parameter per case), builds the
finite-difference Jacobian of the engineering stress at the experiment
markers from the finished FE runs, takes a damped (Levenberg-Marquardt) step
in log-parameter space clipped to the bounds of decks/cp_study/taylor_screen.i,
and writes the next center deck and the next GN deck. Markers below 0.1%
strain are excluded (the FE curve starts at the first load step).

Writes <out_root>/gauss_newton.json with the residuals, Jacobian, step,
predicted RMSE, and the per-marker errors of the center.

Usage: cp_gauss_newton.py decks/cp_study/cal_gn1.i [--factor 1.2]
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
from cp_parameter_study import _markers, _read_study, cases_of  # noqa: E402
from polycrystal_plate import _resolve  # noqa: E402

KEYMAP = {"rho0": "rho0", "k1": "k1", "y_c": "annihilation_diameter", "tau_c0": "tau_c0"}


def _curve(case_dir: Path):
    rows = list(csv.DictReader((case_dir / "results.csv").open()))
    e = np.array([float(r["nominal_strain"]) for r in rows]); s = np.array([float(r["eng_stress"]) for r in rows])
    return e, s


def _deck_value(text: str, key: str) -> float:
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", text, re.M)
    if not m:
        raise SystemExit(f"key {key} not in deck")
    return float(m.group(1))


def main(argv):
    deck = Path(argv[0]).resolve()
    factor = float(argv[argv.index("--factor") + 1]) if "--factor" in argv else 1.2
    cfg = _read_study(deck)
    base = _resolve(str(cfg["study"]["base"]), deck)
    out_root = _resolve(str(cfg["study"]["out_root"]), deck)
    ex_e, ex_s = _markers(_resolve(str(cfg["study"]["markers"]), deck))
    keep = ex_e >= 1.0e-3
    ex_e, ex_s = ex_e[keep], ex_s[keep]
    base_text = base.read_text()
    names = [k for k in ("rho0", "k1", "y_c", "tau_c0")]
    p0 = np.array([_deck_value(base_text, k) for k in names])
    # bounds from the surrogate deck
    from taylor_surrogate import read_deck
    sur = read_deck(STUDY / "decks/cp_study/taylor_screen.i")
    bounds = {}
    for k, v in sur["free"].items():
        lo, hi, _ = (float(t) for t in (v if isinstance(v, list) else str(v).split()))
        bounds[k] = (lo, hi)
    lo = np.log([max(bounds[KEYMAP[k]][0], 1e-9) for k in names]); hi = np.log([bounds[KEYMAP[k]][1] for k in names])

    e0, s0 = _curve(out_root / "base")
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
        dz = np.log(v / p0[j])
        J[:, j] = (f - f0) / dz
        pert[k] = {"value": v, "curve": f.tolist()}
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
    print(f"iteration {n}: center rmse {rmse0:.1f} MPa; predicted after step {rmse_pred:.1f} MPa")
    print("  center:", {k: round(v, 5) for k, v in zip(names, p0)})
    print("  next:  ", {k: round(v, 5) for k, v in zip(names, p1)})
    print("  sensitivity dS/dln p at 2% / 7% / 17% (MPa):")
    for j, k in enumerate(names):
        idx = [int(np.argmin(np.abs(ex_e - t))) for t in (0.02, 0.07, 0.17)]
        print(f"    {k:>7}: " + "  ".join(f"{J[i, j]:7.1f}" for i in idx))
    # next decks
    center_text = base_text
    for k, v in zip(names, p1):
        center_text = re.sub(rf"^(\s*{re.escape(k)}\s*=\s*)\S+", rf"\g<1>{v:.6g}", center_text, count=1, flags=re.M)
    center_text = re.sub(r"^(\s*out_dir\s*=\s*)\S+", rf"\g<1>../../results/cp_study/cal/center{n+1}", center_text, count=1, flags=re.M)
    center_text = re.sub(r"CENTER \(iteration \d+\)", f"CENTER (iteration {n+1})", center_text)
    (STUDY / f"decks/cp_study/cal_center{n+1}.i").write_text(center_text)
    gn = ["# FE Gauss-Newton iteration %d at h40: center from iteration %d's step and one" % (n + 1, n),
          "# multiplicative perturbation per free parameter.", "[Study]",
          f"  base       = cal_center{n+1}.i", f"  out_root   = ../../results/cp_study/cal/gn{n+1}",
          "  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv", "  concurrent = 5", "[]", "[Cases]"]
    for k, v in zip(names, p1):
        f = 1.5 if k == "tau_c0" else factor
        gn += [f"  [./{k}_p]", f"    {k} = {v * f:.6g}", "  []"]
    gn += ["[]", ""]
    (STUDY / f"decks/cp_study/cal_gn{n+1}.i").write_text("\n".join(gn))
    print(f"wrote decks/cp_study/cal_center{n+1}.i and cal_gn{n+1}.i")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
