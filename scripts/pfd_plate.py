#!/usr/bin/env python3
"""Perforated-plate PFD pipeline: pipeline deck -> PINNACLE deck -> run -> curve.

One pipeline deck (PINNACLE ``.i`` syntax) describes a J2+PFD or CP+PFD plate
case; the script renders the matching template under ``decks/pfd_study/
templates/``, runs it and reduces the history to a results file.

  [Plate]
    out_dir      = ../../results/pfd_study/j2/cal_l5_gc1   # the PINNACLE run lands in <out_dir>/solver
    mechanics    = j2                        # j2 | cp
    mesh         = ../../results/pfd_study/mesh/h40_hole5/mesh.e
    points       = ../../results/cp_study/cal/gn2/base/orientations_points.csv   # cp only
    width = 100
    height = 50                              # gross section = height x thickness
    thickness = 1
    strain       = 0.20                      # nominal strain at the end of the path
    t_end        = 200.0                     # physical time of the path (CP rate, damage dt)
    n_steps      = 256                       # staggered load steps
    substeps     = 1                         # mechanics sub-steps per load step
    max_iterations = 3                       # staggered alternations per step ...
    tolerance = 1.0e-4                       # ... and their test
    snapshots_every = 4
    stop_damage = 0.95                       # end the path once a node reaches this damage (0 = full path)
    backend = torch
    use_gpu = false
    threads = 2
    reaction_set = x+                        # hexaform plate sets
    damage_set = grips
    newton_tol = 1.0e-3
    newton_rel_tol = 0
    max_iter = 50
    max_cutbacks = 10                        # load-increment halvings on a Newton failure
    line_search = backtracking               # Newton globalisation (none | backtracking)
    sparse_threshold = 500
    damage_threshold = 0.8
    # J2 (mechanics = j2): either the power law or a piecewise table from j2_fit.json
    j2_fit       = ../../results/pfd_study/j2_fit/j2_fit.json
    j2_law       = piecewise                 # power | piecewise (law of j2_fit.json)
    finite_strain = true
    # CP (mechanics = cp): the Part I constants, as in decks/cp_study/cal_center2.i
    C11 = 108.0e3 ... y_c = 0.00343758        # one key per line: C11 C12 C44 mu gamma0 nexp tau_c0 q_lat alpha0 burgers rho0 k1 y_c
    # PFD
    Gc = 10.0
    l0 = 5.0
    eta = 0.0
    zeta = 1.0
    W_c = 0.0
    damage_model = at2
  []

Stages:
  deck    render <out_dir>/pfd.i (paths relative to it, so the deck runs on any host
          with the same tree) and write <out_dir>/case.json;
  run     `pinnacle run pfd.i` in <out_dir> with OMP/MKL threads set (skipped when
          solver/solver_history.csv exists);
  curve   <out_dir>/results.csv from solver/solver_history.csv: nominal strain
          (load factor x strain), engineering stress (reaction_x over the gross
          section), the volume-mean stress, max damage, initiation counters;
          <out_dir>/fracture.json: fracture strain (first step with a node at or
          above damage_threshold, and the element-mean criterion), peak stress,
          fracture energy (integral of the curve to the end), steps and iterations.

Usage: pfd_plate.py PIPELINE.i deck|run|curve|all
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = {
    "j2": STUDY_ROOT / "decks/pfd_study/templates/j2_pfd_plate_template.i",
    "cp": STUDY_ROOT / "decks/pfd_study/templates/cp_pfd_plate_template.i",
}
PINNACLE = os.environ.get("PINNACLE_BIN", "pinnacle")
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")


def read_pipeline(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "plate" not in cfg:
        raise SystemExit(f"{path}: missing [Plate] block")
    return cfg["plate"]


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _num(v) -> str:
    return f"{float(v):.12g}"


def _bool(v) -> str:
    return "true" if str(v).lower() in ("true", "1", "yes") else "false"


def _list(vals) -> str:
    return "[" + ", ".join(f"{float(v):.10g}" for v in vals) + "]"


def stage_deck(pl: dict, pipe: Path) -> Path:
    mech = str(pl.get("mechanics", "j2")).lower()
    if mech not in TEMPLATES:
        raise SystemExit(f"mechanics = {mech}: use j2 or cp")
    out = _resolve(str(pl["out_dir"]), pipe)
    out.mkdir(parents=True, exist_ok=True)
    rel = lambda q: os.path.relpath(q, out)  # noqa: E731  # host-independent deck
    n_steps = int(pl.get("n_steps", 256))
    t_end = float(pl.get("t_end", 200.0))
    strain = float(pl.get("strain", 0.2))
    width = float(pl.get("width", 100.0))
    sub = {
        "BACKEND": str(pl.get("backend", "torch")),
        "USEGPU": _bool(pl.get("use_gpu", False)),
        "OUTDIR": "solver",
        "MESH": rel(_resolve(str(pl["mesh"]), pipe)),
        "STRAIN": _num(strain),
        "TEND": _num(t_end),
        "PULL": _num(strain * width),
        "NSTEPS": str(n_steps),
        "SUBSTEPS": str(int(pl.get("substeps", 1))),
        "MAXALT": str(int(pl.get("max_iterations", 3))),
        "ALTTOL": _num(pl.get("tolerance", 1.0e-4)),
        "SNAPSHOTS": str(int(pl.get("snapshots_every", 0))),
        "STOPDAMAGE": _num(pl.get("stop_damage", 0.0)),
        "MAXITER": str(int(pl.get("max_iter", 50))),
        "NEWTONTOL": _num(pl.get("newton_tol", 1.0e-3)),
        "NEWTONRELTOL": _num(pl.get("newton_rel_tol", 0.0)),
        "MAXCUTBACKS": str(int(pl.get("max_cutbacks", 10))),
        "LINESEARCH": str(pl.get("line_search", "backtracking")),
        "REACTIONSET": str(pl.get("reaction_set", "x+")),
        "DAMAGESET": str(pl.get("damage_set", "grips")),
        "SPARSETHRESHOLD": str(int(pl.get("sparse_threshold", 500))),
        "DTHRESHOLD": _num(pl.get("damage_threshold", 0.8)),
        "DT": _num(t_end / n_steps),
        "GC": _num(pl["Gc"]),
        "L0": _num(pl["l0"]),
        "ETA": _num(pl.get("eta", 0.0)),
        "ZETA": _num(pl.get("zeta", 1.0)),
        "WC": _num(pl.get("W_c", 0.0)),
        "DAMAGEMODEL": str(pl.get("damage_model", "at2")),
        "FINITESTRAIN": _bool(pl.get("finite_strain", True)),
    }
    if mech == "j2":
        fit = json.loads(_resolve(str(pl["j2_fit"]), pipe).read_text())
        law = str(pl.get("j2_law", "piecewise"))
        entry = fit["laws"][law]
        sub.update({"E": _num(fit["E_MPa"]), "NU": _num(fit["nu"]), "HARDENING": law})
        if law == "power":
            sub.update({"SIGMAY0": _num(entry["sigma_y0_MPa"]), "KHARD": _num(entry["K_MPa"]), "NHARD": _num(entry["n"]),
                        "ALPHAPOINTS": "[]", "SIGMAPOINTS": "[]"})
        else:
            sub.update({"SIGMAY0": _num(entry["yield_stress_points"][0]), "KHARD": "0.0", "NHARD": "1.0",
                        "ALPHAPOINTS": _list(entry["plastic_strain_points"]), "SIGMAPOINTS": _list(entry["yield_stress_points"])})
    else:
        sub["POINTFILE"] = rel(_resolve(str(pl["points"]), pipe))
        for key, name in (("C11", "C11"), ("C12", "C12"), ("C44", "C44"), ("MU", "mu"), ("GAMMA0", "gamma0"), ("NEXP", "nexp"),
                          ("TAUC0", "tau_c0"), ("QLAT", "q_lat"), ("ALPHA0", "alpha0"), ("BURGERS", "burgers"),
                          ("RHO0", "rho0"), ("K1", "k1"), ("YC", "y_c")):
            sub[key] = _num(pl[name])
    src = TEMPLATES[mech].read_text()
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    leftover = sorted(set(re.findall(r"@([A-Z0-9_]+)@", re.sub(r"#.*", "", src))))
    if leftover:
        raise SystemExit(f"unsubstituted placeholders: {leftover}")
    (out / "pfd.i").write_text(src)
    (out / "case.json").write_text(json.dumps({"pipeline": str(pipe), "mechanics": mech, "substitutions": sub,
                                              "plate": {k: v for k, v in pl.items()}}, indent=1, default=str))
    print("wrote", out / "pfd.i")
    return out


def stage_run(pl: dict, pipe: Path) -> None:
    out = _resolve(str(pl["out_dir"]), pipe)
    if (out / "solver" / "solver_history.csv").exists():
        print("exists, skipping run:", out / "solver" / "solver_history.csv")
        return
    threads = str(int(pl.get("threads", 2)))
    env = dict(os.environ, OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads)
    with (out / "run.log").open("w") as log:
        rc = subprocess.run([PINNACLE, "run", "pfd.i"], cwd=out, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
    print("pinnacle run exit", rc, "->", out)


def _converged(v) -> bool:
    t = str(v).strip().lower()
    return t == "true" or (t not in ("", "false") and float(t) == 1.0)


def _f(r: dict, k: str, default=np.nan) -> float:
    v = r.get(k)
    return float(v) if v not in (None, "", "nan") else default


def stage_curve(pl: dict, pipe: Path) -> dict:
    out = _resolve(str(pl["out_dir"]), pipe)
    hist = out / "solver" / "solver_history.csv"
    rows = [r for r in csv.DictReader(hist.open()) if r.get("executioner", r.get("kind", "")) in ("", "staggered", "monolithic")]
    rows = [r for r in rows if "step_factor" in r and r["step_factor"] not in ("", "nan")]
    strain = float(pl.get("strain", 0.2))
    area = float(pl.get("height", 50.0)) * float(pl.get("thickness", 1.0))
    thr = float(pl.get("damage_threshold", 0.8))
    res = []
    for r in rows:
        lf = _f(r, "load_factor", _f(r, "step_factor"))
        res.append({
            "step": int(float(r["step"])),
            "nominal_strain": lf * strain,
            "time": _f(r, "time", lf * float(pl.get("t_end", 200.0))),
            "eng_stress": _f(r, "reaction_x") / area,
            "mean_stress_xx": _f(r, "mean_stress_00"),   # PINNACLE indexes tensor entries 00, 01, ..., 22
            "max_damage": _f(r, "max_damage"),
            "max_element_damage": _f(r, "max_element_damage"),
            "nodes_above": int(_f(r, "nodes_above_threshold", 0)),
            "elements_above": int(_f(r, "elements_above_threshold", 0)),
            "alternations": int(_f(r, "iterations", 0)),
            "newton_iterations": int(_f(r, "mech_iterations", 0)),
            "converged": int(_converged(r.get("converged", "1"))),
            "mech_converged": int(_converged(r.get("mech_converged", "1"))),
        })
    with (out / "results.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    e = np.array([x["nominal_strain"] for x in res]); s = np.array([x["eng_stress"] for x in res])
    first_node = next((x for x in res if x["nodes_above"] > 0), None)
    first_elem = next((x for x in res if x["elements_above"] > 0), None)
    ok = np.isfinite(s)
    frac = {
        "fracture_strain_node": first_node["nominal_strain"] if first_node else None,
        "fracture_strain_element": first_elem["nominal_strain"] if first_elem else None,
        "damage_threshold": thr,
        "peak_eng_stress": float(np.nanmax(s)) if ok.any() else None,
        "strain_at_peak": float(e[np.nanargmax(s)]) if ok.any() else None,
        "fracture_energy_MPa": float(np.trapezoid(s[ok], e[ok])) if ok.sum() > 1 else None,
        "final_strain": float(e[-1]), "final_eng_stress": float(s[-1]), "final_max_damage": res[-1]["max_damage"],
        "stopped_early": len(res) < int(pl.get("n_steps", 256)),
        "mechanics_failed": int(any(not x["mech_converged"] for x in res)),
        "steps": len(res), "steps_converged": int(sum(x["converged"] for x in res)),
        "newton_iterations": int(sum(x["newton_iterations"] for x in res)),
        "alternations": int(sum(x["alternations"] for x in res)),
    }
    (out / "fracture.json").write_text(json.dumps(frac, indent=1))
    print(json.dumps(frac))
    return frac


def main(argv):
    if len(argv) < 3:
        print(__doc__); return 2
    pipe = Path(argv[1]).resolve(); pl = read_pipeline(pipe)
    stages = ["deck", "run", "curve"] if argv[2] == "all" else argv[2:]
    for st in stages:
        {"deck": stage_deck, "run": stage_run, "curve": stage_curve}[st](pl, pipe)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
