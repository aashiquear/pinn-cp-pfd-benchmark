#!/usr/bin/env python3
"""Ni single-crystal dog-bone specimen pipeline: specimen deck -> PINNACLE deck -> run -> curve.

Every stage is driven by ONE specimen deck (PINNACLE ``.i`` syntax):

  [Specimen]
    out_dir         = ../../../results/ni_study/cp/check_L1_h040
    mesh            = ../../../results/ni_study/mesh/h040/L.1/mesh.e
    orientations    = ../../../data/extracted/ni_table1_orientations.csv
    group           = L                      # orientation group of Table 1 (Bunge angles)
    thickness       = ../../../data/extracted/ni_table2_thickness.csv
    specimen        = L.1                    # row of Table 2: the mean initial gage section
    gage_width      = 1.85                   # mm (data/extracted/ni_fig2_planform.csv)
    length          = 25.4                   # mm, the pulled bar length
    strain          = 0.25                   # grip strain at t_end (ramp over the bar length)
    t_end           = 250.0                  # s: the grip strain rate strain / t_end
    n_steps         = 256
    snapshots_every = 0
    backend         = torch     use_gpu = false     threads = 2     memory_cap_gb = 16   # resident-memory cap per run, cgroup (0 = none)
    newton_tol = 1.0e-3   newton_rel_tol = 0   max_cutbacks = 10   linear_method = direct   sparse_threshold = 500
    C11 = 253000  C12 = 152000  C44 = 124000  mu = 26000  gamma0 = 1.0e-5  nexp = 8.3333  tau_c0 = 0.26
    q_lat = 1.4  alpha0 = 0.4  burgers = 2.5e-7  rho0 = 2.35e6  k1 = 0.0295  y_c = 0.79e-6
    predictor = none | linear | fem_solver_grain   (corrector, corrector_tol, predictor_device, predictor_gpu as in pfd_plate.py)
    state_update = implicit   state_iterations = 2   # PINNACLE ADR 0033: implicit (backward Euler) is the default here; explicit = Parts I-II scheme

Stages:
  deck    render decks/ni/templates/cp_specimen_template.i into <out_dir>/cp.i (+ case.json)
  run     `pinnacle run cp.i` in <out_dir> with OMP/MKL threads set (skipped when the history exists)
  curve   results.csv per accepted step: gauge strain (PINNACLE gauge_node_sets, the DIC analogue),
          grip strain, engineering stress = reaction on x+ / mean initial gage section
          (gage width x mean of the five station thicknesses), Newton count; summary.json

Usage: ni_specimen.py SPECIMEN.i deck|run|curve|all
"""
from __future__ import annotations

import csv
import json
import os
import re
import resource
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = STUDY_ROOT / "decks/ni/templates/cp_specimen_template.i"
PINNACLE = os.environ.get("PINNACLE_BIN", "pinnacle")
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")

CP_KEYS = (("C11", "C11"), ("C12", "C12"), ("C44", "C44"), ("MU", "mu"), ("GAMMA0", "gamma0"), ("NEXP", "nexp"),
           ("TAUC0", "tau_c0"), ("QLAT", "q_lat"), ("ALPHA0", "alpha0"), ("BURGERS", "burgers"),
           ("RHO0", "rho0"), ("K1", "k1"), ("YC", "y_c"))


def read_specimen(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "specimen" not in cfg:
        raise SystemExit(f"{path}: missing [Specimen] block")
    return cfg["specimen"]


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _num(v) -> str:
    return f"{float(v):.12g}"


def _bool(v) -> str:
    return "true" if str(v).lower() in ("true", "1", "yes") else "false"


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))


def gage_section(sp: dict, deck: Path) -> tuple[float, float]:
    """Mean initial gage section (mm^2) and mean thickness of the specimen's Table 2 row."""
    row = next(r for r in _rows(_resolve(str(sp["thickness"]), deck)) if r["specimen"] == str(sp["specimen"]))
    t = np.mean([float(row[f"t{k}_mm"]) for k in range(1, 6)])
    return float(sp.get("gage_width", 1.85)) * float(t), float(t)


def euler_of(sp: dict, deck: Path) -> str:
    row = next(r for r in _rows(_resolve(str(sp["orientations"]), deck)) if r["group"] == str(sp["group"]))
    return f"{float(row['phi1_deg']):g} {float(row['Phi_deg']):g} {float(row['phi2_deg']):g}"


def stage_deck(sp: dict, deck: Path) -> Path:
    out = _resolve(str(sp["out_dir"]), deck); out.mkdir(parents=True, exist_ok=True)
    rel = lambda q: os.path.relpath(q, out)  # noqa: E731
    sub = {
        "BACKEND": str(sp.get("backend", "torch")), "USEGPU": _bool(sp.get("use_gpu", False)), "OUTDIR": "solver",
        "MESH": rel(_resolve(str(sp["mesh"]), deck)), "STRAIN": _num(sp.get("strain", 0.25)), "LENGTH": _num(sp.get("length", 25.4)),
        "TEND": _num(sp.get("t_end", 250.0)), "NSTEPS": str(int(sp.get("n_steps", 256))),
        "NEWTONTOL": _num(sp.get("newton_tol", 1.0e-3)), "NEWTONRELTOL": _num(sp.get("newton_rel_tol", 0.0)),
        "MAXCUTBACKS": str(int(sp.get("max_cutbacks", 10))), "LINEAR": str(sp.get("linear_method", "direct")),
        "SPARSETHRESHOLD": str(int(sp.get("sparse_threshold", 500))), "SNAPSHOTS": str(int(sp.get("snapshots_every", 0))),
        "EULER": euler_of(sp, deck), "CORRECTOR": str(sp.get("corrector", "newton")), "CORRECTORTOL": _num(sp.get("corrector_tol", 0.05)),
        "REJECTRATIO": _num(sp.get("predictor_reject_ratio", 10)),
        "STATEUPDATE": str(sp.get("state_update", "implicit")),
        "STATEITER": str(int(sp.get("state_iterations", 2))),
    }
    for key, name in CP_KEYS:
        sub[key] = _num(sp[name])
    pred = str(sp.get("predictor", "none")).lower()
    if pred == "none":
        sub["UTILITIES"] = ""
    elif pred == "linear":
        sub["UTILITIES"] = "[Utilities]\n  [./predictor]\n    kind         = load_step_predictor\n    mode         = linear\n    warmup_steps = 2\n  []\n[]\n"
    elif pred == "fem_solver_grain":
        sub["UTILITIES"] = "\n".join([
            "[Utilities]", "  [./predictor]", "    kind             = load_step_predictor", "    mode             = fem_residual",
            "    baseline         = linear", "    warmup_steps     = 2", "    backend          = torch",
            f"    device           = {sp.get('predictor_device', 'auto')}", f"    use_gpu          = {_bool(sp.get('predictor_gpu', True))}",
            "    precision        = float64", "    hidden_layers    = 3", "    neurons          = 64", "    activation       = tanh",
            "    fourier_features = false", "    epochs_initial   = 5", "    epochs_per_step  = 4", "    lr               = 1",
            "    optimizer        = lbfgs", f"    residual_target  = {_num(sp.get('predictor_target', 0.05))}", "    reset_each_step  = false",
            "    seed             = 0", "    node_features    = 'orientation grain_boundary'", "  []", "[]", ""])
    else:
        raise SystemExit(f"predictor = {pred}: use none, linear or fem_solver_grain")
    src = TEMPLATE.read_text()
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    leftover = sorted(set(re.findall(r"@([A-Z0-9_]+)@", re.sub(r"#.*", "", src))))
    if leftover:
        raise SystemExit(f"unsubstituted placeholders: {leftover}")
    (out / "cp.i").write_text(src)
    area, t_mean = gage_section(sp, deck)
    (out / "case.json").write_text(json.dumps({"specimen_deck": str(deck), "substitutions": sub, "gage_area_mm2": area,
                                              "mean_thickness_mm": t_mean, "specimen": {k: v for k, v in sp.items()}}, indent=1, default=str))
    print("wrote", out / "cp.i")
    return out


def _capped(cmd: list, cap_gb: float):
    """(command, preexec_fn) running ``cmd`` under a resident-memory cap of
    ``cap_gb`` GB: a user cgroup scope (MemoryMax, no swap; the kernel kills
    the run at the cap) when systemd-run works, else an address-space rlimit."""
    if cap_gb <= 0:
        return cmd, None
    if shutil.which("systemd-run") and subprocess.run(
            ["systemd-run", "--user", "--scope", "-q", "true"], capture_output=True, check=False).returncode == 0:
        return ["systemd-run", "--user", "--scope", "-q", f"-pMemoryMax={cap_gb:g}G", "-pMemorySwapMax=0"] + cmd, None
    cap = int(cap_gb * 1024**3)
    return cmd, lambda: resource.setrlimit(resource.RLIMIT_AS, (cap, cap))


def stage_run(sp: dict, deck: Path) -> None:
    out = _resolve(str(sp["out_dir"]), deck)
    # solver/summary.json is written only by a finished run; the per-step
    # history is rewritten every step, so it does not mean the run is done.
    if (out / "solver" / "summary.json").exists():
        print("exists, skipping run:", out / "solver" / "summary.json"); return
    threads = str(int(sp.get("threads", 2)))
    env = dict(os.environ, OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads, PYTHONFAULTHANDLER="1")
    # Memory cap per run (memory_cap_gb, 0 = none): a runaway allocation after
    # a backward-Euler fallback once took 52 GB and froze the host (docs/10
    # §10). The cap is on resident memory (a cgroup scope, systemd-run): an
    # address-space cap counts the ~11 GB the ROCm torch libraries map but
    # barely touch, and failed runs at 6 GB resident (docs/12 §4.3).
    cmd, limit = _capped([PINNACLE, "run", "cp.i"], float(sp.get("memory_cap_gb", 16)))
    with (out / "run.log").open("w") as log:
        rc = subprocess.run(cmd, cwd=out, stdout=log, stderr=subprocess.STDOUT, env=env,
                            preexec_fn=limit).returncode
    (out / "run_exit.txt").write_text(f"{rc}\n")
    print("pinnacle run exit", rc, "->", out)


def _converged(v) -> bool:
    t = str(v).strip().lower()
    return t == "true" or (t not in ("", "false") and float(t) == 1.0)


def _f(r: dict, k: str, default=np.nan) -> float:
    v = r.get(k)
    return float(v) if v not in (None, "", "nan") else default


def stage_curve(sp: dict, deck: Path) -> dict:
    out = _resolve(str(sp["out_dir"]), deck)
    rows = [r for r in csv.DictReader((out / "solver" / "solver_history.csv").open()) if _converged(r.get("converged", "1"))]
    area, t_mean = gage_section(sp, deck)
    strain = float(sp.get("strain", 0.25)); t_end = float(sp.get("t_end", 250.0))
    res = []
    for r in rows:
        lf = _f(r, "load_factor")
        res.append({"step": int(float(r["step"])), "grip_strain": lf * strain, "gauge_strain": _f(r, "gauge_strain"),
                    "time": _f(r, "time", lf * t_end), "eng_stress": _f(r, "reaction_x") / area,
                    "mean_stress_xx": _f(r, "mean_stress_00"), "newton_iterations": int(_f(r, "iterations", 0)),
                    "cutbacks": int(_f(r, "cutbacks", 0))})
    with (out / "results.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    e = np.array([x["gauge_strain"] for x in res]); s = np.array([x["eng_stress"] for x in res])
    summary = {"steps": len(res), "final_gauge_strain": float(e[-1]), "final_grip_strain": res[-1]["grip_strain"],
               "gauge_to_grip_ratio": float(e[-1] / res[-1]["grip_strain"]) if res[-1]["grip_strain"] else None,
               "peak_eng_stress": float(np.nanmax(s)), "strain_at_peak": float(e[int(np.nanargmax(s))]),
               "gage_area_mm2": area, "mean_thickness_mm": t_mean,
               "newton_iterations": int(sum(x["newton_iterations"] for x in res))}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return summary


def main(argv):
    if len(argv) < 3:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); sp = read_specimen(deck)
    stages = ["deck", "run", "curve"] if argv[2] == "all" else argv[2:]
    for st in stages:
        {"deck": stage_deck, "run": stage_run, "curve": stage_curve}[st](sp, deck)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
