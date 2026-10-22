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
  deck    render decks/ni/templates/cp_specimen_template.i into <out_dir>/cp.i (+ case.json); with
          Gc in the deck, CP+PFD: decks/ni/templates/cp_pfd_specimen_template.i into <out_dir>/pfd.i
          (Gc, l0, eta = 0, zeta = 1, W_c = 0, damage_model = at2, max_iterations = 3 staggered
          alternations, tolerance = 1e-3, damage_threshold = 0.8, stop_damage, snapshots_every,
          checkpoint_every = 8, as the Al2219 plate pipeline scripts/pfd_plate.py); the run stops at
          through-section failure (PINNACLE stop_separation 'x- x+', separation_damage = 0.85)
  run     `pinnacle run cp.i|pfd.i` in <out_dir> with OMP/MKL threads set (skipped when the history exists)
  curve   results.csv per converged step, up to the first unconverged one: gauge strain (PINNACLE gauge_node_sets, the DIC analogue),
          grip strain, engineering stress = reaction on x+ / mean initial gage section
          (gage width x mean of the five station thicknesses), Newton count; summary.json

Usage: ni_specimen.py SPECIMEN.i deck|run|curve|all
"""
from __future__ import annotations

import csv
import itertools
import json
import os
import pickle
import re
import resource
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = STUDY_ROOT / "decks/ni/templates/cp_specimen_template.i"
TEMPLATE_PFD = STUDY_ROOT / "decks/ni/templates/cp_pfd_specimen_template.i"
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
    if is_pfd(sp):
        n_steps, t_end = int(sp.get("n_steps", 256)), float(sp.get("t_end", 250.0))
        sub.update({
            "SUBSTEPS": str(int(sp.get("substeps", 1))), "MAXALT": str(int(sp.get("max_iterations", 3))),
            "ALTTOL": _num(sp.get("tolerance", 1.0e-3)), "STOPDAMAGE": _num(sp.get("stop_damage", 0.0)),
            "CHECKPOINT": str(int(sp.get("checkpoint_every", 8))), "DTHRESHOLD": _num(sp.get("damage_threshold", 0.8)),
            "DT": _num(t_end / n_steps), "GC": _num(sp["Gc"]), "L0": _num(sp["l0"]), "ETA": _num(sp.get("eta", 0.0)),
            "ZETA": _num(sp.get("zeta", 1.0)), "WC": _num(sp.get("W_c", 0.0)),
            "DAMAGEMODEL": str(sp.get("damage_model", "at2")),
            "SEPDAMAGE": _num(sp.get("separation_damage", 0.85)),
        })
    src = (TEMPLATE_PFD if is_pfd(sp) else TEMPLATE).read_text()
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    leftover = sorted(set(re.findall(r"@([A-Z0-9_]+)@", re.sub(r"#.*", "", src))))
    if leftover:
        raise SystemExit(f"unsubstituted placeholders: {leftover}")
    (out / deck_name(sp)).write_text(src)
    area, t_mean = gage_section(sp, deck)
    (out / "case.json").write_text(json.dumps({"specimen_deck": str(deck), "substitutions": sub, "gage_area_mm2": area,
                                              "mean_thickness_mm": t_mean, "specimen": {k: v for k, v in sp.items()}}, indent=1, default=str))
    print("wrote", out / deck_name(sp))
    return out


def is_pfd(sp: dict) -> bool:
    """CP+PFD when the specimen deck sets a fracture toughness."""
    return "Gc" in sp


def deck_name(sp: dict) -> str:
    return "pfd.i" if is_pfd(sp) else "cp.i"


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
    cmd, limit = _capped([PINNACLE, "run", deck_name(sp)], float(sp.get("memory_cap_gb", 16)))
    with (out / "run.log").open("w") as log:
        rc = subprocess.run(cmd, cwd=out, stdout=log, stderr=subprocess.STDOUT, env=env,
                            preexec_fn=limit).returncode
    (out / "run_exit.txt").write_text(f"{rc}\n")
    # CP+PFD restart checkpoints: the newest is kept when the run died (its curve), all go when it finished
    ckpts = _checkpoints(out)
    for c in (ckpts if (out / "solver" / "summary.json").exists() else ckpts[:-1]):
        c.unlink()
    print("pinnacle run exit", rc, "->", out)


def _checkpoints(out: Path) -> list[Path]:
    """Restart checkpoints of a CP+PFD run, oldest first."""
    return sorted((out / "solver").glob("restart_step*.pkl"))


def _converged(v) -> bool:
    t = str(v).strip().lower()
    return t == "true" or (t not in ("", "false") and float(t) == 1.0)


def _f(r: dict, k: str, default=np.nan) -> float:
    v = r.get(k)
    return float(v) if v not in (None, "", "nan") else default


def stage_curve(sp: dict, deck: Path) -> dict:
    if is_pfd(sp):
        return _curve_pfd(sp, deck)
    out = _resolve(str(sp["out_dir"]), deck)
    # Converged steps up to the first unconverged one: a recovered cutback leaves
    # no row, so an unconverged row means the run-total cutback budget was spent;
    # the executioner then moves on from a state that is not a solution, and a
    # step that converges after it is not part of the curve.
    rows = list(itertools.takewhile(lambda r: _converged(r.get("converged", "1")),
                                    csv.DictReader((out / "solver" / "solver_history.csv").open())))
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


def _curve_pfd(sp: dict, deck: Path) -> dict:
    """results.csv of a CP+PFD run (one row per staggered load step), from the
    step history, or from the newest restart checkpoint when the run died."""
    out = _resolve(str(sp["out_dir"]), deck)
    died = not (out / "solver" / "summary.json").exists()
    if died:
        ckpts = _checkpoints(out)
        if not ckpts:
            raise SystemExit(f"{out}: the run did not finish and left no checkpoint")
        with ckpts[-1].open("rb") as fh:
            rows = [{k: ("" if v is None else str(float(v)) if isinstance(v, (bool, int, float, np.floating, np.integer)) else str(v))
                     for k, v in r.items()} for r in pickle.load(fh)["history"]]
        print("run died: curve from", ckpts[-1].name)
    else:
        rows = [r for r in csv.DictReader((out / "solver" / "solver_history.csv").open())
                if r.get("executioner", r.get("kind", "")) in ("", "staggered", "monolithic")]
    rows = [r for r in rows if r.get("step_factor") not in (None, "", "nan")]
    area, t_mean = gage_section(sp, deck)
    strain = float(sp.get("strain", 0.25)); t_end = float(sp.get("t_end", 250.0))
    res = []
    for r in rows:
        lf = _f(r, "load_factor", _f(r, "step_factor"))
        res.append({"step": int(float(r["step"])), "grip_strain": lf * strain, "gauge_strain": _f(r, "gauge_strain"),
                    "time": _f(r, "time", lf * t_end), "eng_stress": _f(r, "reaction_x") / area,
                    "max_damage": _f(r, "max_damage"), "nodes_above": int(_f(r, "nodes_above_threshold", 0)),
                    "elements_above": int(_f(r, "elements_above_threshold", 0)),
                    "alternations": int(_f(r, "iterations", 0)), "newton_iterations": int(_f(r, "mech_iterations_total", _f(r, "mech_iterations", 0))),
                    "cutbacks": int(_f(r, "cutbacks", 0)), "separated": int(_f(r, "separated", 0)),
                    "converged": int(_converged(r.get("converged", "1"))), "mech_converged": int(_converged(r.get("mech_converged", "1")))})
    with (out / "results.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    # the curve is valid up to the first step whose mechanics did not converge
    valid = list(itertools.takewhile(lambda x: x["mech_converged"], res))
    e = np.array([x["gauge_strain"] for x in valid]); s = np.array([x["eng_stress"] for x in valid])
    first_node = next((x for x in valid if x["nodes_above"] > 0), None)
    # Failure: the first separated step (PINNACLE stop_separation). Separation can come
    # in the same step whose mechanics stops converging (the section loses its load);
    # that step has no gauge strain, so the failure strain is the last converged one.
    sep = next((x for x in res if x["separated"]), None)
    sep_strain = None
    if sep is not None:
        sep_strain = sep["gauge_strain"] if np.isfinite(sep["gauge_strain"]) else (
            valid[-1]["gauge_strain"] if valid else None)
    ok = np.isfinite(s)
    summary = {"steps": len(res), "steps_valid": len(valid), "run_died": died,
               "final_gauge_strain": float(e[-1]) if len(e) else None, "final_grip_strain": valid[-1]["grip_strain"] if valid else None,
               "peak_eng_stress": float(np.nanmax(s)) if ok.any() else None,
               "strain_at_peak": float(e[int(np.nanargmax(s))]) if ok.any() else None,
               "initiation_gauge_strain_node": first_node["gauge_strain"] if first_node else None,
               "failure_gauge_strain": sep_strain,   # PINNACLE stop_separation
               "failure_step": sep["step"] if sep else None,
               "damage_threshold": float(sp.get("damage_threshold", 0.8)),
               "final_max_damage": valid[-1]["max_damage"] if valid else None,
               "gage_area_mm2": area, "mean_thickness_mm": t_mean,
               "newton_iterations": int(sum(x["newton_iterations"] for x in res)),
               "alternations": int(sum(x["alternations"] for x in res))}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return summary


def stage_fracture(sp: dict, deck: Path) -> dict:
    """The paper's failure criterion on the damage snapshots: the first snapshot in
    which the nodes with damage above ``separation_damage`` (0.85, the paper's c >
    0.85 across the cross-section) cut the specimen in two, i.e. no path of intact
    nodes (element edges) joins the x- grip face to the x+ grip face. Also the
    crack's location along the gage (mean x of the cut nodes, as a fraction of the
    gage from gage- to gage+) and its angle to the loading axis in the x-y plane
    (principal direction of the cut nodes). Writes <out_dir>/fracture.json."""
    from vtk_legacy import read_vtk
    out = _resolve(str(sp["out_dir"]), deck)
    thr = float(sp.get("separation_damage", 0.85))
    pvd = out / "solver" / "fields.pvd"
    frames = [(float(t), out / "solver" / f) for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', pvd.read_text())]
    # The .pvd stamps a coupled run's snapshots with their ordinal k: snapshot k
    # is load step (k + 1) * snapshots_every - 1, and a run stopped by
    # stop_separation adds one more at its failure step (summary.json).
    gauge_of = {int(r["step"]): float(r["gauge_strain"]) for r in _rows(out / "results.csv")}
    every = max(int(sp.get("snapshots_every", 1)), 1)
    summ = json.loads((out / "summary.json").read_text()) if (out / "summary.json").exists() else {}
    fstep = summ.get("failure_step")

    def step_of(k: int, n: int) -> int:
        if fstep is not None and k == n - 1 and (k + 1) * every - 1 != int(fstep):
            return int(fstep)
        return (k + 1) * every - 1
    edges = None
    found = None
    for kf, (_t, path) in enumerate(frames):
        v = read_vtk(path)
        d = np.asarray(v["point"]["damage"], dtype=float).reshape(-1)
        X = v["nodes"]
        if edges is None:
            cells = np.asarray(v["cells"], dtype=int)
            hexe = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
            edges = np.unique(np.sort(np.concatenate([cells[:, [a, b]] for a, b in hexe]), axis=1), axis=0)
            x_lo, x_hi = X[:, 0].min(), X[:, 0].max()
            left = np.where(np.isclose(X[:, 0], x_lo))[0]; right = set(np.where(np.isclose(X[:, 0], x_hi))[0].tolist())
        intact = d <= thr
        adj: dict[int, list[int]] = {}
        for a, b in edges:
            if intact[a] and intact[b]:
                adj.setdefault(int(a), []).append(int(b)); adj.setdefault(int(b), []).append(int(a))
        seen = {int(k) for k in left if intact[k]}; stack = list(seen); reached = False
        while stack and not reached:
            k = stack.pop()
            for m in adj.get(k, ()):
                if m not in seen:
                    if m in right:
                        reached = True; break
                    seen.add(m); stack.append(m)
        if not reached:
            cut = np.where(~intact)[0]
            found = (step_of(kf, len(frames)), path, X[cut])
            break
    res = {"separation_damage": thr, "separated": found is not None}
    if found is not None:
        step, path, P = found
        gx = _gage_x(sp, deck)
        c = P[:, :2] - P[:, :2].mean(axis=0)
        _, vec = np.linalg.eigh(c.T @ c)
        crack_dir = vec[:, -1]
        res.update({"step": step, "snapshot": path.name,
                    "failure_gauge_strain": gauge_of.get(step) if np.isfinite(gauge_of.get(step, np.nan))
                    else summ.get("failure_gauge_strain"),
                    "location_fraction": float((P[:, 0].mean() - gx[0]) / (gx[1] - gx[0])),
                    "angle_to_axis_deg": float(np.degrees(np.arccos(abs(crack_dir[0])))),
                    "cut_nodes": int(len(P))})
    (out / "fracture.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res))
    return res


def _gage_x(sp: dict, deck: Path) -> tuple[float, float]:
    """x (mm) of the first and last thickness station (the gage- / gage+ planes): the
    gage is centered on the bar (hexaform dogbone, mesh.json: length, gage_length)."""
    meta = json.loads(_resolve(str(sp["mesh"]), deck).with_name("mesh.json").read_text())
    x0 = 0.5 * (float(meta["length"]) - float(meta["gage_length"]))
    return x0, x0 + float(meta["gage_length"])


def main(argv):
    if len(argv) < 3:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); sp = read_specimen(deck)
    stages = ["deck", "run", "curve"] if argv[2] == "all" else argv[2:]
    for st in stages:
        {"deck": stage_deck, "run": stage_run, "curve": stage_curve, "fracture": stage_fracture}[st](sp, deck)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
