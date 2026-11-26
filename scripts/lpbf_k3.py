#!/usr/bin/env python3
"""Part V, case K3: tensile test of a polycrystal cube and its crystal-plasticity calibration.

usage: scripts/lpbf_k3.py run <deck.i> [key=value ...]
       scripts/lpbf_k3.py calibrate <start deck.i> <target.csv> [iterations] [tag]
       scripts/lpbf_k3.py figure <target.csv> <out.png> <label=run_dir> ...
       scripts/lpbf_k3.py cube <column run dir> <out.csv> <x0> <y0> <z0> <edge>

  run        run the deck (PINNACLE's cp_quasi_static through the Runner) and write, into
             its output_dir, curve.csv (time, engineering strain, mean Cauchy stress along
             the load axis, engineering stress), fields.npz (final displacement, element
             stress and plastic strain, snapshots when the deck keeps them), wall.json.
             Engineering stress = mean Cauchy stress / (1 + engineering strain): the
             cross-section shrinks at constant volume (the elastic volume change, about
             0.2 %, is neglected).
  calibrate  Gauss-Newton fit of tau_c0, k1 and annihilation_diameter (in the logarithm
             of the parameters) to the target engineering curve between 1 % strain and
             the end of the run. Each iteration runs the center deck and one deck per
             parameter raised by 5 % (four processes at once); the Jacobian is their
             finite difference; the step
             is damped (Levenberg-Marquardt, lambda = 1e-3 of the diagonal) and limited
             to a factor of 2 per parameter. Decks cal<tag>_<i>_<name>.i are written next
             to the start deck, results under <output_dir of the start deck>/../cal<tag>/.
             Stops when the root-mean-square error changes by less than 2 % of itself.
  figure     engineering curves of the given runs against the target (mean and the
             standard deviation between the specimens).

  cube       cut a cube out of a build-column run of scripts/lpbf_column.py and write
             its grain structure as an orientation point file (one row per node of the
             column inside the cube: position relative to the cube corner, Bunge angles
             of the order parameter that is largest there, grain id = that order
             parameter). Nodes that are not solid (sum of eta^2 below 0.25) are left out;
             the tensile deck takes the nearest listed point.

PINNACLE must be importable (PYTHONPATH).
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

FITTED = ("tau_c0", "k1", "annihilation_diameter")
LOAD_COMPONENT = 2


def load_deck(deck: Path, overrides: list[str]):
    from pinnacle.io.loader import load_input

    cfg = load_input(str(deck))
    applied = {}
    for item in overrides:
        dotted, value = item.split("=", 1)
        node = cfg
        parts = dotted.split(".")
        for part in parts[:-1]:
            node = node[part]
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError:
            parsed = value
        node[parts[-1]] = parsed
        applied[dotted] = parsed
    return cfg, applied


def out_dir(cfg: dict, deck: Path) -> Path:
    q = Path(str(cfg["global"]["output_dir"]))
    return q if q.is_absolute() else (deck.parent / q).resolve()


def cmd_run(deck: Path, overrides: list[str]) -> Path:
    from pinnacle.plugins import discover_components
    from pinnacle.runner import Runner

    discover_components(force=True)
    cfg, applied = load_deck(deck, overrides)
    target = out_dir(cfg, deck)
    target.mkdir(parents=True, exist_ok=True)
    cfg["global"]["output_dir"] = str(target)
    here = Path.cwd()
    os.chdir(deck.parent)
    try:
        runner = Runner(cfg)
        runner.build()
        tic = time.time()
        result = runner.model.train()["composed"]
        wall = time.time() - tic
    finally:
        os.chdir(here)
    mesh = runner.model.composed.mesh
    nodes = np.asarray(mesh.nodes, dtype=float)
    length = float(nodes[:, LOAD_COMPONENT].max() - nodes[:, LOAD_COMPONENT].min())
    pull = max(
        (float(bc["values"][0]) for bc in cfg["bcs"].values()
         if bc.get("components") == [LOAD_COMPONENT] and float(bc["values"][0]) != 0.0),
        key=abs,
    )
    history = result["history"]["cp_quasi_static"]
    t_end = float(cfg["executioner"]["t_end"])
    key = f"mean_stress_{LOAD_COMPONENT}{LOAD_COMPONENT}"
    rows = [(0.0, 0.0, 0.0, 0.0)]
    for row in history:
        strain = pull / length * float(row["time"]) / t_end
        cauchy = float(row[key])
        rows.append((float(row["time"]), strain, cauchy, cauchy / (1.0 + strain)))
    np.savetxt(target / "curve.csv", np.array(rows), delimiter=",", comments="",
               header="time_s,eng_strain,mean_cauchy_stress_Pa,eng_stress_Pa")
    n_e = int(mesh.n_elements)
    stress = np.asarray(result["stress"], dtype=float).reshape(n_e, -1, 3, 3).mean(axis=1)
    eqps = np.asarray(result["equivalent_plastic_strain"], dtype=float).reshape(n_e, -1).mean(axis=1)
    arrays = dict(nodes=nodes, connectivity=np.asarray(mesh.connectivity, dtype=np.int32),
                  displacement=np.asarray(result["displacement"], dtype=float).reshape(-1, 3),
                  stress_element=stress, eqps_element=eqps)
    frames = result.get("snapshots", []) or []
    if frames:
        arrays["snapshot_displacement"] = np.stack(
            [np.asarray(f["displacement"], dtype=np.float32).reshape(-1, 3) for f in frames if "displacement" in f]
        ) if "displacement" in frames[0] else np.zeros((0, 0, 3), dtype=np.float32)
    np.savez_compressed(target / "fields.npz", **arrays)
    summary = {
        "wall_s": round(wall, 1), "steps": len(history), "elements": n_e,
        "newton_iterations": int(sum(r["iterations"] for r in history)),
        "cutbacks": int(history[-1]["cutbacks"]) if history else 0,
        "unconverged_steps": int(sum(not r["converged"] for r in history)),
        "overrides": applied,
    }
    (target / "wall.json").write_text(json.dumps(summary, indent=1))
    print(target.name, json.dumps(summary))
    return target


def read_target(path: Path):
    data = np.genfromtxt(path, delimiter=",", skip_header=2)
    return data[:, 0], data[:, 1] * 1e6, data[:, 2] * 1e6


def residual(run_dir: Path, strain_t: np.ndarray, stress_t: np.ndarray):
    curve = np.genfromtxt(run_dir / "curve.csv", delimiter=",", skip_header=1)
    use = (strain_t >= 0.01) & (strain_t <= curve[-1, 1])
    model = np.interp(strain_t[use], curve[:, 1], curve[:, 3])
    return model - stress_t[use], strain_t[use]


def deck_values(text: str) -> dict[str, float]:
    return {k: float(re.search(rf"^\s*{k}\s*=\s*([-+0-9.eE]+)", text, re.M).group(1)) for k in FITTED}


def write_deck(template: str, values: dict[str, float], output: str, path: Path) -> None:
    text = template
    for key, value in values.items():
        text = re.sub(rf"^(\s*{key}\s*=\s*)[-+0-9.eE]+", lambda m, v=value: f"{m.group(1)}{v:.6e}", text,
                      count=1, flags=re.M)
    text = re.sub(r"^(\s*output_dir\s*=\s*).*$", lambda m: f"{m.group(1)}{output}", text, count=1, flags=re.M)
    path.write_text(text)


def cmd_calibrate(start: Path, target_csv: Path, iterations: int, tag: str) -> None:
    strain_t, stress_t, _sd = read_target(target_csv)
    template = start.read_text()
    values = deck_values(template)
    cfg, _ = load_deck(start, [])
    root_rel = Path(str(cfg["global"]["output_dir"])).parent / f"cal{tag}"
    log = []
    previous = None
    for it in range(iterations):
        runs = {}
        jobs = []
        for name in ("center",) + FITTED:
            trial = dict(values)
            if name != "center":
                trial[name] *= 1.05
            deck = start.parent / f"cal{tag}_{it}_{name}.i"
            write_deck(template, trial, f"{root_rel}/it{it}_{name}", deck)
            runs[name] = out_dir(load_deck(deck, [])[0], deck)
            runs[name].mkdir(parents=True, exist_ok=True)
            # the four runs of an iteration are independent: one process each
            env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("K3_THREADS", "3"))
            jobs.append(subprocess.Popen(
                [sys.executable, "-W", "ignore", str(Path(__file__).resolve()), "run", str(deck)],
                env=env, stdout=(runs[name] / "run.log").open("w"), stderr=subprocess.STDOUT))
        codes = [job.wait() for job in jobs]
        if any(codes):
            raise SystemExit(f"iteration {it}: a run failed (exit codes {codes}); see run.log")
        r0, used = residual(runs["center"], strain_t, stress_t)
        rms = float(np.sqrt(np.mean(r0**2)))
        columns = []
        for name in FITTED:
            r1, _ = residual(runs[name], strain_t, stress_t)
            columns.append((r1 - r0) / np.log(1.05))
        jac = np.array(columns).T
        normal = jac.T @ jac
        step = -np.linalg.solve(normal + 1.0e-3 * np.diag(np.diag(normal)), jac.T @ r0)
        step = np.clip(step, -np.log(2.0), np.log(2.0))
        predicted = float(np.sqrt(np.mean((r0 + jac @ step) ** 2)))
        entry = {
            "iteration": it, "values": dict(values), "rms_MPa": rms / 1e6,
            "max_abs_MPa": float(np.abs(r0).max() / 1e6), "points": int(used.size),
            "log_step": {k: float(s) for k, s in zip(FITTED, step)},
            "predicted_rms_MPa": predicted / 1e6,
            "jacobian_condition": float(np.linalg.cond(jac)),
        }
        log.append(entry)
        print(json.dumps(entry), flush=True)
        (out_dir(cfg, start).parent / f"cal{tag}").mkdir(parents=True, exist_ok=True)
        (out_dir(cfg, start).parent / f"cal{tag}" / "calibration.json").write_text(json.dumps(log, indent=1))
        if previous is not None and abs(previous - rms) < 0.02 * rms:
            break
        previous = rms
        values = {k: values[k] * float(np.exp(s)) for k, s in zip(FITTED, step)}
    best = min(log, key=lambda e: e["rms_MPa"])
    print("best iteration", best["iteration"], "rms %.2f MPa" % best["rms_MPa"], best["values"])


def cmd_figure(target_csv: Path, out: Path, pairs: list[str]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    strain_t, stress_t, sd = read_target(target_csv)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.fill_between(100 * strain_t, (stress_t - sd) / 1e6, (stress_t + sd) / 1e6, color="0.8",
                    label="measured, mean $\\pm$ s.d. of 8 specimens")
    ax.plot(100 * strain_t, stress_t / 1e6, color="k", linewidth=1.2)
    for pair in pairs:
        label, run_dir = pair.split("=", 1)
        curve = np.genfromtxt(Path(run_dir) / "curve.csv", delimiter=",", skip_header=1)
        ax.plot(100 * curve[:, 1], curve[:, 3] / 1e6, label=label, linewidth=1.2)
    ax.set_xlabel("engineering strain (%)")
    ax.set_ylabel("engineering stress (MPa)")
    ax.set_xlim(0, 22)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out)


def cmd_cube(run_dir: Path, out: Path, corner: list[float], edge: float) -> None:
    from pinnacle.materials.rotations import matrix_to_bunge

    data = np.load(run_dir / "fields.npz")
    nodes, eta, rotations = data["nodes"], data["eta"], data["rotations"]
    low = np.array(corner)
    inside = np.all((nodes >= low - 1e-12) & (nodes <= low + edge + 1e-12), axis=1)
    solid = (eta.astype(float) ** 2).sum(axis=1) >= 0.25
    keep = inside & solid
    label = eta[keep].argmax(axis=1)
    angles = matrix_to_bunge(rotations[label]) % 360.0
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(f"# grain structure of {run_dir.name}: cube corner {corner}, edge {edge}; "
                 f"{int(keep.sum())} of {int(inside.sum())} nodes solid; "
                 f"{np.unique(label).size} orientations\n")
        fh.write("x,y,z,phi1,Phi,phi2,grain_id\n")
        for point, angle, grain in zip(nodes[keep] - low, angles, label):
            fh.write("%.9e,%.9e,%.9e,%.6f,%.6f,%.6f,%d\n" % (*point, *angle, grain))
    print(f"wrote {out}: {int(keep.sum())} points, {np.unique(label).size} orientations, "
          f"solid fraction of the cube {keep.sum() / max(inside.sum(), 1):.3f}")


def cmd_ebsd(volume: Path, out: Path, corner_um: list[float], edge_um: float, stride: int) -> None:
    """Cut a cube out of the NIST 3-D EBSD volume (build coordinates) as an orientation point file.

    `corner_um`: lower corner in the build frame, x and y measured from the corner of the
    volume, z from the plate top (z = 0 of the record). Every `stride`-th voxel is kept.
    """
    import h5py

    with h5py.File(volume, "r") as h:
        container = h["DataContainers/EBSDDataContainer"]
        geometry = container["_SIMPL_GEOMETRY"]
        origin = np.asarray(geometry["ORIGIN"][()], dtype=float)
        spacing = np.asarray(geometry["SPACING"][()], dtype=float)
        low = np.array([corner_um[0], corner_um[1], corner_um[2] - origin[2]])
        first = np.ceil(low / spacing - 0.5).astype(int)
        count = np.floor(edge_um / spacing).astype(int)
        sl = tuple(slice(int(first[d]), int(first[d] + count[d]), stride) for d in (2, 1, 0))
        euler = np.asarray(container["CellData/EulerAngles"][sl], dtype=float)
        sample = np.asarray(container["CellData/SampleVoxels"][sl])[..., 0] > 0
    if not sample.all():
        raise SystemExit(f"the cube is not inside the sample: {sample.mean():.3f} of its voxels are")
    index = [np.arange(first[d], first[d] + count[d], stride) for d in (0, 1, 2)]
    centers = [((index[d] + 0.5) * spacing[d] - low[d]) * 1e-6 for d in (0, 1, 2)]
    zz, yy, xx = np.meshgrid(centers[2], centers[1], centers[0], indexing="ij")
    angles = np.degrees(euler.reshape(-1, 3)) % 360.0
    phi2, big = euler[..., 2].ravel(), euler[..., 1].ravel()
    # crystal direction along the build direction (Bunge, sample -> crystal: third column)
    direction = np.abs(np.stack([np.sin(phi2) * np.sin(big), np.cos(phi2) * np.sin(big), np.cos(big)], axis=1))
    direction.sort(axis=1)
    poles = {"100": np.array([0, 0, 1.0]), "110": np.array([0, 1, 1.0]) / np.sqrt(2), "111": np.ones(3) / np.sqrt(3)}
    texture = {k: float((direction @ v >= np.cos(np.radians(15.0))).mean()) for k, v in poles.items()}
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(f"# 3-D EBSD of NIST mds2-2767 ({volume.name}), build frame; cube corner x, y = "
                 f"{corner_um[0]}, {corner_um[1]} um from the corner of the volume, z = {corner_um[2]} um above "
                 f"the plate top; edge {edge_um} um; every {stride}. voxel ({spacing[0] * stride:g} x "
                 f"{spacing[1] * stride:g} x {spacing[2] * stride:g} um); Bunge angles in degrees\n")
        fh.write("x,y,z,phi1,Phi,phi2,grain_id\n")
        for i, (x, y, z, a) in enumerate(zip(xx.ravel(), yy.ravel(), zz.ravel(), angles)):
            fh.write("%.6e,%.6e,%.6e,%.4f,%.4f,%.4f,%d\n" % (x, y, z, *a, i))
    print(f"wrote {out}: {angles.shape[0]} voxels; build direction within 15 deg of "
          f"<100>, <110>, <111>: {texture['100']:.3f}, {texture['110']:.3f}, {texture['111']:.3f}")


def main() -> None:
    args = sys.argv[1:]
    if len(args) >= 2 and args[0] == "run":
        cmd_run(Path(args[1]).resolve(), args[2:])
    elif len(args) >= 3 and args[0] == "calibrate":
        cmd_calibrate(Path(args[1]).resolve(), Path(args[2]).resolve(),
                      int(args[3]) if len(args) > 3 else 6, args[4] if len(args) > 4 else "")
    elif len(args) == 7 and args[0] == "cube":
        cmd_cube(Path(args[1]), Path(args[2]), [float(v) for v in args[3:6]], float(args[6]))
    elif len(args) in (7, 8) and args[0] == "ebsd":
        cmd_ebsd(Path(args[1]), Path(args[2]), [float(v) for v in args[3:6]], float(args[6]),
                 int(args[7]) if len(args) == 8 else 1)
    elif len(args) >= 3 and args[0] == "figure":
        cmd_figure(Path(args[1]), Path(args[2]), args[3:])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
