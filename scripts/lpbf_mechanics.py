#!/usr/bin/env python3
"""Part V: run a process-stress deck (crystal plasticity under a temperature history).

usage: scripts/lpbf_mechanics.py <deck.i> run [key=value ...] | summary
       scripts/lpbf_mechanics.py cooldown <series.npz> <out.npz> <temperature K> <time s>

  run      run the PINNACLE deck in this process and write, into its output_dir:
           fields.npz (nodes, connectivity, displacement, the Cauchy stress and the
           equivalent plastic strain at the quadrature points of the last step, element
           means of both, and the snapshots the executioner kept), history.csv (per load
           step: time, Newton iterations, cutbacks, mean stress, temperature range, molten
           points), wall.json. `key=value` overrides a deck key (recorded in wall.json).
  summary  print the final mean and extreme stresses from fields.npz
  cooldown copy a temperature series and append one frame of uniform temperature at the
           given time: the cool-down of the block to the plate temperature, which the
           thermal run does not step through

PINNACLE must be importable (PYTHONPATH).
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


def load(deck: Path, overrides: list[str]):
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


def von_mises(stress: np.ndarray) -> np.ndarray:
    dev = stress - np.trace(stress, axis1=-2, axis2=-1)[..., None, None] * np.eye(3) / 3.0
    return np.sqrt(1.5 * np.sum(dev * dev, axis=(-2, -1)))


def cmd_run(deck: Path, overrides: list[str]) -> None:
    from pinnacle.plugins import discover_components
    from pinnacle.runner import Runner

    discover_components(force=True)
    cfg, applied = load(deck, overrides)
    target = out_dir(cfg, deck)
    target.mkdir(parents=True, exist_ok=True)
    cfg["global"]["output_dir"] = str(target)
    os.chdir(deck.parent)
    runner = Runner(cfg)
    runner.build()
    tic = time.time()
    result = runner.model.train()["composed"]
    wall = time.time() - tic
    mesh = runner.model.composed.mesh
    history = result["history"]["cp_quasi_static"]
    stress = np.asarray(result["stress"], dtype=float)
    n_e = int(mesh.n_elements)
    stress = stress.reshape(n_e, -1, 3, 3)
    eqps = np.asarray(result["equivalent_plastic_strain"], dtype=float).reshape(n_e, -1)
    frames = result.get("snapshots", []) or []
    np.savez_compressed(
        target / "fields.npz",
        nodes=np.asarray(mesh.nodes, dtype=float),
        connectivity=np.asarray(mesh.connectivity, dtype=np.int32),
        displacement=np.asarray(result["displacement"], dtype=float).reshape(-1, 3),
        stress_qp=stress.astype(np.float32),
        stress_element=stress.mean(axis=1),
        von_mises_element=von_mises(stress.mean(axis=1)),
        eqps_element=eqps.mean(axis=1),
        snapshot_keys=np.array(sorted({k for f in frames for k in f}), dtype=object),
    )
    keys = ["step", "time", "dt", "iterations", "converged", "cutbacks", "mean_stress_00",
            "mean_stress_11", "mean_stress_22", "mean_stress_01", "mean_stress_02", "mean_stress_12",
            "heat_temperature_max", "heat_temperature_min", "heat_molten_points", "newton_time_s"]
    with (target / "history.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for row in history:
            w.writerow([row.get(k) for k in keys])
    summary = {
        "wall_s": round(wall, 1), "steps": len(history),
        "s_per_step": round(wall / max(len(history), 1), 2),
        "elements": n_e, "quadrature_points": int(stress.shape[0] * stress.shape[1]),
        "newton_iterations": int(sum(r["iterations"] for r in history)),
        "cutbacks": int(history[-1]["cutbacks"]) if history else 0,
        "unconverged_steps": int(sum(not r["converged"] for r in history)),
        "time_end": float(history[-1]["time"]) if history else 0.0,
        "overrides": applied, "backend": str(cfg["global"].get("backend")),
    }
    (target / "wall.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    cmd_summary(deck, overrides)


def cmd_summary(deck: Path, overrides: list[str]) -> None:
    cfg, _ = load(deck, overrides)
    data = np.load(out_dir(cfg, deck) / "fields.npz", allow_pickle=True)
    stress = data["stress_element"]
    vm = data["von_mises_element"]
    report = {
        "mean_stress_MPa": {k: float(stress[:, i, j].mean() / 1e6) for k, (i, j) in
                            {"xx": (0, 0), "yy": (1, 1), "zz": (2, 2), "xy": (0, 1)}.items()},
        "von_mises_MPa": {"mean": float(vm.mean() / 1e6), "p95": float(np.percentile(vm, 95) / 1e6),
                          "max": float(vm.max() / 1e6)},
        "sigma_xx_MPa": {"min": float(stress[:, 0, 0].min() / 1e6), "max": float(stress[:, 0, 0].max() / 1e6)},
        "eqps": {"mean": float(data["eqps_element"].mean()), "max": float(data["eqps_element"].max())},
    }
    print(json.dumps(report, indent=1))
    (out_dir(cfg, deck) / "stress_summary.json").write_text(json.dumps(report, indent=1))


def cmd_cooldown(series: Path, out: Path, temperature: float, time: float) -> None:
    """Append one frame of uniform temperature to a stored series (the cool-down to the plate)."""
    data = np.load(series)
    times, values = data["times"], data["values"]
    if time <= times[-1]:
        raise SystemExit(f"the cool-down time {time} is not after the last frame {times[-1]}")
    extra = {k: data[k] for k in data.files if k not in ("times", "values")}
    np.savez_compressed(
        out,
        times=np.append(times, time),
        values=np.concatenate([values, np.full((1, values.shape[1]), temperature, dtype=values.dtype)]),
        **extra,
    )
    print(f"wrote {out}: {times.size + 1} frames, last frame {temperature} K at {time} s "
          f"(last computed frame: {times[-1]:.6g} s, {values[-1].min():.1f} to {values[-1].max():.1f} K)")


def main() -> None:
    if len(sys.argv) == 6 and sys.argv[1] == "cooldown":
        cmd_cooldown(Path(sys.argv[2]), Path(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
        return
    if len(sys.argv) < 3 or sys.argv[2] not in ("run", "summary"):
        raise SystemExit(__doc__)
    deck = Path(sys.argv[1]).resolve()
    (cmd_run if sys.argv[2] == "run" else cmd_summary)(deck, sys.argv[3:])


if __name__ == "__main__":
    main()
