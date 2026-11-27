#!/usr/bin/env python3
"""Part V, cases K2 and after: run a build-column deck (temperature + grains) and keep its fields.

usage: scripts/lpbf_column.py <deck.i> run [key=value ...] | grains | thermal | figure | dwell

  run     run the PINNACLE deck in this process and write, into the deck's output_dir:
          fields.npz (nodes, final temperature and order parameters, the orientation of
          every order parameter as rotation matrices, snapshot times and grain maps),
          history.csv (per step: time, pool measures, nuclei), wall.json.
          `key=value` overrides a deck key for a check run, e.g. utilities.scan.layers=1
          executioner.end_time=5.6e-3 global.output_dir=... (the override is recorded).
  grains  segment the final order parameters into grains (connected regions of one order
          parameter) and write grains.csv and grain_summary.json: equivalent diameter,
          extent along the build direction and across it, and the fractions of the built
          volume whose build direction lies within 15 degrees of <100>, <110>, <111>
  thermal time above melting and cooling rate at freezing of the last layer's surface
          elements, beside the thermography medians of the NIST record
  figure  sections of the grain map
  dwell   from history.csv: how many elements are still liquid when each track and each
          layer starts (the two dwell times of the scan path are long enough when all
          entries are 0), and when the last liquid froze

Several GPUs: with global.rank_devices="cuda:0 cuda:2 ..." (in the deck or as an override
of `run`) the mesh is cut into one part per device and every part runs in its own
process (PINNACLE ADR 0057). The files are the same. Every stored frame is reduced on
its rank to the grain map, the solid share, the temperature and the velocity, so the
order parameters of the frames never reach the host memory of one process.

PINNACLE must be importable (PYTHONPATH). Figures need matplotlib (not on the host).
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent


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


def frame_of_column(frame: dict) -> dict:
    """Snapshot filter: what the animations need of a frame, in small types."""
    eta = np.asarray(frame["eta"])
    out = {
        "step": frame["step"],
        "time": frame["time"],
        "grain": eta.argmax(axis=1).astype(np.int16),
        "solid": (eta**2).sum(axis=1).astype(np.float16),
        "temperature": np.asarray(frame["temperature"], dtype=np.float32).ravel(),
    }
    if "velocity" in frame:
        out["velocity"] = np.asarray(frame["velocity"], dtype=np.float16)
    return out


def cmd_run(deck: Path, overrides: list[str]) -> None:
    import os

    from pinnacle.plugins import discover_components
    from pinnacle.runner import Runner

    discover_components(force=True)
    cfg, applied = load(deck, overrides)
    target = out_dir(cfg, deck)
    target.mkdir(parents=True, exist_ok=True)
    cfg["global"]["output_dir"] = str(target)
    for name, util in cfg.get("utilities", {}).items():
        if "record_file" in util and not applied.get(f"utilities.{name}.record_file"):
            util["record_file"] = str(target / "process_state.npz")
        for key, name_ in (("series_file", "temperature_series.npz"), ("state_series_file", "consolidated_series.npz")):
            if key in util and not applied.get(f"utilities.{name}.{key}"):
                util[key] = str(target / name_)
    os.chdir(deck.parent)
    ranks = None
    if int(cfg["global"].get("ranks", 1)) > 1 or cfg["global"].get("rank_devices"):
        from pinnacle.mesh.builder import use_mesh
        from pinnacle.mesh.generators import box_mesh
        from pinnacle.parallel import run_decomposed

        tic = time.time()
        out = run_decomposed(cfg, snapshot_filter=frame_of_column)
        ranks = out["ranks"]
        wall = max(r["wall_s"] for r in ranks)
        total = time.time() - tic
        result = {"time": out["time"], "snapshots": out["snapshots"], **out["fields"]}
        history = out["history"]
        nodes = np.asarray(out["nodes"], dtype=float)
        n_elements = int(out["n_elements"])
        # the orientation of each order parameter is a property of the deck: read it from
        # the kernel of the same deck built on a mesh of one element
        small = dict(cfg, **{"global": {k: v for k, v in cfg["global"].items() if k not in ("ranks", "rank_devices", "partition")}})
        small["global"]["use_gpu"] = False
        with use_mesh(box_mesh(bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)), divisions=(1, 1, 1))):
            probe = Runner(small)
            probe.build()
        kernel = next(k for k in probe.model.composed.kernels if hasattr(k, "rotations"))
    else:
        runner = Runner(cfg)
        runner.build()
        runner.model.composed.executioner.snapshot_filter = frame_of_column
        tic = time.time()
        result = runner.model.train()["composed"]
        wall = total = time.time() - tic
        problem = runner.model.composed
        nodes = np.asarray(problem.mesh.nodes, dtype=float)
        n_elements = int(problem.mesh.n_elements)
        history = result["history"]["pf_fem_transient"]
        kernel = next(k for k in problem.kernels if hasattr(k, "rotations"))
    eta = np.asarray(result["eta"], dtype=np.float32)
    frames = result.get("snapshots", [])
    np.savez_compressed(
        target / "fields.npz",
        nodes=nodes,
        temperature=np.asarray(result["temperature"], dtype=np.float32).ravel(),
        eta=eta,
        rotations=kernel.rotations(eta.shape[1]),
        snapshot_times=np.array([f["time"] for f in frames]),
        snapshot_grain=np.stack([np.asarray(f["grain"], dtype=np.int16) for f in frames])
        if frames else np.zeros((0, 0), dtype=np.int16),
        snapshot_solid=np.stack([np.asarray(f["solid"], dtype=np.float16) for f in frames])
        if frames else np.zeros((0, 0), dtype=np.float16),
        snapshot_temperature=np.stack([np.asarray(f["temperature"], dtype=np.float32).ravel() for f in frames])
        if frames else np.zeros((0, 0), dtype=np.float32),
            **(
            {
                "velocity": np.asarray(result["velocity"], dtype=np.float32),
                "snapshot_velocity": np.stack(
                    [np.asarray(f["velocity"], dtype=np.float16) for f in frames]
                ),
            }
            if frames and "velocity" in frames[0]
            else {}
        ),
    )
    with (target / "history.csv").open("w", newline="") as fh:
        keys = ["step", "t", "dt", "temperature_max", "eta_max"]
        service_keys = sorted({f"{s}.{k}" for row in history[:1] for s, d in row.get("services", {}).items() for k in d})
        w = csv.writer(fh)
        w.writerow(keys + service_keys)
        for row in history:
            services = row.get("services", {})
            w.writerow([row.get(k) for k in keys] + [services.get(s.split(".")[0], {}).get(s.split(".", 1)[1]) for s in service_keys])
    summary = {
        "wall_s": round(wall, 1),
        "steps": len(history),
        "s_per_step": round(wall / max(len(history), 1), 4),
        "elements": n_elements,
        "nodes": int(nodes.shape[0]),
        "order_parameters": int(eta.shape[1]),
        "time_end": float(result["time"]),
        "overrides": applied,
        "run_end": result.get("services", {}),
        "nucleated_total": history[-1].get("services", {}).get("nuclei", {}).get("nucleated_total"),
        "backend": str(cfg["global"].get("backend")),
        "use_gpu": bool(cfg["global"].get("use_gpu")),
    }
    if ranks is not None:
        summary["wall_with_start_and_merge_s"] = round(total, 1)
        summary["ranks"] = ranks
    (target / "wall.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


def grid_of(nodes: np.ndarray):
    axes = [np.unique(np.round(nodes[:, d], 12)) for d in range(3)]
    order = np.lexsort((nodes[:, 0], nodes[:, 1], nodes[:, 2]))
    shape = (axes[2].size, axes[1].size, axes[0].size)
    return axes, order, shape


def cmd_grains(deck: Path, overrides: list[str]) -> None:
    from scipy import ndimage

    cfg, _ = load(deck, overrides)
    target = out_dir(cfg, deck)
    data = np.load(target / "fields.npz")
    nodes, eta, rot = data["nodes"], data["eta"], data["rotations"]
    axes, order, shape = grid_of(nodes)
    h = float(axes[0][1] - axes[0][0])
    component = eta.argmax(axis=1)[order].reshape(shape)
    solid = ((eta.astype(float) ** 2).sum(axis=1)[order].reshape(shape)) > 0.5
    plate_top = float(cfg["utilities"]["state"].get("substrate_top", 0.0))
    z = axes[2]
    built = solid & (z[:, None, None] > plate_top + 0.5 * h)
    rows = []
    label_total = 0
    labels = np.zeros(shape, dtype=np.int32)
    for comp in np.unique(component[solid]):
        lab, count = ndimage.label((component == comp) & solid)
        lab[lab > 0] += label_total
        labels += lab
        label_total += count
    index = np.arange(1, label_total + 1)
    volume = ndimage.sum(np.ones(shape), labels, index) * h**3
    slices = ndimage.find_objects(labels)
    for gid, sl in zip(index, slices):
        if sl is None:
            continue
        zs = z[sl[0]]
        if zs[-1] <= plate_top + 0.5 * h:
            continue  # entirely in the plate
        comp = int(component[sl][labels[sl] == gid][0])
        rows.append({
            "grain": int(gid), "order_parameter": comp, "volume_um3": volume[gid - 1] * 1e18,
            "d_eq_um": 1e6 * (6.0 * volume[gid - 1] / np.pi) ** (1.0 / 3.0),
            "height_um": 1e6 * (zs[-1] - zs[0] + h),
            "width_um": 1e6 * h * max(sl[1].stop - sl[1].start, sl[2].stop - sl[2].start),
            "z_min_um": 1e6 * zs[0], "touches_plate": bool(zs[0] <= plate_top + 0.5 * h),
        })
    with (target / "grains.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()})
    # texture of the built volume: build direction (z) in the crystal frame of each voxel
    bd = rot[:, 2, :]  # row 2 of R (crystal to sample): sample z expressed in the crystal frame
    a = np.sort(np.abs(bd), axis=1)[:, ::-1]
    cos100 = a[:, 0]
    cos110 = (a[:, 0] + a[:, 1]) / np.sqrt(2.0)
    cos111 = a.sum(axis=1) / np.sqrt(3.0)
    limit = np.cos(np.deg2rad(15.0))
    comp_built = component[built]
    fractions = {name: float(np.mean((c >= limit)[comp_built])) for name, c in (("100", cos100), ("110", cos110), ("111", cos111))}
    # a grain needs a few elements to be told from a fragment of a diffuse boundary
    min_cells = 8
    kept = [r for r in rows if r["volume_um3"] >= min_cells * (1e6 * h) ** 3 * (1 - 1e-9)]
    fragments = len(rows) - len(kept)
    rows = kept
    d = np.array([r["d_eq_um"] for r in rows])
    new = np.array([not r["touches_plate"] for r in rows])
    summary = {
        "mesh_um": 1e6 * h, "min_grain_cells": min_cells, "fragments_below_min": fragments,
        "grains_in_build": len(rows),
        "grains_not_touching_plate": int(new.sum()),
        "d_eq_um_median": float(np.median(d)), "d_eq_um_mean": float(d.mean()),
        "height_um_median": float(np.median([r["height_um"] for r in rows])),
        "aspect_median": float(np.median([r["height_um"] / r["width_um"] for r in rows])),
        "solid_fraction_above_plate": float(
            built.sum() / max(int((z > plate_top + 0.5 * h).sum()) * shape[1] * shape[2], 1)
        ),
        "texture_fraction_within_15deg": fractions,
        "measured": {"d_eq_um_median": 11.5, "texture_fraction_within_15deg": {"100": 0.10, "110": 0.21, "111": 0.10},
                     "source": "3-D EBSD of schwalbach 2024 sectioningdata as analysed in ambench-data/agent (not re-derived)"},
    }
    (target / "grain_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


def cmd_thermal(deck: Path, overrides: list[str]) -> None:
    """Time above melting and cooling rate at freezing of the top surface of the last layer."""
    cfg, _ = load(deck, overrides)
    target = out_dir(cfg, deck)
    d = np.load(target / "process_state.npz")
    c = d["centroid"]
    z = np.unique(np.round(c[:, 2], 12))
    built = d["melt_count"] > 0
    top = z[np.nonzero([np.any(built & (np.abs(c[:, 2] - zz) < 1e-12)) for zz in z])[0][-1]]
    surface = built & (np.abs(c[:, 2] - top) < 1e-12)
    tam = d["time_above_melting"][surface]
    rate = d["freeze_cooling_rate"][surface]
    rate = rate[np.isfinite(rate)]
    grad = d["freeze_gradient"][built]
    grad = grad[np.isfinite(grad)]
    cool = d["freeze_cooling_rate"][built]
    cool = cool[np.isfinite(cool)]
    summary = {
        "surface_height_um": float(top * 1e6), "surface_elements": int(surface.sum()),
        "time_above_melting_s": {k: float(np.percentile(tam, q)) for k, q in (("p25", 25), ("p50", 50), ("p75", 75))},
        "cooling_rate_at_freezing_K_s": {k: float(np.percentile(rate, q)) for k, q in (("p25", 25), ("p50", 50), ("p75", 75))},
        "all_molten_elements": {
            "gradient_K_m_p50": float(np.median(grad)), "cooling_rate_K_s_p50": float(np.median(cool)),
            "front_speed_m_s_p50": float(np.median(cool / np.maximum(grad[: cool.size], 1e-30))) if cool.size == grad.size else None,
        },
        "measured_leg9": {"time_above_melting_s_p50": 1.276e-3, "solid_cooling_rate_K_s_p50": 1.894e5,
                          "source": "data/lpbf/in718/k2_thermography_leg9_summary.csv"},
    }
    (target / "thermal_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


def cmd_figure(deck: Path, overrides: list[str]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfg, _ = load(deck, overrides)
    target = out_dir(cfg, deck)
    data = np.load(target / "fields.npz")
    nodes, eta, rot = data["nodes"], data["eta"], data["rotations"]
    axes, order, shape = grid_of(nodes)
    component = eta.argmax(axis=1)[order].reshape(shape)
    solid = (eta.astype(float) ** 2).sum(axis=1)[order].reshape(shape)
    # inverse-pole-figure color of the build direction per order parameter
    bd = np.abs(rot[:, 2, :]); bd = np.sort(bd, axis=1)
    rgb = np.stack([bd[:, 2] - bd[:, 1], bd[:, 1] - bd[:, 0], bd[:, 0] * np.sqrt(3.0)], axis=1)
    rgb /= rgb.max(axis=1, keepdims=True)
    image = rgb[component]
    image[solid < 0.5] = 1.0
    um = [1e6 * a for a in axes]
    plate_top = 1e6 * float(cfg["utilities"]["state"].get("substrate_top", 0.0))
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3))
    jy = shape[1] // 2
    ax[0].imshow(image[:, jy, :], origin="lower", extent=[um[0][0], um[0][-1], um[2][0], um[2][-1]], interpolation="nearest")
    ax[0].axhline(plate_top, color="k", lw=0.6, ls="--")
    ax[0].set_xlabel("x ($\\mu$m)"); ax[0].set_ylabel("z, build direction ($\\mu$m)"); ax[0].set_title("section y = %.0f $\\mu$m" % um[1][jy])
    ix = shape[2] // 2
    ax[1].imshow(image[:, :, ix], origin="lower", extent=[um[1][0], um[1][-1], um[2][0], um[2][-1]], interpolation="nearest")
    ax[1].axhline(plate_top, color="k", lw=0.6, ls="--")
    ax[1].set_xlabel("y ($\\mu$m)"); ax[1].set_ylabel("z ($\\mu$m)"); ax[1].set_title("section x = %.0f $\\mu$m" % um[0][ix])
    filled = np.nonzero((solid > 0.5).mean(axis=(1, 2)) > 0.2)[0]
    top_built = um[2][filled[-1]] if filled.size else um[2][-1]
    kz = int(np.argmin(np.abs(um[2] - (plate_top + 0.5 * max(top_built - plate_top, 0.0)))))
    ax[2].imshow(image[kz], origin="lower", extent=[um[0][0], um[0][-1], um[1][0], um[1][-1]], interpolation="nearest")
    ax[2].set_xlabel("x ($\\mu$m)"); ax[2].set_ylabel("y ($\\mu$m)"); ax[2].set_title("layer plane z = %.0f $\\mu$m" % um[2][kz])
    fig.suptitle("Grains colored by the crystal direction along the build direction (inverse pole figure colors)", fontsize=10)
    fig.tight_layout()
    out = REPO / "figures" / "lpbf" / f"k2_grains_{target.name}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out.relative_to(REPO))


def cmd_dwell(deck: Path, overrides: list[str]) -> None:
    """Liquid elements at the start of every track and layer; the time of the last liquid."""
    from pinnacle.components.additive.scan_path import ScanPath

    cfg, _ = load(deck, overrides)
    target = out_dir(cfg, deck)
    scan_cfg = next(u for u in cfg["utilities"].values() if u.get("kind") == "scan_path")
    scan = ScanPath({k: v for k, v in scan_cfg.items() if k not in ("kind", "type")})
    with (target / "history.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["t"]) for r in rows])
    liquid = np.array([float(r["state.liquid_elements"]) for r in rows])
    # a run with `hold_scan_until_solid` moved its way-points: replay the holds
    held = np.array([float(r.get("state.scan_delay") or 0.0) for r in rows])
    for k in np.nonzero(np.diff(np.concatenate([[0.0], held])) > 0.0)[0]:
        scan.delay(t[k] * (1.0 - 1e-12), held[k] - (held[k - 1] if k else 0.0))
    report = {"tracks": [], "recoat": scan.recoat, "scan_delay_ms": round(1e3 * float(held[-1]), 4)}
    previous_z = None
    for segment in scan.segments():
        if segment["on"] <= 0.0 or segment["t0"] > t[-1]:
            continue
        # the state at the end of the last step before the beam comes on
        index = int(np.searchsorted(t, segment["t0"] * (1.0 + 1e-12), side="right") - 1)
        left = int(liquid[index]) if index >= 0 else 0
        report["tracks"].append({
            "start_ms": round(1e3 * segment["t0"], 4),
            "first_of_layer": bool(previous_z is None or segment["z"] > previous_z + 1e-12),
            "liquid_elements_at_start": left,
        })
        previous_z = segment["z"]
    wet = np.nonzero(liquid > 0)[0]
    report["last_liquid_ms"] = round(1e3 * float(t[wet[-1]]), 4) if wet.size else None
    report["end_ms"] = round(1e3 * float(t[-1]), 4)
    report["liquid_at_end"] = int(liquid[-1])
    report["tracks_started_on_liquid"] = int(sum(1 for r in report["tracks"] if r["liquid_elements_at_start"] > 0))
    (target / "dwell_check.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    deck = Path(sys.argv[1]).resolve()
    command, rest = sys.argv[2], sys.argv[3:]
    {"run": cmd_run, "grains": cmd_grains, "thermal": cmd_thermal, "figure": cmd_figure, "dwell": cmd_dwell}.get(command, lambda *_: (_ for _ in ()).throw(SystemExit(__doc__)))(deck, rest)


if __name__ == "__main__":
    main()
