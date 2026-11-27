#!/usr/bin/env python3
"""Part V, defect indicators: keyhole numbers of the single tracks, the depth law, lack of fusion.

usage: scripts/lpbf_defects.py law <gan supplementary data 1 .xlsx> <out.json>
       scripts/lpbf_defects.py tracks <k1 targets.csv> <out.csv> <absorptivity> [<absorptivity> ...]
       scripts/lpbf_defects.py figure <gan .xlsx> <tracks.csv> <out.png>
       scripts/lpbf_defects.py compare <out.png> <label>=<k1_cross_sections.csv> ...
       scripts/lpbf_defects.py column <column run dir> <x0> <y0> <x1> <y1> <plate top> [hatch layer width depth]
       scripts/lpbf_defects.py cracking <run dir> <out.png> <out.csv> <name> [<name> ...]
       scripts/lpbf_defects.py pockets <run dir> <plate surface height> [<out.json>]
       scripts/lpbf_defects.py pores <run dir> <out.png> <surface height>

  law      recompute the keyhole number of every row of the published keyhole-depth data
           (gan 2021 keyhole, Supplementary Data 1: power, measured absorptivity, speed,
           spot radius, properties, keyhole depth) with PINNACLE's `keyhole_number` and
           compare the depth law e / r0 = 0.4 (Ke - 1.4) with the measured depths: the
           check that the criterion is coded as published
  tracks   for the NIST AMB2022-03 single tracks: keyhole number and regime, normalized
           enthalpy, depth of the vapor depression by the law, for each stated
           absorptivity, beside the measured depth, width and width over depth
  figure   keyhole depth over spot radius against the keyhole number: published data,
           the law, and the measured melt depths of the NIST tracks
  compare  depth and width of the single tracks, computed against measured, for several
           studies of scripts/lpbf_k1.py (fitted depth, keyhole-number depth)
  column   lack of fusion in a build-column run: deposited elements of the scanned patch
           (x0..x1, y0..y1, above the plate top) that never melted, as connected voids;
           with hatch, layer, pool width and depth also the overlap index

  cracking for a run with `cracking_indicators` services that wrote cracking_<name>.json
           and cracking_<name>.npz: one table row per solid-fraction curve (end of the
           mushy zone, steepness index, time ratio, pressure drops) and a figure with
           the curves, the maps of the pressure drop and of the time ratio in the
           section across the tracks, and the distribution of the time ratio

  pockets  for a free-surface run (fields.npz with the phase field `phi`, and
           process_state.npz): the gas enclosed below the surface at the end of the run
           (PINNACLE `gas_pockets`: connected gas elements without a path through the
           gas to the top), what is left of the depression, and the largest rise of the
           surface above the plate

  pores    for a run with a `pore_tracers` service (pores.npz, pores.json,
           process_state.npz): where the trapped pores sit in the melted material (two
           sections), and how long the lost ones lived

Material values of alloy 718 are those of the decks (docs/16_lpbf_notebook.md, section
2): volumetric capacity and liquid conductivity at the liquidus. PINNACLE must be
importable (PYTHONPATH); `law` and `figure` need openpyxl, `figure` matplotlib.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

LIQUIDUS, PLATE = 1610.15, 296.65
CAPACITY, K_LIQUID = 6.3473e6, 26.80      # J/(m3 K), W/(m K) at the liquidus
BOILING = 3186.0


def gan_rows(path: Path) -> list[dict]:
    import openpyxl

    sheet = openpyxl.load_workbook(path, data_only=True)["Sheet1"]
    rows, label = [], None
    for r in list(sheet.iter_rows(values_only=True))[1:]:
        if r[0]:
            label = str(r[0]).replace("μ", "u")
        if r[1] is None:
            continue
        p, eta, v, d, r0, k, rho, cp, alpha, _, _, tl_t0, _, e = (float(x) for x in r[1:15])
        rows.append({"set": label, "power": p, "absorptivity": eta, "speed": v, "radius": r0,
                     "capacity": rho * cp, "diffusivity": alpha, "rise": tl_t0, "depth": e})
    return rows


def with_ke(rows: list[dict]) -> np.ndarray:
    from pinnacle.components.additive.defects import keyhole_number

    return np.array([
        keyhole_number(r["absorptivity"], r["power"], r["speed"], r["radius"], r["rise"], 0.0,
                       r["capacity"], r["diffusivity"]) for r in rows])


def cmd_law(xlsx: Path, out: Path) -> None:
    rows = gan_rows(xlsx)
    ke = with_ke(rows)
    aspect = np.array([r["depth"] / r["radius"] for r in rows])
    law = np.maximum(0.4 * (ke - 1.4), 0.0)
    slope, offset = np.polyfit(ke, aspect, 1)
    report = {
        "rows": len(rows),
        "sets": {s: sum(1 for r in rows if r["set"] == s) for s in sorted({r["set"] for r in rows})},
        "keyhole_number_range": [float(ke.min()), float(ke.max())],
        "absorptivity_range": [min(r["absorptivity"] for r in rows), max(r["absorptivity"] for r in rows)],
        "law": "e / r0 = 0.4 (Ke - 1.4)",
        "r2_of_the_law": float(1.0 - ((aspect - law) ** 2).sum() / ((aspect - aspect.mean()) ** 2).sum()),
        "rms_of_the_law_in_e_over_r0": float(np.sqrt(((aspect - law) ** 2).mean())),
        "least_squares_line": {"slope": float(slope), "onset": float(-offset / slope)},
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


def track_rows(targets: Path, absorptivities: list[float]) -> list[dict]:
    from pinnacle.components.additive.defects import keyhole_depth, keyhole_number, normalized_enthalpy

    alpha = K_LIQUID / CAPACITY
    out = []
    with targets.open() as fh:
        reader = csv.DictReader(line for line in fh if not line.startswith(('"#', "#")))
        for r in reader:
            power, speed = float(r["power_W"]), 1e-3 * float(r["speed_mm_s"])
            radius = 0.5e-6 * float(r["d4sigma_um"])
            row = {"line": r["line"], "power_W": power, "speed_mm_s": 1e3 * speed,
                   "radius_um": 1e6 * radius, "depth_meas_um": float(r["depth_um"]),
                   "width_meas_um": float(r["width_um"]),
                   "width_over_depth_meas": float(r["width_um"]) / float(r["depth_um"]),
                   "depth_over_radius_meas": float(r["depth_um"]) / (1e6 * radius)}
            for eta in absorptivities:
                ke = keyhole_number(eta, power, speed, radius, LIQUIDUS, PLATE, CAPACITY, alpha)
                tag = f"{eta:g}"
                row[f"Ke_{tag}"] = ke
                row[f"regime_{tag}"] = "conduction" if ke < 1.4 else "transition" if ke < 6.0 else "keyhole"
                row[f"keyhole_depth_um_{tag}"] = 1e6 * keyhole_depth(ke, radius)
                row[f"normalized_enthalpy_{tag}"] = normalized_enthalpy(
                    eta, power, speed, radius, LIQUIDUS, CAPACITY, alpha)
            out.append(row)
    return out


def cmd_tracks(targets: Path, out: Path, absorptivities: list[float]) -> None:
    rows = track_rows(targets, absorptivities)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"normalized-enthalpy threshold pi T_b / T_m = {np.pi * BOILING / LIQUIDUS:.2f}")
    for r in rows:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", out)


def cmd_figure(xlsx: Path, tracks: Path, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = gan_rows(xlsx)
    ke = with_ke(rows)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    markers = {"Ti64": "o", "Ti6": "s", "Al6061": "^", "SS316": "D"}
    for label in sorted({r["set"] for r in rows}):
        pick = np.array([r["set"] == label for r in rows])
        ax.scatter(ke[pick], [r["depth"] / r["radius"] for r, p in zip(rows, pick) if p], s=16,
                   marker=markers.get(label.split(",")[0], "o"), facecolors="none", edgecolors="0.45",
                   label=f"keyhole depth, {label} (gan 2021)")
    line = np.linspace(0.0, 34.0, 50)
    ax.plot(line, np.maximum(0.4 * (line - 1.4), 0.0), color="k", lw=1.0, label="e / r$_0$ = 0.4 (Ke $-$ 1.4)")
    with tracks.open() as fh:
        data = list(csv.DictReader(fh))
    etas = sorted({k[3:] for k in data[0] if k.startswith("Ke_")}, key=float)
    for eta, color in zip(etas, ("tab:blue", "tab:red", "tab:green")):
        ax.scatter([float(r[f"Ke_{eta}"]) for r in data], [float(r["depth_over_radius_meas"]) for r in data],
                   s=42, color=color, marker="*", label=f"melt depth, NIST alloy 718 tracks, absorptivity {eta}")
    ax.axvline(6.0, color="0.7", lw=0.7, ls="--")
    ax.set_xlabel("keyhole number Ke"); ax.set_ylabel("depth / beam radius")
    ax.set_xlim(0.0, 34.0); ax.set_ylim(0.0, 13.0)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out)


def cmd_column(run_dir: Path, box: list[float], plate_top: float, rule: list[float]) -> None:
    from scipy import ndimage

    data = np.load(run_dir / "process_state.npz")
    centroid, volume, state = data["centroid"], data["volume"], data["state"]
    axes = [np.unique(np.round(centroid[:, d], 12)) for d in range(3)]
    index = tuple(np.searchsorted(axes[d], np.round(centroid[:, d], 12)) for d in range(3))
    grid = np.zeros(tuple(a.size for a in axes), dtype=bool)
    patch = ((centroid[:, 0] >= box[0]) & (centroid[:, 0] <= box[2]) & (centroid[:, 1] >= box[1])
             & (centroid[:, 1] <= box[3]) & (centroid[:, 2] > plate_top))
    built = patch & (state > 0.5)
    unmelted = built & (state < 1.5)
    grid[index[0][unmelted], index[1][unmelted], index[2][unmelted]] = True
    label, count = ndimage.label(grid)          # face neighbours
    element = float(np.median(volume))
    sizes = ndimage.sum(grid, label, index=np.arange(1, count + 1)) if count else np.zeros(0)
    diameters = (6.0 * sizes * element / np.pi) ** (1.0 / 3.0)
    layers = {}
    for z in np.unique(np.round(centroid[built, 2], 12)):
        row = built & np.isclose(centroid[:, 2], z)
        layers[f"{1e6 * z:.0f}"] = float(volume[row & unmelted].sum() / volume[row].sum())
    report = {
        "run": run_dir.name,
        "patch": box, "plate_top": plate_top,
        "built_elements": int(built.sum()), "unmelted_elements": int(unmelted.sum()),
        "unmelted_fraction": float(volume[unmelted].sum() / max(volume[built].sum(), 1e-300)),
        "voids": int(count),
        "void_equivalent_diameter_um": {
            "max": float(1e6 * diameters.max()) if count else 0.0,
            "median": float(1e6 * np.median(diameters)) if count else 0.0},
        "largest_void_elements": int(sizes.max()) if count else 0,
        "unmelted_fraction_by_element_layer_um": layers,
    }
    if len(rule) == 4:
        from pinnacle.components.additive.defects import overlap_index

        report["overlap_index"] = overlap_index(*[rule[i] for i in (0, 2, 1, 3)])
        report["overlap_rule_inputs"] = {"hatch": rule[0], "layer": rule[1], "width": rule[2], "depth": rule[3]}
    (run_dir / "lack_of_fusion.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


def cmd_compare(out: Path, pairs: list[str]) -> None:
    """Depth and width, computed against measured, for several K1 studies (label=csv)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.1))
    markers = ("o", "s", "^", "D")
    colors = ("tab:blue", "tab:red", "tab:green", "tab:purple")
    for (label, path), marker, color in zip((p.split("=", 1) for p in pairs), markers, colors):
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        for ax, key in zip(axes, ("depth", "width")):
            meas = np.array([float(r[f"{key}_meas_um"]) for r in rows])
            sim = np.array([float(r[f"{key}_um"]) for r in rows])
            rms = float(np.sqrt((((sim / meas) - 1.0) ** 2).mean()) * 100.0)
            ax.scatter(meas, sim, s=30, marker=marker, color=color, label=f"{label} (rms {rms:.1f} %)")
    for ax, key in zip(axes, ("depth", "width")):
        lim = [80.0, 280.0] if key == "depth" else [80.0, 180.0]
        ax.plot(lim, lim, color="0.3", lw=0.8)
        ax.plot(lim, [0.9 * v for v in lim], color="0.6", lw=0.6, ls="--")
        ax.plot(lim, [1.1 * v for v in lim], color="0.6", lw=0.6, ls="--")
        ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
        ax.set_xlabel(f"measured {key} ($\\mu$m)"); ax.set_ylabel(f"computed {key} ($\\mu$m)")
        ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out)


def cmd_cracking(run_dir: Path, out: Path, table: Path, names: list[str]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    reports = {n: json.loads((run_dir / f"cracking_{n}.json").read_text()) for n in names}
    maps = {n: np.load(run_dir / f"cracking_{n}.npz") for n in names}
    rows = []
    for n in names:
        r, sf = reports[n], reports[n]["solid_fraction"]
        rows.append({
            "curve": n, "elements": r["elements"],
            "end_temperature_K": round(sf["end_temperature"], 1), "end_fraction": round(sf["end_fraction"], 3),
            "steepness_index_K": round(r["steepness_index"], 0),
            "time_ratio_uniform_cooling": round(r["time_ratio_at_constant_cooling_rate"], 3),
            "time_ratio_median": round(r["time_ratio"]["median"], 3),
            "time_ratio_p05": round(r["time_ratio"]["p05"], 3), "time_ratio_p95": round(r["time_ratio"]["p95"], 3),
            "time_vulnerable_median_us": round(1e6 * r["time_vulnerable"]["median"], 2),
            "time_relief_median_us": round(1e6 * r["time_relief"]["median"], 2),
            "dp_shrinkage_median_MPa": round(1e-6 * r["pressure_drop_shrinkage"]["median"], 3),
            "dp_deformation_median_MPa": round(1e-6 * r["pressure_drop_deformation"]["median"], 4),
            "dp_total_p05_MPa": round(1e-6 * r["pressure_drop"]["p05"], 3),
            "dp_total_p95_MPa": round(1e-6 * r["pressure_drop"]["p95"], 2),
            "strain_rate_median_per_s": round(r["strain_rate"]["median"], 2),
            "max_strain_rate_median_per_s": round(r["max_strain_rate"]["median"], 2),
            "share_above_cavitation": round(r["share_above_cavitation"], 3),
        })
    table.parent.mkdir(parents=True, exist_ok=True)
    with table.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    for row in rows:
        print(json.dumps(row))

    fig, axes = plt.subplots(len(names), 3, figsize=(12.5, 3.0 * len(names)), squeeze=False)
    centroid = maps[names[0]]["centroid"]
    size = float(np.min(np.diff(np.unique(np.round(centroid[:, 0], 9)))))
    plane = np.unique(np.round(centroid[:, 1], 9))
    y0 = plane[np.argmin(np.abs(plane - 0.5 * (plane[0] + plane[-1])))]
    pick = np.abs(centroid[:, 1] - y0) < 0.25 * size
    xs, zs = np.unique(np.round(centroid[pick, 0], 9)), np.unique(np.round(centroid[pick, 2], 9))
    extent = [1e6 * (xs[0] - size / 2), 1e6 * (xs[-1] + size / 2), 1e6 * (zs[0] - size / 2), 1e6 * (zs[-1] + size / 2)]

    def image(values: np.ndarray) -> np.ndarray:
        grid = np.full((zs.size, xs.size), np.nan)
        ix = np.searchsorted(xs, np.round(centroid[pick, 0], 9))
        iz = np.searchsorted(zs, np.round(centroid[pick, 2], 9))
        grid[iz, ix] = values[pick]
        return grid

    drops = [maps[n]["pressure_drop_deformation"] + maps[n]["pressure_drop_shrinkage"] for n in names]
    finite = np.concatenate([d[np.isfinite(d)] for d in drops]) * 1e-6
    norm = LogNorm(vmin=float(np.percentile(finite, 1)), vmax=float(np.percentile(finite, 99)))
    ratios = [maps[n]["time_ratio"] for n in names]
    top = float(np.nanpercentile(np.concatenate(ratios), 97))
    colors = ("tab:blue", "tab:red", "tab:green")
    for i, n in enumerate(names):
        ax = axes[i, 0]
        im = ax.imshow(image(drops[i] * 1e-6), origin="lower", extent=extent, norm=norm, cmap="viridis", aspect="equal")
        fig.colorbar(im, ax=ax, shrink=0.8, label="pressure drop (MPa)")
        ax.set_title(f"{n}: pressure drop", fontsize=9)
        ax = axes[i, 1]
        im = ax.imshow(image(ratios[i]), origin="lower", extent=extent, vmin=0.0, vmax=top, cmap="magma", aspect="equal")
        fig.colorbar(im, ax=ax, shrink=0.8, label="$t_V / t_R$")
        ax.set_title(f"{n}: time ratio", fontsize=9)
        for ax in axes[i, :2]:
            ax.set_xlabel("x across the tracks (um)"); ax.set_ylabel("z (um)")
    ax = axes[0, 2]
    for n, color in zip(names, colors):
        sf = reports[n]["solid_fraction"]
        values = (sf["temperature"], sf["fraction"])
        ax.plot(values[0], values[1], color=color, label=f"{n} (ends at {sf['end_temperature']:.0f} K, $f_s$ {sf['end_fraction']:.2f})")
    ax.set_xlabel("temperature (K)"); ax.set_ylabel("solid fraction $f_s$"); ax.legend(fontsize=7)
    ax.set_title("solid-fraction curves", fontsize=9)
    if len(names) > 1:
        ax = axes[1, 2]
        bins = np.linspace(0.0, top, 60)
        for n, r, color in zip(names, ratios, colors):
            ax.hist(r[np.isfinite(r)], bins=bins, histtype="step", color=color, label=n)
            ax.axvline(reports[n]["time_ratio_at_constant_cooling_rate"], color=color, lw=0.7, ls="--")
        ax.set_xlabel("$t_V / t_R$ (dashed: uniform cooling rate)"); ax.set_ylabel("elements"); ax.legend(fontsize=7)
        ax.set_title("distribution of the time ratio", fontsize=9)
        for extra in axes[2:, 2]:
            extra.axis("off")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out, "section at y =", y0)


def cmd_pockets(run_dir: Path, surface: float, out: Path | None) -> dict:
    from scipy.spatial import cKDTree

    from pinnacle.components.additive.pores import gas_pockets

    fields = np.load(run_dir / "fields.npz")
    record = np.load(run_dir / "process_state.npz")
    nodes, phi = fields["nodes"], np.asarray(fields["phi"], dtype=float).ravel()
    temperature = np.asarray(fields["temperature"], dtype=float).ravel()
    centroid, volume = record["centroid"], record["volume"]
    size = float(np.min(np.diff(np.unique(np.round(centroid[:, 2], 12)))))
    # the phase of an element: the mean over its eight nodes (the nearest to its centre)
    _, corner = cKDTree(nodes).query(centroid, k=2 ** centroid.shape[1])
    element_phi = phi[corner].mean(axis=1)
    element_t = temperature[corner].mean(axis=1)
    is_gas = element_phi < 0.5
    top = float(nodes[:, -1].max())
    report = gas_pockets(centroid, volume, is_gas, size, top, element_t < LIQUIDUS)
    below = is_gas & (centroid[:, -1] < surface)
    above = ~is_gas & (centroid[:, -1] > surface)
    report["run"] = run_dir.name
    report["element_um"] = 1e6 * size
    report["remaining_depression_um"] = 1e6 * float(surface - (centroid[below, -1].min() - size / 2)) if below.any() else 0.0
    report["largest_rise_um"] = 1e6 * float(centroid[above, -1].max() + size / 2 - surface) if above.any() else 0.0
    report["metal_volume_change_pct"] = 100.0 * float((volume * element_phi).sum() / volume[centroid[:, -1] < surface].sum() - 1.0)
    # during the run: enclosed gas in every stored frame
    if "snapshot_phi" in fields.files:
        frames = fields["snapshot_phi"]
        frames = frames.reshape(frames.shape[0], -1)
        times = fields["snapshot_times"]
        seen = []
        for k in range(frames.shape[0]):
            gas_k = frames[k][corner].mean(axis=1) < 0.5
            found = gas_pockets(centroid, volume, gas_k, size, top)
            if found["pockets"]:
                seen.append({"time_us": 1e6 * float(times[k]), "pockets": found["pockets"],
                             "diameter_max_um": 1e6 * found["pocket_diameter_max"],
                             "depth_um": [1e6 * (surface - e["center"][-1]) for e in found["pocket_list"]]})
        report["frames"] = int(frames.shape[0])
        report["frames_with_pockets"] = len(seen)
        report["transient_pockets"] = seen
        report["transient_diameter_max_um"] = max((f["diameter_max_um"] for f in seen), default=0.0)
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k not in ("pocket_list", "transient_pockets")}))
    return report


def cmd_pores(run_dir: Path, out: Path, surface: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    table = np.load(run_dir / "pores.npz")
    report = json.loads((run_dir / "pores.json").read_text())
    record = np.load(run_dir / "process_state.npz")
    position, status, diameter = table["position"], table["status"], table["diameter"]
    trapped, lost = status == 1, status == 2
    melted = record["melt_count"] > 0
    centroid = record["centroid"]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    for ax, axis, label in ((axes[0], 0, "x"), (axes[1], 1, "y")):
        ax.scatter(1e6 * centroid[melted, axis], 1e6 * centroid[melted, 2], s=2, color="0.88", marker="s",
                   label="melted material")
        ax.scatter(1e6 * position[trapped, axis], 1e6 * position[trapped, 2], s=1.0e11 * diameter[trapped] ** 2,
                   color="tab:red", edgecolors="k", linewidths=0.4, label=f"trapped pores ({int(trapped.sum())})")
        ax.axhline(1e6 * surface, color="k", lw=0.6, ls="--")
        ax.set_xlabel(f"{label} (um)"); ax.set_ylabel("z (um)")
        ax.set_ylim(1e6 * surface - 180.0, 1e6 * surface + 10.0)
    axes[0].legend(fontsize=7, loc="lower left")
    life = 1e6 * (table["settled_at"] - table["released_at"])[lost]
    axes[2].hist(life[np.isfinite(life)], bins=np.linspace(0.0, max(400.0, float(np.nanmax(life))), 40), color="tab:blue")
    axes[2].set_xlabel("time from release to loss at the surface (us)"); axes[2].set_ylabel("pores")
    axes[2].set_title(f"{int(lost.sum())} of {report['released']} pores lost at the surface", fontsize=9)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    depth = 1e6 * (surface - position[trapped, 2])
    print(json.dumps({"released": report["released"], "lost": int(lost.sum()), "trapped": int(trapped.sum()),
                      "trapped_depth_um_median": float(np.median(depth)) if trapped.any() else None,
                      "trapped_depth_um_max": float(depth.max()) if trapped.any() else None,
                      "life_lost_us_median": float(np.nanmedian(life)) if lost.any() else None,
                      "gas_porosity": report.get("gas_porosity")}))
    print("wrote", out)


def main() -> None:
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == "law":
        cmd_law(Path(a[1]), Path(a[2]))
    elif len(a) >= 4 and a[0] == "tracks":
        cmd_tracks(Path(a[1]), Path(a[2]), [float(v) for v in a[3:]])
    elif len(a) == 4 and a[0] == "figure":
        cmd_figure(Path(a[1]), Path(a[2]), Path(a[3]))
    elif len(a) in (7, 11) and a[0] == "column":
        cmd_column(Path(a[1]), [float(v) for v in a[2:6]], float(a[6]), [float(v) for v in a[7:]])
    elif len(a) == 4 and a[0] == "pores":
        cmd_pores(Path(a[1]), Path(a[2]), float(a[3]))
    elif len(a) in (3, 4) and a[0] == "pockets":
        cmd_pockets(Path(a[1]), float(a[2]), Path(a[3]) if len(a) == 4 else None)
    elif len(a) >= 5 and a[0] == "cracking":
        cmd_cracking(Path(a[1]), Path(a[2]), Path(a[3]), a[4:])
    elif len(a) >= 3 and a[0] == "compare":
        cmd_compare(Path(a[1]), a[2:])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
