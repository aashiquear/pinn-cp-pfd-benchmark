#!/usr/bin/env python3
"""Part V: figures of the report that are not drawn by the case scripts.

usage: scripts/lpbf_figures.py track <process_state.npz> <depth_um> <width_um> <out.png>
       scripts/lpbf_figures.py stress <mechanics output dir> <thermal process_state.npz> <out.png>

  track   cross section of a single track: peak temperature reached by every element in
          the plane across the track at mid-length, the liquidus and solidus isotherms,
          and the measured depth and width as a box
  stress  process stress of a track: sigma_xx (along the track) and the von Mises stress
          of the element means on the section across the track at mid-length and on the
          top face, the outline of the molten zone, and the history of the mean stress
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def grid(centroid: np.ndarray, values: np.ndarray):
    axes = [np.unique(np.round(centroid[:, d], 12)) for d in range(3)]
    order = np.lexsort((centroid[:, 0], centroid[:, 1], centroid[:, 2]))
    return axes, values[order].reshape(axes[2].size, axes[1].size, axes[0].size)


def edges(axis: np.ndarray) -> list[float]:
    h = axis[1] - axis[0]
    return [1e6 * (axis[0] - h / 2), 1e6 * (axis[-1] + h / 2)]


def cmd_track(record: str, depth: str, width: str, out: str) -> None:
    d = np.load(record)
    axes, peak = grid(d["centroid"], d["peak_temperature"])
    x, y, z = axes
    section = peak[:, :, x.size // 2]
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    im = ax.imshow(np.clip(section, None, 3500.0), origin="lower", extent=edges(y) + edges(z),
                   cmap="inferno", vmin=300.0, vmax=3500.0)
    ax.contour(1e6 * y, 1e6 * z, section, levels=[1528.15, 1610.15], colors=["w", "c"], linewidths=1.0)
    top = edges(z)[1]
    yc = 1e6 * y.mean()
    w, dep = float(width), float(depth)
    ax.plot([yc - w / 2, yc - w / 2, yc + w / 2, yc + w / 2], [top, top - dep, top - dep, top],
            color="lime", lw=1.2, ls="--")
    ax.set_xlabel("across the track ($\\mu$m)")
    ax.set_ylabel("height ($\\mu$m)")
    ax.set_ylim(top - 2.2 * dep, top)
    fig.colorbar(im, ax=ax, label="peak temperature (K), cut at 3500")
    ax.plot([], [], color="w", label="solidus isotherm")
    ax.plot([], [], color="c", label="liquidus isotherm")
    ax.plot([], [], color="lime", ls="--", label="measured width and depth")
    ax.legend(loc="lower right", fontsize=7, facecolor="0.4", labelcolor="w")
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print("wrote", out)


def cmd_stress(directory: str, record: str, out: str) -> None:
    base = Path(directory)
    data = np.load(base / "fields.npz", allow_pickle=True)
    nodes, conn = data["nodes"], data["connectivity"]
    centroid = nodes[conn].mean(axis=1)
    sxx = data["stress_element"][:, 0, 0] / 1e6
    vm = data["von_mises_element"] / 1e6
    thermal = np.load(record)
    axes, s_grid = grid(centroid, sxx)
    _, v_grid = grid(centroid, vm)
    _, melt = grid(thermal["centroid"], (thermal["melt_count"] > 0).astype(float))
    x, y, z = axes
    ix, kz = x.size // 2, z.size - 1
    fig, ax = plt.subplots(2, 3, figsize=(13.0, 6.6))
    lim = float(np.percentile(np.abs(sxx), 99))
    panels = [
        (ax[0, 0], s_grid[:, :, ix], edges(y) + edges(z), "RdBu_r", -lim, lim, "$\\sigma_{xx}$ (MPa), section across the track", melt[:, :, ix], (y, z), ("y ($\\mu$m)", "z ($\\mu$m)")),
        (ax[0, 1], s_grid[kz], edges(x) + edges(y), "RdBu_r", -lim, lim, "$\\sigma_{xx}$ (MPa), top layer of elements", melt[kz], (x, y), ("x, scan direction ($\\mu$m)", "y ($\\mu$m)")),
        (ax[1, 0], v_grid[:, :, ix], edges(y) + edges(z), "viridis", 0.0, float(np.percentile(vm, 99)), "von Mises stress (MPa), section across the track", melt[:, :, ix], (y, z), ("y ($\\mu$m)", "z ($\\mu$m)")),
        (ax[1, 1], v_grid[kz], edges(x) + edges(y), "viridis", 0.0, float(np.percentile(vm, 99)), "von Mises stress (MPa), top layer of elements", melt[kz], (x, y), ("x, scan direction ($\\mu$m)", "y ($\\mu$m)")),
    ]
    for a, field, extent, cmap, lo, hi, title, mask, (u, v), labels in panels:
        im = a.imshow(field, origin="lower", extent=extent, cmap=cmap, vmin=lo, vmax=hi)
        if mask.min() < 0.5 < mask.max():
            a.contour(1e6 * u, 1e6 * v, mask, levels=[0.5], colors="k", linewidths=0.8)
        a.set_title(title, fontsize=9)
        a.set_xlabel(labels[0]); a.set_ylabel(labels[1])
        fig.colorbar(im, ax=a, shrink=0.85)
    rows = list(csv.DictReader((base / "history.csv").open()))
    t = np.array([float(r["time"]) for r in rows]) * 1e3
    for key, label in (("mean_stress_00", "$\\sigma_{xx}$"), ("mean_stress_11", "$\\sigma_{yy}$"), ("mean_stress_22", "$\\sigma_{zz}$")):
        ax[0, 2].plot(t, [float(r[key]) / 1e6 for r in rows], label=label)
    ax[0, 2].set_xlabel("time (ms)"); ax[0, 2].set_ylabel("mean stress of the volume (MPa)")
    ax[0, 2].legend(frameon=False, fontsize=8); ax[0, 2].set_xscale("log")
    ax[1, 2].plot(t, [float(r["heat_temperature_max"]) for r in rows], color="C3", label="maximum")
    ax[1, 2].plot(t, [float(r["heat_temperature_min"]) for r in rows], color="C0", label="minimum")
    ax[1, 2].set_xlabel("time (ms)"); ax[1, 2].set_ylabel("temperature of the quadrature points (K)")
    ax[1, 2].set_xscale("log"); ax[1, 2].set_yscale("log"); ax[1, 2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out)


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    if sys.argv[1] == "track" and len(sys.argv) == 6:
        cmd_track(*sys.argv[2:6])
    elif sys.argv[1] == "stress" and len(sys.argv) == 5:
        cmd_stress(*sys.argv[2:5])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
