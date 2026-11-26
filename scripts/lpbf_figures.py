#!/usr/bin/env python3
"""Part V: figures of the report that are not drawn by the case scripts.

usage: scripts/lpbf_figures.py track <process_state.npz> <depth_um> <width_um> <out.png>
       scripts/lpbf_figures.py stress <mechanics output dir> <thermal process_state.npz> <out.png>
       scripts/lpbf_figures.py column_stress <mechanics output dir> <column process_state.npz> <out.png>

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
    ]
    for a, field, extent, cmap, lo, hi, title, mask, (u, v), labels in panels:
        im = a.imshow(field, origin="lower", extent=extent, cmap=cmap, vmin=lo, vmax=hi)
        if mask.min() < 0.5 < mask.max():
            a.contour(1e6 * u, 1e6 * v, mask, levels=[0.5], colors="k", linewidths=0.8)
        a.set_title(title, fontsize=9)
        a.set_xlabel(labels[0]); a.set_ylabel(labels[1])
        fig.colorbar(im, ax=a, shrink=0.85)
    # mean stress of every element layer against the height (the paper's Figure 15d)
    _, y_grid = grid(centroid, data["stress_element"][:, 1, 1] / 1e6)
    ax[1, 1].plot(1e6 * z, s_grid.mean(axis=(1, 2)), label="$\\sigma_{xx}$")
    ax[1, 1].plot(1e6 * z, y_grid.mean(axis=(1, 2)), label="$\\sigma_{yy}$")
    ax[1, 1].axhline(0.0, color="0.6", linewidth=0.6)
    ax[1, 1].set_xlabel("z ($\\mu$m)"); ax[1, 1].set_ylabel("mean stress of the element layer (MPa)")
    ax[1, 1].legend(frameon=False, fontsize=8)
    print("layer means z_um, sxx, syy:")
    for zz, a_, b_ in zip(1e6 * z, s_grid.mean(axis=(1, 2)), y_grid.mean(axis=(1, 2))):
        print(f"  {zz:6.1f} {a_:9.1f} {b_:9.1f}")
    rows = list(csv.DictReader((base / "history.csv").open()))
    t = np.array([float(r["time"]) for r in rows]) * 1e3
    for key, label in (("mean_stress_00", "$\\sigma_{xx}$"), ("mean_stress_11", "$\\sigma_{yy}$"), ("mean_stress_22", "$\\sigma_{zz}$")):
        ax[0, 2].plot(t, [float(r[key]) / 1e6 for r in rows], label=label)
    ax[0, 2].set_xlabel("time (ms)"); ax[0, 2].set_ylabel("mean stress of the volume (MPa)")
    ax[0, 2].legend(frameon=False, fontsize=8)
    ax[1, 2].plot(t, [float(r["heat_temperature_max"]) for r in rows], color="C3", label="maximum")
    ax[1, 2].plot(t, [float(r["heat_temperature_min"]) for r in rows], color="C0", label="minimum")
    ax[1, 2].set_xlabel("time (ms)"); ax[1, 2].set_ylabel("temperature of the quadrature points (K)")
    ax[1, 2].set_yscale("log"); ax[1, 2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print("wrote", out)


def cmd_column_stress(directory: str, record: str, out: str) -> None:
    """Residual stress of a build column (case K4): sections, layer means, history."""
    import csv

    data = np.load(Path(directory) / "fields.npz", allow_pickle=True)
    nodes, conn = data["nodes"], data["connectivity"]
    stress, mises = data["stress_element"] / 1e6, data["von_mises_element"] / 1e6
    center = nodes[conn].mean(axis=1)
    axes = [np.unique(np.round(center[:, d], 12)) for d in range(3)]
    size = [float(np.diff(a).mean()) for a in axes]
    index = [np.searchsorted(axes[d], np.round(center[:, d], 12)) for d in range(3)]
    shape = tuple(a.size for a in axes)
    # consolidated share of every mechanics element from the record of the thermal run
    rec = np.load(record)
    cell = [np.clip(np.floor(rec["centroid"][:, d] / size[d]).astype(int), 0, shape[d] - 1) for d in range(3)]
    flat = (cell[2] * shape[1] + cell[1]) * shape[0] + cell[0]
    count = np.bincount(flat, minlength=int(np.prod(shape)))
    share = np.bincount(flat, weights=(rec["state"] > 1.5).astype(float), minlength=count.size) / np.maximum(count, 1)
    solid = share[(index[2] * shape[1] + index[1]) * shape[0] + index[0]] > 0.5

    def field(values: np.ndarray) -> np.ndarray:
        grid_ = np.full(shape, np.nan)
        grid_[index[0][solid], index[1][solid], index[2][solid]] = values[solid]
        return grid_

    sxx, vm = field(stress[:, 0, 0]), field(mises)
    um = [a * 1e6 for a in axes]
    half = [0.5e6 * v for v in size]
    lim = 800.0
    fig, ax = plt.subplots(2, 2, figsize=(11.0, 7.6), layout="constrained")
    jy, kz = shape[1] // 2, min(shape[2] - 1, int(np.searchsorted(axes[2], 250e-6)))
    ext_xz = [um[0][0] - half[0], um[0][-1] + half[0], um[2][0] - half[2], um[2][-1] + half[2]]
    ext_xy = [um[0][0] - half[0], um[0][-1] + half[0], um[1][0] - half[1], um[1][-1] + half[1]]
    im = ax[0, 0].imshow(sxx[:, jy, :].T, origin="lower", extent=ext_xz, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal")
    ax[0, 0].set_title(f"$\\sigma_{{xx}}$ (MPa), section y = {um[1][jy]:.0f} $\\mu$m (gray: not consolidated)", fontsize=9)
    ax[0, 0].set_xlabel("x ($\\mu$m)"); ax[0, 0].set_ylabel("z, build direction ($\\mu$m)")
    ax[0, 0].set_facecolor("0.8")
    fig.colorbar(im, ax=ax[0, 0], shrink=0.8)
    im = ax[0, 1].imshow(sxx[:, :, kz].T, origin="lower", extent=ext_xy, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal")
    ax[0, 1].set_title(f"$\\sigma_{{xx}}$ (MPa), element layer z = {um[2][kz]:.0f} $\\mu$m", fontsize=9)
    ax[0, 1].set_xlabel("x ($\\mu$m)"); ax[0, 1].set_ylabel("y ($\\mu$m)")
    ax[0, 1].set_facecolor("0.8")
    fig.colorbar(im, ax=ax[0, 1], shrink=0.8)
    print("layer means of the consolidated elements: z_um, sxx, syy, szz, von Mises")
    rows = []
    for k in range(shape[2]):
        pick = solid & (index[2] == k)
        if pick.any():
            rows.append((um[2][k], stress[pick, 0, 0].mean(), stress[pick, 1, 1].mean(),
                         stress[pick, 2, 2].mean(), mises[pick].mean()))
            print("%6.0f %8.1f %8.1f %8.1f %8.1f" % rows[-1])
    rows_a = np.array(rows)
    ax[1, 0].plot(rows_a[:, 1], rows_a[:, 0], "o-", label="$\\sigma_{xx}$")
    ax[1, 0].plot(rows_a[:, 2], rows_a[:, 0], "s-", label="$\\sigma_{yy}$")
    ax[1, 0].plot(rows_a[:, 3], rows_a[:, 0], "^-", color="0.5", label="$\\sigma_{zz}$")
    ax[1, 0].axvline(0.0, color="k", lw=0.6); ax[1, 0].axhline(200.0, color="k", lw=0.6, ls="--")
    ax[1, 0].set_xlabel("mean stress of the element layer (MPa)"); ax[1, 0].set_ylabel("z ($\\mu$m); dashed: top of the plate")
    ax[1, 0].legend(fontsize=8)
    with (Path(directory) / "solver_history.csv").open() as fh:
        history = list(csv.DictReader(fh))
    t = np.array([float(r["time"]) for r in history]) * 1e3
    for key, label in (("mean_stress_00", "$\\sigma_{xx}$"), ("mean_stress_11", "$\\sigma_{yy}$")):
        ax[1, 1].plot(t, np.array([float(r[key]) for r in history]) / 1e6, label=label)
    ax[1, 1].set_xlabel("time (ms)"); ax[1, 1].set_ylabel("mean stress of the block (MPa)")
    ax[1, 1].axhline(0.0, color="k", lw=0.6); ax[1, 1].legend(fontsize=8)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170)
    print("wrote", out, f"({int(solid.sum())} of {solid.size} elements consolidated; von Mises mean "
          f"{mises[solid].mean():.0f}, max {mises[solid].max():.0f} MPa)")


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    if sys.argv[1] == "track" and len(sys.argv) == 6:
        cmd_track(*sys.argv[2:6])
    elif sys.argv[1] == "stress" and len(sys.argv) == 5:
        cmd_stress(*sys.argv[2:5])
    elif sys.argv[1] == "column_stress" and len(sys.argv) == 5:
        cmd_column_stress(*sys.argv[2:5])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
