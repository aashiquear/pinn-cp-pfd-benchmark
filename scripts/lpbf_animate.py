#!/usr/bin/env python3
"""Part V: animations from the field series of a process run (scripts/lpbf_run.py).

usage: scripts/lpbf_animate.py track <run_dir> <out.gif> [liquidus_K] [solidus_K]
       scripts/lpbf_animate.py grains <column_run_dir> <out.gif>

  track   one frame per stored snapshot of fields.npz: the temperature in the plane along
          the track (through the track axis) and in the top layer of nodes, with the
          liquidus and solidus isotherms; when the run has a velocity field, arrows of
          the in-plane velocity and the largest speed. The frames are also written as
          PNG files next to the GIF (<out>_frames/), the stack the figures of the report
          are taken from.
  grains  the same for a build-column run of scripts/lpbf_column.py: grains (colored by
          order parameter) and the temperature in a vertical section, per snapshot.

Needs matplotlib and Pillow. Colors: temperature on a perceptually uniform map with a
fixed range over all frames, so that frames can be compared.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def lattice(nodes: np.ndarray):
    axes = [np.unique(np.round(nodes[:, d], 12)) for d in range(nodes.shape[1])]
    order = np.lexsort(tuple(nodes[:, d] for d in range(nodes.shape[1])))
    shape = tuple(a.size for a in reversed(axes))
    return axes, order, shape


def save_gif(frames: list[Path], out: Path, ms: int = 180) -> None:
    images = [Image.open(f).convert("P", palette=Image.ADAPTIVE) for f in frames]
    images[0].save(out, save_all=True, append_images=images[1:], duration=ms, loop=0)
    print("wrote", out, f"({len(images)} frames; PNG stack in {frames[0].parent})")


def cmd_track(run_dir: str, out: str, liquidus: float = 1610.15, solidus: float = 1528.15) -> None:
    data = np.load(Path(run_dir) / "fields.npz")
    nodes = data["nodes"]
    (x, y, z), order, shape = lattice(nodes)
    times = data["snapshot_times"]
    raw = data["snapshot_temperature"]
    temp = raw.reshape(raw.shape[0], -1)[:, order].reshape((-1,) + shape)
    has_flow = "snapshot_velocity" in data.files
    if has_flow:
        vel = data["snapshot_velocity"][:, order, :].reshape((-1,) + shape + (3,))
    jc = y.size // 2
    um = 1e6
    t_hi = 3200.0
    frame_dir = Path(out).with_suffix("")
    frame_dir = frame_dir.parent / (frame_dir.name + "_frames")
    frame_dir.mkdir(parents=True, exist_ok=True)
    files = []
    zlo = max(0, z.size - 22)
    for k, t in enumerate(times):
        fig, ax = plt.subplots(2, 1, figsize=(8.6, 6.4), gridspec_kw={"height_ratios": [1, 1.25]})
        side = temp[k, zlo:, jc, :]
        im = ax[0].pcolormesh(x * um, z[zlo:] * um, side, cmap="inferno", vmin=300.0, vmax=t_hi, shading="gouraud")
        ax[0].contour(x * um, z[zlo:] * um, side, levels=[solidus, liquidus], colors=["c", "w"], linewidths=0.9)
        top = temp[k, -1, :, :]
        ax[1].pcolormesh(x * um, y * um, top, cmap="inferno", vmin=300.0, vmax=t_hi, shading="gouraud")
        ax[1].contour(x * um, y * um, top, levels=[solidus, liquidus], colors=["c", "w"], linewidths=0.9)
        title = f"t = {t * 1e6:6.0f} us"
        if has_flow:
            speed = np.linalg.norm(vel[k], axis=-1)
            title += f"   largest speed {speed.max():.2f} m/s"
            s = 2
            u, w = vel[k, zlo:, jc, :, 0], vel[k, zlo:, jc, :, 2]
            mask = side >= solidus
            if mask.any():
                ax[0].quiver(x[::s] * um, z[zlo:][::s] * um, np.where(mask, u, np.nan)[::s, ::s],
                             np.where(mask, w, np.nan)[::s, ::s], color="w", scale=30.0, width=0.003)
            ut, vt = vel[k, -1, :, :, 0], vel[k, -1, :, :, 1]
            mask = top >= solidus
            if mask.any():
                ax[1].quiver(x[::s] * um, y[::s] * um, np.where(mask, ut, np.nan)[::s, ::s],
                             np.where(mask, vt, np.nan)[::s, ::s], color="w", scale=30.0, width=0.003)
        ax[0].set_title(title + "   (section along the track axis)", fontsize=9)
        ax[0].set_xlabel("x, scan direction ($\\mu$m)"); ax[0].set_ylabel("z ($\\mu$m)")
        ax[0].set_aspect("equal")
        ax[1].set_title("top surface (white: liquidus, cyan: solidus)", fontsize=9)
        ax[1].set_xlabel("x, scan direction ($\\mu$m)"); ax[1].set_ylabel("y ($\\mu$m)")
        ax[1].set_aspect("equal")
        fig.colorbar(im, ax=ax, shrink=0.8, label="temperature (K), color range capped at 3,200 K")
        path = frame_dir / f"frame_{k:03d}.png"
        fig.savefig(path, dpi=110)
        plt.close(fig)
        files.append(path)
    save_gif(files, Path(out))


def cmd_grains(run_dir: str, out: str) -> None:
    data = np.load(Path(run_dir) / "fields.npz")
    nodes = data["nodes"]
    (x, y, z), order, shape = lattice(nodes)
    times = data["snapshot_times"]
    if times.size == 0:
        raise SystemExit("no snapshots in this run (set snapshots_every in [Executioner])")
    grain = data["snapshot_grain"][:, order].reshape((-1,) + shape).astype(float)
    solid = data["snapshot_solid"][:, order].reshape((-1,) + shape).astype(float)
    temp = data["snapshot_temperature"][:, order].reshape((-1,) + shape)
    jc = y.size // 2
    um = 1e6
    frame_dir = Path(out).with_suffix("")
    frame_dir = frame_dir.parent / (frame_dir.name + "_frames")
    frame_dir.mkdir(parents=True, exist_ok=True)
    files = []
    for k, t in enumerate(times):
        fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6))
        shown = np.where(solid[k, :, jc, :] > 0.25, grain[k, :, jc, :], np.nan)
        ax[0].pcolormesh(x * um, z * um, shown, cmap="tab20", vmin=0, vmax=max(grain.max(), 1), shading="nearest")
        ax[0].set_title(f"grains (order parameter), t = {t * 1e3:6.2f} ms", fontsize=9)
        im = ax[1].pcolormesh(x * um, z * um, temp[k, :, jc, :], cmap="inferno", vmin=300.0, vmax=3200.0, shading="gouraud")
        ax[1].contour(x * um, z * um, temp[k, :, jc, :], levels=[1610.15], colors="w", linewidths=0.9)
        ax[1].set_title("temperature (K), white: liquidus", fontsize=9)
        for a in ax:
            a.set_xlabel("x ($\\mu$m)"); a.set_ylabel("z, build direction ($\\mu$m)"); a.set_aspect("equal")
        fig.colorbar(im, ax=ax[1], shrink=0.8)
        path = frame_dir / f"frame_{k:03d}.png"
        fig.savefig(path, dpi=110)
        plt.close(fig)
        files.append(path)
    save_gif(files, Path(out))


def main() -> None:
    if len(sys.argv) >= 4 and sys.argv[1] == "track":
        cmd_track(sys.argv[2], sys.argv[3], *[float(v) for v in sys.argv[4:6]])
    elif len(sys.argv) == 4 and sys.argv[1] == "grains":
        cmd_grains(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
