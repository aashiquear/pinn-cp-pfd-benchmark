#!/usr/bin/env python3
"""Part V: animations from the field series of a process run (scripts/lpbf_run.py).

usage: scripts/lpbf_animate.py track <run_dir> <out.gif> [liquidus_K] [solidus_K]
       scripts/lpbf_animate.py grains <column_run_dir> <out.gif>
       scripts/lpbf_animate.py column <column_run_dir> <out_prefix> [keep_every]

  track   one frame per stored snapshot of fields.npz: the temperature in the plane along
          the track (through the track axis) and in the top layer of nodes, with the
          liquidus and solidus isotherms; when the run has a velocity field, arrows of
          the in-plane velocity and the largest speed. The frames are also written as
          PNG files next to the GIF (<out>_frames/), the stack the figures of the report
          are taken from.
  grains  the same for a build-column run of scripts/lpbf_column.py: grains (colored by
          order parameter) and the temperature in a vertical section, per snapshot.

  column  two animations of a build-column run, <out_prefix>_2d.gif and
          <out_prefix>_3d.gif. 2-D: grains (inverse pole figure colors of the build
          direction; white: liquid; gray: powder or empty) and temperature, each in a
          vertical section and in the plane of the current build surface; arrows of the
          velocity when the run has a flow. 3-D: the same two fields on the faces of the
          block in an isometric view, cut at the edges of the scanned patch and at the
          current build surface. All frames are kept in the run directory; every
          `keep_every`-th (default 10) is copied to <out_prefix>_{2d,3d}_frames/.

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


def cmd_column(run_dir: str, prefix: str, keep_every: int = 10, liquidus: float = 1610.15,
               plate_top: float = 200.0e-6, crop: tuple[float, float] = (120.0e-6, 480.0e-6)) -> None:
    import shutil

    run = Path(run_dir)
    data = np.load(run / "fields.npz")
    nodes = data["nodes"]
    (x, y, z), order, shape = lattice(nodes)
    times = data["snapshot_times"]
    if times.size == 0:
        raise SystemExit("no snapshots in this run (set snapshots_every in [Executioner])")
    grain = data["snapshot_grain"][:, order].reshape((-1,) + shape)
    solid = data["snapshot_solid"][:, order].reshape((-1,) + shape).astype(np.float32)
    temp = data["snapshot_temperature"][:, order].reshape((-1,) + shape)
    vel = None
    if "snapshot_velocity" in data.files:
        vel = data["snapshot_velocity"][:, order, :].reshape((-1,) + shape + (3,)).astype(np.float32)
    rot = data["rotations"]
    bd = np.sort(np.abs(rot[:, 2, :]), axis=1)
    rgb = np.stack([bd[:, 2] - bd[:, 1], bd[:, 1] - bd[:, 0], bd[:, 0] * np.sqrt(3.0)], axis=1)
    rgb /= rgb.max(axis=1, keepdims=True)
    um = 1e6
    t_lo, t_hi = 350.0, 3200.0
    cmap = plt.get_cmap("inferno")
    # index of the build surface per frame: the highest plane of nodes whose temperature has moved
    moved = np.abs(temp - temp[0, -1].min()) > 1.0e-2
    planes = moved.reshape(moved.shape[0], shape[0], -1).any(axis=2)
    top_now = np.array([np.nonzero(p)[0].max() if p.any() else 0 for p in planes])
    top_index = np.maximum.accumulate(np.maximum(top_now, int(np.searchsorted(z, plate_top - 1e-9))))
    jc, ic = y.size // 2, x.size // 2
    i0, i1 = np.searchsorted(x, crop[0] - 1e-9), np.searchsorted(x, crop[1] - 1e-9)
    j0, j1 = np.searchsorted(y, crop[0] - 1e-9), np.searchsorted(y, crop[1] - 1e-9)

    def grain_image(k: int, index: tuple) -> np.ndarray:
        image = rgb[grain[(k,) + index]].copy()
        image[solid[(k,) + index] < 0.25] = 0.80          # powder or empty
        image[temp[(k,) + index] >= liquidus] = 1.0       # liquid
        return image

    def heat_image(k: int, index: tuple) -> np.ndarray:
        return cmap(np.clip((temp[(k,) + index] - t_lo) / (t_hi - t_lo), 0.0, 1.0))[..., :3]

    dirs = {tag: run / f"animation_frames_{tag}" for tag in ("2d", "3d")}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    files: dict[str, list[Path]] = {"2d": [], "3d": []}
    scalar = plt.cm.ScalarMappable(norm=plt.Normalize(t_lo, t_hi), cmap=cmap)
    for k, t in enumerate(times):
        kt = int(top_index[k])
        section, plane = (slice(None), jc, slice(None)), (kt, slice(None), slice(None))
        # ---------------- 2-D ----------------
        fig, ax = plt.subplots(2, 2, figsize=(11.0, 8.0), layout="constrained",
                               gridspec_kw={"width_ratios": [1.45, 1.0]})
        ext_xz = [x[0] * um, x[-1] * um, z[0] * um, z[-1] * um]
        ext_xy = [x[0] * um, x[-1] * um, y[0] * um, y[-1] * um]
        ax[0, 0].imshow(grain_image(k, section), origin="lower", extent=ext_xz, interpolation="nearest")
        ax[0, 1].imshow(grain_image(k, plane), origin="lower", extent=ext_xy, interpolation="nearest")
        ax[1, 0].imshow(heat_image(k, section), origin="lower", extent=ext_xz, interpolation="bilinear")
        ax[1, 1].imshow(heat_image(k, plane), origin="lower", extent=ext_xy, interpolation="bilinear")
        for a_, index, axes_ in ((ax[0, 0], section, (x, z)), (ax[1, 0], section, (x, z)),
                                 (ax[0, 1], plane, (x, y)), (ax[1, 1], plane, (x, y))):
            field = temp[(k,) + index]
            if field.max() >= liquidus:
                a_.contour(axes_[0] * um, axes_[1] * um, field, levels=[liquidus],
                           colors="k" if a_ in (ax[0, 0], ax[0, 1]) else "w", linewidths=0.9)
        note = ""
        if vel is not None:
            speed = np.linalg.norm(vel[k], axis=-1)
            note = f"   largest speed {speed.max():.2f} m/s"
            u, w = vel[k, :, jc, :, 0], vel[k, :, jc, :, 2]
            mask = temp[k, :, jc, :] >= liquidus
            if mask.any():
                ax[1, 0].quiver(x * um, z * um, np.where(mask, u, np.nan), np.where(mask, w, np.nan),
                                color="c", scale=40.0, width=0.003)
            u, v = vel[k, kt, :, :, 0], vel[k, kt, :, :, 1]
            mask = temp[k, kt] >= liquidus
            if mask.any():
                ax[1, 1].quiver(x * um, y * um, np.where(mask, u, np.nan), np.where(mask, v, np.nan),
                                color="c", scale=40.0, width=0.003)
        for a_ in (ax[0, 0], ax[1, 0]):
            a_.axhline(plate_top * um, color="k" if a_ is ax[0, 0] else "w", lw=0.6, ls="--")
            a_.set_xlabel("x ($\\mu$m)"); a_.set_ylabel("z, build direction ($\\mu$m)")
        for a_ in (ax[0, 1], ax[1, 1]):
            a_.set_xlabel("x ($\\mu$m)"); a_.set_ylabel("y ($\\mu$m)")
        ax[0, 0].set_title(f"grains, section y = {y[jc] * um:.0f} $\\mu$m", fontsize=9)
        ax[0, 1].set_title(f"grains, build surface z = {z[kt] * um:.0f} $\\mu$m", fontsize=9)
        ax[1, 0].set_title("temperature, same section", fontsize=9)
        ax[1, 1].set_title("temperature, build surface", fontsize=9)
        fig.suptitle(f"t = {t * 1e3:6.2f} ms{note}    grains: inverse pole figure color of the build direction; "
                     "white: liquid; gray: powder or empty; line: liquidus", fontsize=9)
        fig.colorbar(scalar, ax=ax, shrink=0.6, label="temperature (K), color range capped at 3,200 K")
        path = dirs["2d"] / f"frame_{k:03d}.png"
        fig.savefig(path, dpi=90)
        plt.close(fig)
        files["2d"].append(path)
        # ---------------- 3-D ----------------
        fig = plt.figure(figsize=(11.0, 5.4))
        xs, ys, zs = x[i0:i1 + 1] * um, y[j0:j1 + 1] * um, z[:kt + 1] * um
        for col, (image_of, title) in enumerate(((grain_image, "grains"), (heat_image, "temperature"))):
            a3 = fig.add_subplot(1, 2, col + 1, projection="3d")
            top = image_of(k, (kt, slice(j0, j1 + 1), slice(i0, i1 + 1)))
            front = image_of(k, (slice(0, kt + 1), j0, slice(i0, i1 + 1)))
            right = image_of(k, (slice(0, kt + 1), slice(j0, j1 + 1), i1))
            gx, gy = np.meshgrid(xs, ys)
            a3.plot_surface(gx, gy, np.full_like(gx, zs[-1]), facecolors=top[:-1, :-1], rstride=1, cstride=1,
                            shade=False, antialiased=False, linewidth=0)
            gx, gz = np.meshgrid(xs, zs)
            a3.plot_surface(gx, np.full_like(gx, ys[0]), gz, facecolors=front[:-1, :-1], rstride=1, cstride=1,
                            shade=False, antialiased=False, linewidth=0)
            gy, gz = np.meshgrid(ys, zs)
            a3.plot_surface(np.full_like(gy, xs[-1]), gy, gz, facecolors=right[:-1, :-1], rstride=1, cstride=1,
                            shade=False, antialiased=False, linewidth=0)
            a3.set_xlim(xs[0], xs[-1]); a3.set_ylim(ys[0], ys[-1]); a3.set_zlim(0.0, z[-1] * um)
            a3.set_box_aspect((xs[-1] - xs[0], ys[-1] - ys[0], z[-1] * um))
            a3.view_init(elev=24.0, azim=-52.0)
            a3.set_xlabel("x ($\\mu$m)"); a3.set_ylabel("y ($\\mu$m)"); a3.set_zlabel("z ($\\mu$m)")
            a3.set_title(title, fontsize=10)
        fig.suptitle(f"t = {t * 1e3:6.2f} ms{note}    block cut at the edges of the scanned patch and at the "
                     f"build surface (z = {z[kt] * um:.0f} $\\mu$m)", fontsize=9)
        fig.colorbar(scalar, ax=fig.axes, shrink=0.7, label="temperature (K)")
        path = dirs["3d"] / f"frame_{k:03d}.png"
        fig.savefig(path, dpi=90)
        plt.close(fig)
        files["3d"].append(path)
    for tag, paths in files.items():
        out = Path(f"{prefix}_{tag}.gif")
        out.parent.mkdir(parents=True, exist_ok=True)
        save_gif(paths, out, ms=120)
        kept = out.parent / f"{out.stem}_frames"
        kept.mkdir(exist_ok=True)
        for path in paths[::keep_every]:
            shutil.copy2(path, kept / path.name)


def main() -> None:
    if len(sys.argv) >= 4 and sys.argv[1] == "track":
        cmd_track(sys.argv[2], sys.argv[3], *[float(v) for v in sys.argv[4:6]])
    elif len(sys.argv) == 4 and sys.argv[1] == "grains":
        cmd_grains(sys.argv[2], sys.argv[3])
    elif len(sys.argv) in (4, 5) and sys.argv[1] == "column":
        cmd_column(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) == 5 else 10)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
