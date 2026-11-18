#!/usr/bin/env python3
"""Case F4 (quaternion orientation field in three dimensions): collection, measures, figures.

Usage: python3 scripts/cryst_quaternion_figures.py <tag> <run dir> [<run dir> ...]
  run dirs relative to the repo, each with the VTK series of a PINNACLE run (phi, c, q).

Per run: frames.npz (phi, c, q at the saved steps, float32) and a row of
results/crystallization/quaternion/<tag>_measures.csv:
  solid fraction, equivalent radius (cells) of the solid volume, sphericity of the solid
  (surface of the equal-volume sphere / counted surface of the phi = 0.5 voxels x 2/3),
  mean misorientation angle between neighboring crystal cells (degrees, 2 acos|q.q'|),
  boundary fraction (neighbor pairs above 15 degrees), number of orientation domains
  (connected crystal regions whose neighbor misorientation stays below 15 degrees, > 30 cells).
Figures: figures/crystallization/<tag>_views.png: for every run a surface view along z (height
of the topmost crystal cell as shading, colored by the tilt of the crystal z axis against the
laboratory z axis, 0-90 degrees, as in the paper's figures) and the mid-plane sections of the
composition and of the tilt; <tag>_<run>.gif: the surface view over time.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from scipy import ndimage  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk, series_files  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures/crystallization"
TILT_CMAP = plt.get_cmap("turbo")


def collect(run: Path) -> dict:
    cache = run / "frames.npz"
    if cache.exists():
        return dict(np.load(cache))
    files = series_files(run)
    first = read_vtk(files[0][1])
    nodes = first["nodes"]
    axes = [np.unique(np.round(nodes[:, a], 9)) for a in range(3)]
    idx = [np.searchsorted(axes[a], np.round(nodes[:, a], 9)) for a in range(3)]
    shape = tuple(len(a) for a in axes)
    out: dict[str, list] = {"phi": [], "c": [], "q": []}
    for _, path in files:
        data = read_vtk(path)["point"]

        def cube(key: str) -> np.ndarray:
            a = np.zeros(shape, dtype=np.float32)
            a[idx[0], idx[1], idx[2]] = data[key]
            return a[:-1, :-1, :-1]  # drop the periodic images

        out["phi"].append(cube("field_phi"))
        out["c"].append(cube("field_c"))
        out["q"].append(np.stack([cube(f"field_q_{i}") for i in range(4)], axis=-1))
    hist = list(csv.DictReader((run / "solver_history.csv").open()))
    n_steps = len(hist)
    steps = np.linspace(0, n_steps, len(files)).round().astype(int)
    result = {"phi": np.array(out["phi"]), "c": np.array(out["c"]), "q": np.array(out["q"]), "steps": steps}
    np.savez_compressed(cache, **result)
    return result


def tilt(q: np.ndarray) -> np.ndarray:
    """Angle (degrees, folded to 0-90) between the crystal and laboratory z axes."""
    r33 = q[..., 0] ** 2 - q[..., 1] ** 2 - q[..., 2] ** 2 + q[..., 3] ** 2
    angle = np.degrees(np.arccos(np.clip(r33, -1.0, 1.0)))
    return np.minimum(angle, 180.0 - angle)


def misorientation(q: np.ndarray, solid: np.ndarray) -> np.ndarray:
    """Angles (degrees) between periodic neighbors that are both crystal, all three axes."""
    out = []
    for axis in range(3):
        dot = np.abs(np.sum(q * np.roll(q, -1, axis=axis), axis=-1))
        both = solid & np.roll(solid, -1, axis=axis)
        out.append(np.degrees(2.0 * np.arccos(np.clip(dot[both], 0.0, 1.0))))
    return np.concatenate(out)


def domains(q: np.ndarray, solid: np.ndarray, limit: float = 15.0, minimum: int = 30) -> int:
    """Connected crystal regions bounded by neighbor misorientations above ``limit``."""
    interior = solid.copy()
    for axis in range(3):
        dot = np.abs(np.sum(q * np.roll(q, -1, axis=axis), axis=-1))
        sharp = np.degrees(2.0 * np.arccos(np.clip(dot, 0.0, 1.0))) > limit
        interior &= ~sharp
        interior &= ~np.roll(sharp, 1, axis=axis)
    labels, count = ndimage.label(interior)
    if count == 0:
        return 0
    sizes = ndimage.sum(np.ones_like(labels), labels, index=np.arange(1, count + 1))
    return int(np.sum(sizes > minimum))


def surface_view(phi: np.ndarray, q: np.ndarray) -> np.ndarray:
    """RGB image looking down the z axis: tilt color of the topmost crystal cell, depth shading."""
    solid = phi < 0.5
    n = phi.shape[2]
    top = np.where(solid.any(axis=2), n - 1 - np.argmax(solid[:, :, ::-1], axis=2), -1)
    angle = tilt(q)
    ii, jj = np.meshgrid(np.arange(phi.shape[0]), np.arange(phi.shape[1]), indexing="ij")
    color = TILT_CMAP(angle[ii, jj, np.clip(top, 0, n - 1)] / 90.0)[..., :3]
    height = np.where(top >= 0, top, 0).astype(float)
    gx, gy = np.gradient(ndimage.gaussian_filter(height, 0.7))
    shade = np.clip(0.75 + 0.12 * (gx - gy), 0.35, 1.0) * (0.55 + 0.45 * height / max(height.max(), 1.0))
    image = color * shade[..., None]
    image[top < 0] = 0.0
    return np.transpose(image, (1, 0, 2))


def main() -> None:
    tag = sys.argv[1]
    runs = [ROOT / a for a in sys.argv[2:]]
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    rows = []
    fig, axes = plt.subplots(3, len(runs), figsize=(3.3 * len(runs), 9.4), squeeze=False)
    for j, run in enumerate(runs):
        d = collect(run)
        phi, c, q = d["phi"][-1].astype(float), d["c"][-1].astype(float), d["q"][-1].astype(float)
        solid = phi < 0.5
        volume = float(solid.sum())
        faces = sum(float(np.sum(solid != np.roll(solid, 1, axis=a))) for a in range(3))
        radius = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0)
        mis = misorientation(q, solid)
        rows.append({"run": run.name, "step": int(d["steps"][-1]), "solid_fraction": volume / phi.size,
                     "equivalent_radius": radius,
                     "sphericity": 4.0 * np.pi * radius**2 / (faces * 2.0 / 3.0) if faces else float("nan"),
                     "mean_misorientation_deg": float(mis.mean()) if mis.size else float("nan"),
                     "boundary_fraction": float((mis > 15.0).mean()) if mis.size else float("nan"),
                     "orientation_domains": domains(q, solid)})
        print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in rows[-1].items()))
        mid = phi.shape[2] // 2
        axes[0, j].imshow(surface_view(phi, q), origin="lower")
        axes[0, j].set_title(f"{run.name}\nstep {int(d['steps'][-1])}", fontsize=7)
        im1 = axes[1, j].imshow(c[:, :, mid].T, origin="lower", cmap="viridis", vmin=0.394, vmax=0.471)
        section = np.ma.masked_where(~solid[:, :, mid].T, tilt(q)[:, :, mid].T)
        axes[2, j].set_facecolor("k")
        im2 = axes[2, j].imshow(section, origin="lower", cmap=TILT_CMAP, vmin=0, vmax=90)
        for a in axes[:, j]:
            a.set_xticks([]); a.set_yticks([])
    axes[0, 0].set_ylabel("view along $z$")
    axes[1, 0].set_ylabel("mid-plane, composition $c$")
    axes[2, 0].set_ylabel("mid-plane, tilt of the crystal $z$ axis")
    fig.colorbar(im1, ax=axes[1, :], shrink=0.8, pad=0.01)
    fig.colorbar(im2, ax=axes[2, :], shrink=0.8, pad=0.01, label="degrees")
    fig.savefig(FIG / f"{tag}_views.png", bbox_inches="tight")
    plt.close(fig)
    out = ROOT / "results/crystallization/quaternion"
    with (out / f"{tag}_measures.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for run in runs:
        d = collect(run)
        fig, ax = plt.subplots(figsize=(3.4, 3.5))
        im = ax.imshow(surface_view(d["phi"][0].astype(float), d["q"][0].astype(float)), origin="lower")
        ax.set_xticks([]); ax.set_yticks([])
        title = ax.set_title("")

        def draw(i, d=d, im=im, title=title):
            im.set_data(surface_view(d["phi"][i].astype(float), d["q"][i].astype(float)))
            title.set_text(f"step {int(d['steps'][i])}")
            return im, title

        FuncAnimation(fig, draw, frames=len(d["steps"])).save(FIG / f"{tag}_{run.name}.gif",
                                                              writer=PillowWriter(fps=3), dpi=80)
        plt.close(fig)
    print("wrote", f"figures/crystallization/{tag}_views.png")


if __name__ == "__main__":
    main()
