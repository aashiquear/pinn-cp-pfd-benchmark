#!/usr/bin/env python3
"""Damage and plastic-strain maps of PFD runs on the perforated plate.

Deck-driven (PINNACLE ``.i`` syntax):

  [Fields]
    figures    = ../../../figures/pfd
    prefix     = cp_cal_h40
    strains    = '5 10 15 17.2 19'     # nominal strains (%) of the maps (nearest snapshot)
    strain_max = 0.20
    threshold  = 0.8                   # crack path: elements with mean damage >= threshold
    field      = equivalent_plastic_strain   # second row (cell or point field of the VTK); optional
    vmax       = 0.6
    gif        = true                        # animation of the second field with the cracked elements (d >= threshold) solid
    gif_every  = 1                           # every k-th snapshot
    columns    = 4                           # crack figure and GIF: runs per row
  [Runs]
    [./cp]  dir = ../../../results/pfd_study/cp/cal_h40/l5/gc1106.12  label = 'CP+PFD, l0 5' []
    [./j2]  dir = ../../../results/pfd_study/j2/cal_h40/l5/gc2063.43  label = 'J2+PFD, l0 5' []

For every run (``results.csv`` for the strain of each step,
``solver/fields_step*.vtk`` for the fields) and every requested strain the
script draws the nodal damage on the deformed plate (displacement scaled by
``deform``, default 1) and, in a second figure, the second field, each with
a colour bar; the crack figure shows the second field at the end of the
path with the cracked elements (mean damage at or above the threshold)
filled solid. With ``gif = true`` the same crack view of every snapshot is
assembled into ``<prefix>_evolution.gif`` (rows: runs). Writes
``<prefix>_damage.png`` (rows: runs; columns: strains),
``<prefix>_<field>.png``, ``<prefix>_crack.png`` and
``<prefix>_fields.csv`` (per run and strain: max and mean damage, number of
elements above the threshold, max of the second field).

Usage: pfd_fields.py decks/pfd_study/report/fields_cal_h40.i
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(STUDY_ROOT / "scripts"))
from vtk_legacy import read_vtk, series_files  # noqa: E402

STYLE = {"font.size": 8, "axes.labelsize": 9, "mathtext.fontset": "dejavusans"}


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("fields", "runs"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def _levels(v) -> list[float]:
    return [float(x) for x in (v if isinstance(v, (list, tuple)) else str(v).split())]


def _snapshots(run: Path):
    """(strain, vtk path) per snapshot: the snapshot index maps to the step through snapshots_every."""
    rows = list(csv.DictReader((run / "results.csv").open()))
    strains = [float(r["nominal_strain"]) for r in rows]
    files = series_files(run / "solver")
    if not files:
        return []
    # snapshots are written every k steps: k = n_steps / n_snapshots (rounded)
    k = max(1, round(len(strains) / len(files)))
    out = []
    for i, (_, f) in enumerate(files):
        step = min((i + 1) * k - 1, len(strains) - 1)
        out.append((strains[step], f))
    return out


def _quads(v, deform: float):
    pts = np.asarray(v["nodes"], dtype=float)
    disp = v["point"].get("displacement")
    if disp is not None and deform:
        pts = pts + deform * np.asarray(disp, dtype=float)[:, : pts.shape[1]]
    cells = np.asarray(v["cells"], dtype=int)
    return pts[cells[:, :4], :2], cells


def _cell_mean(v, name: str, cells: np.ndarray):
    if name in v["cell"]:
        return np.asarray(v["cell"][name], dtype=float).reshape(len(cells), -1)[:, 0]
    arr = v["point"].get(name)
    if arr is None:
        return None
    arr = np.asarray(arr, dtype=float).reshape(len(v["nodes"]), -1)[:, 0]
    return arr[cells].mean(axis=1)


CRACK = "#ff2d2d"   # solid fill of the cracked elements


def _frame(ax, quads, d_cell, p_cell, thr, vmax2, field2, title):
    """Second field with the cracked elements solid; returns the field collection (for the colour bar)."""
    pc = None
    if p_cell is not None:
        pc = PolyCollection(quads, array=p_cell, cmap="viridis", edgecolors="none", clim=(0, vmax2))
        ax.add_collection(pc)
    else:
        ax.text(0.5, 0.5, f"{field2} not in the series", transform=ax.transAxes, ha="center", va="center", fontsize=7, color="#888888")
    crack = d_cell >= thr
    if crack.any():
        ax.add_collection(PolyCollection(quads[crack], facecolors=CRACK, edgecolors=CRACK, linewidths=0.3))
    ax.set_aspect("equal"); ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=7)
    return pc, int(crack.sum())


def _field_title(field2: str) -> str:
    return {"equivalent_plastic_strain": "equivalent plastic strain", "driving_H": "plastic work density (MPa)"}.get(field2, field2.replace("_", " "))


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    f = cfg["fields"]
    figs = _resolve(str(f["figures"]), deck)
    figs.mkdir(parents=True, exist_ok=True)
    prefix = str(f.get("prefix", "pfd"))
    strains = _levels(f.get("strains", "5 10 15 17.2"))
    thr = float(f.get("threshold", 0.8))
    field2 = str(f.get("field", "equivalent_plastic_strain"))
    vmax2 = float(f.get("vmax", 0.6))
    deform = float(f.get("deform", 1.0))
    gif = str(f.get("gif", "false")).lower() in ("true", "1", "yes")
    gif_every = int(f.get("gif_every", 1))
    ncol = max(1, min(int(f.get("columns", 4)), len(cfg["runs"])))
    nrow = -(-len(cfg["runs"]) // ncol)
    runs = [(n, _resolve(str(r["dir"]), deck), _str(r.get("label", n))) for n, r in cfg["runs"].items()]
    plt.rcParams.update(STYLE)
    records = []
    n_r, n_c = len(runs), len(strains)
    fig_d, axes_d = plt.subplots(n_r, n_c, figsize=(3.2 * n_c + 0.6, 1.75 * n_r + 0.5), dpi=200, squeeze=False, constrained_layout=True)
    fig_p, axes_p = plt.subplots(n_r, n_c, figsize=(3.2 * n_c + 0.6, 1.75 * n_r + 0.5), dpi=200, squeeze=False, constrained_layout=True)
    fig_c, axes_c = plt.subplots(nrow, ncol, figsize=(4.2 * ncol + 0.6, 2.3 * nrow + 0.2), dpi=200, squeeze=False, constrained_layout=True)
    axes_c = axes_c.reshape(1, -1)
    for ax in axes_c[0, n_r:]:
        ax.axis("off")
    pc_d = pc_p = pc_c = None
    all_snaps = []
    for i, (name, run, label) in enumerate(runs):
        snaps = _snapshots(run)
        all_snaps.append(snaps)
        if not snaps:
            print("no series in", run)
            continue
        for j, e_pct in enumerate(strains):
            k = int(np.argmin([abs(100 * s - e_pct) for s, _ in snaps]))
            s_act, path = snaps[k]
            v = read_vtk(path)
            quads, cells = _quads(v, deform)
            d_cell = _cell_mean(v, "damage", cells)
            ax = axes_d[i, j]
            pc_d = PolyCollection(quads, array=d_cell, cmap="inferno", edgecolors="none", clim=(0, 1))
            ax.add_collection(pc_d)
            ax.set_aspect("equal"); ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"{label}: {100 * s_act:.1f}%  max d {np.nanmax(d_cell):.2f}", fontsize=7)
            p_cell = _cell_mean(v, field2, cells)
            ax2 = axes_p[i, j]
            if p_cell is not None:
                pc_p = PolyCollection(quads, array=p_cell, cmap="viridis", edgecolors="none", clim=(0, vmax2))
                ax2.add_collection(pc_p)
            else:
                ax2.text(0.5, 0.5, f"{field2} not in the series", transform=ax2.transAxes, ha="center", va="center", fontsize=7, color="#888888")
            ax2.set_aspect("equal"); ax2.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax2.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
            ax2.set_xticks([]); ax2.set_yticks([])
            ax2.set_title(f"{label}: {100 * s_act:.1f}%" + (f"  max {np.nanmax(p_cell):.2f}" if p_cell is not None else ""), fontsize=7)
            records.append({"run": name, "strain_pct": 100 * s_act, "max_damage": float(np.nanmax(d_cell)), "mean_damage": float(np.nanmean(d_cell)),
                            "elements_above": int(np.sum(d_cell >= thr)), f"max_{field2}": float(np.nanmax(p_cell)) if p_cell is not None else None,
                            "snapshot": path.name})
        # crack path at the last snapshot: the second field with the cracked elements solid
        s_last, path = snaps[-1]
        v = read_vtk(path)
        quads, cells = _quads(v, deform)
        d_cell = _cell_mean(v, "damage", cells); p_cell = _cell_mean(v, field2, cells)
        n_crack = int((d_cell >= thr).sum())
        pc, _ = _frame(axes_c[0, i], quads, d_cell, p_cell, thr, vmax2, field2, f"{label}: {100 * s_last:.1f}%, {n_crack} elements with d >= {thr:g} (solid)")
        pc_c = pc or pc_c
    if pc_d is not None:
        fig_d.colorbar(pc_d, ax=axes_d.ravel().tolist(), shrink=0.8, pad=0.01, label="damage $d$ (element mean)")
    if pc_p is not None:
        fig_p.colorbar(pc_p, ax=axes_p.ravel().tolist(), shrink=0.8, pad=0.01, label=_field_title(field2))
    if pc_c is not None:
        fig_c.colorbar(pc_c, ax=axes_c.ravel().tolist(), shrink=0.8, pad=0.01, label=_field_title(field2))
    for fig, nm in ((fig_d, "damage"), (fig_p, field2), (fig_c, "crack")):
        fig.savefig(figs / f"{prefix}_{nm}.png"); plt.close(fig)
    if gif:
        from PIL import Image
        frames = []
        n_frames = max(len(sn) for sn in all_snaps)
        for k in range(0, n_frames, gif_every):
            fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol + 0.6, 2.3 * nrow + 0.2), dpi=120, squeeze=False, constrained_layout=True)
            axes = axes.reshape(1, -1)
            for ax in axes[0, n_r:]:
                ax.axis("off")
            pc_any = None
            for i, (name, run, label) in enumerate(runs):
                snaps = all_snaps[i]
                kk = min(k, len(snaps) - 1)   # a run that stopped early holds its last frame
                s_act, path = snaps[kk]
                v = read_vtk(path)
                quads, cells = _quads(v, deform)
                d_cell = _cell_mean(v, "damage", cells); p_cell = _cell_mean(v, field2, cells)
                pc, n_crack = _frame(axes[0, i], quads, d_cell, p_cell, thr, vmax2, field2, f"{label}: {100 * s_act:.1f}%, max d {np.nanmax(d_cell):.2f}, {n_crack} cracked elements")
                pc_any = pc or pc_any
            if pc_any is not None:
                fig.colorbar(pc_any, ax=axes.ravel().tolist(), shrink=0.8, pad=0.01, label=_field_title(field2))
            fig.canvas.draw()
            frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
            plt.close(fig)
        gif_path = figs / f"{prefix}_evolution.gif"
        frames[0].save(gif_path, save_all=True, append_images=frames[1:], duration=200, loop=0)
        print("wrote", gif_path, f"({len(frames)} frames)")
    with (STUDY_ROOT / "results" / "pfd_study" / f"{prefix}_fields.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    for r in records:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", figs / f"{prefix}_{{damage,{field2},crack}}.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
