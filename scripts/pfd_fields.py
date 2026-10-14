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
  [Runs]
    [./cp]  dir = ../../../results/pfd_study/cp/cal_h40/l5/gc1106.12  label = 'CP+PFD, l0 5' []
    [./j2]  dir = ../../../results/pfd_study/j2/cal_h40/l5/gc2063.43  label = 'J2+PFD, l0 5' []

For every run (``results.csv`` for the strain of each step,
``solver/fields_step*.vtk`` for the fields) and every requested strain the
script draws the nodal damage on the deformed plate (displacement scaled by
``deform``, default 1) and, below it, the second field; the crack path at
the end of the path (elements whose mean damage is at or above the
threshold) is drawn as a separate figure with the damage contour of the
last snapshot. Writes ``<prefix>_damage.png`` (rows: runs; columns:
strains), ``<prefix>_<field>.png`` and ``<prefix>_crack.png``, and
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
    runs = [(n, _resolve(str(r["dir"]), deck), _str(r.get("label", n))) for n, r in cfg["runs"].items()]
    plt.rcParams.update(STYLE)
    records = []
    n_r, n_c = len(runs), len(strains)
    fig_d, axes_d = plt.subplots(n_r, n_c, figsize=(3.2 * n_c, 1.75 * n_r + 0.5), dpi=200, squeeze=False)
    fig_p, axes_p = plt.subplots(n_r, n_c, figsize=(3.2 * n_c, 1.75 * n_r + 0.5), dpi=200, squeeze=False)
    fig_c, axes_c = plt.subplots(1, n_r, figsize=(4.2 * n_r, 2.4), dpi=200, squeeze=False)
    for i, (name, run, label) in enumerate(runs):
        snaps = _snapshots(run)
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
            pc = PolyCollection(quads, array=d_cell, cmap="inferno", edgecolors="none", clim=(0, 1))
            ax.add_collection(pc)
            ax.set_aspect("equal"); ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"{label}: {100 * s_act:.1f}%  max d {np.nanmax(d_cell):.2f}", fontsize=7)
            p_cell = _cell_mean(v, field2, cells)
            ax2 = axes_p[i, j]
            if p_cell is not None:
                pc2 = PolyCollection(quads, array=p_cell, cmap="viridis", edgecolors="none", clim=(0, vmax2))
                ax2.add_collection(pc2)
            ax2.set_aspect("equal"); ax2.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax2.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
            ax2.set_xticks([]); ax2.set_yticks([])
            ax2.set_title(f"{label}: {100 * s_act:.1f}%" + (f"  max {np.nanmax(p_cell):.2f}" if p_cell is not None else ""), fontsize=7)
            records.append({"run": name, "strain_pct": 100 * s_act, "max_damage": float(np.nanmax(d_cell)), "mean_damage": float(np.nanmean(d_cell)),
                            "elements_above": int(np.sum(d_cell >= thr)), f"max_{field2}": float(np.nanmax(p_cell)) if p_cell is not None else None,
                            "snapshot": path.name})
        # crack path at the last snapshot
        s_last, path = snaps[-1]
        v = read_vtk(path)
        quads, cells = _quads(v, deform)
        d_cell = _cell_mean(v, "damage", cells)
        ax = axes_c[0, i]
        ax.add_collection(PolyCollection(quads, array=d_cell, cmap="inferno", edgecolors="none", clim=(0, 1), alpha=0.9))
        crack = d_cell >= thr
        if crack.any():
            ax.add_collection(PolyCollection(quads[crack], facecolors="none", edgecolors="#00e5ff", linewidths=0.6))
        ax.set_aspect("equal"); ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max()); ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{label}: crack path at {100 * s_last:.1f}% ({int(crack.sum())} elements, d >= {thr:g})", fontsize=7)
    for fig, nm in ((fig_d, "damage"), (fig_p, field2), (fig_c, "crack")):
        fig.tight_layout(); fig.savefig(figs / f"{prefix}_{nm}.png"); plt.close(fig)
    with (figs.parent / "results" / "pfd_study" / f"{prefix}_fields.csv").open("w", newline="") if False else (STUDY_ROOT / "results" / "pfd_study" / f"{prefix}_fields.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    for r in records:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", figs / f"{prefix}_{{damage,{field2},crack}}.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
