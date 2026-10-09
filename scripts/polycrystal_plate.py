#!/usr/bin/env python3
"""Polycrystal pseudo-3D plate pipeline: phase-field grains -> CPFE deck.

Every stage is driven by ONE pipeline deck (PINNACLE ``.i`` syntax, parsed by
PINNACLE's own parser) so a structure and its mechanics deck can be rebuilt
from versioned inputs alone. The stages, in order:

  label   read the final ``eta`` field of a PINNACLE grain-growth run
          (``fields.vtk``), label grains (argmax colour + connected
          components, since several non-adjacent grains share one order
          parameter), drop specks, and write ``grains.npz`` + morphology
          metrics (``morphology.json``, via scripts/grain_morphology.py).
  orient  draw one Bunge Euler triple per grain, uniform on SO(3) from a
          recorded seed, and write the PINNACLE point file
          (``orientations_points.csv``: x,y,z,phi1,Phi,phi2,grain_id -- one row
          per growth-grid node, at mid-thickness) and ``grains.csv``.
  deck    fill ``decks/polycrystal_plate/cpfe_plate_template.i`` with the
          pipeline deck's [Mechanics] values and the point-file path, and
          write ``cpfe.i`` next to the point file. Run it with
          ``pinnacle run cpfe.i``.
  curve   after the run: nominal stress-strain curve from
          ``solver_history.csv`` (``results.csv``).

Pipeline deck blocks (see decks/polycrystal_plate/*.i for a filled example):

  [Growth]        run_dir  = <dir holding fields.vtk of the growth run>
                  deck     = <the growth deck, recorded for provenance>
  [Labeling]      threshold = 0.5   min_cells = 4
  [Orientations]  seed = <int>
  [Mechanics]     out_dir, strain, t_end, n_steps, nx, ny, thickness,
                  newton_tol, max_cutbacks, and the CP material values
                  (C11 C12 C44 mu gamma0 nexp tau_c0 q_lat alpha0 burgers
                  rho0 k1 y_c)

Usage: polycrystal_plate.py DECK.i STAGE [STAGE ...]   (STAGE in label orient deck curve all)
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

STUDY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(STUDY / "scripts"))

TEMPLATE = STUDY / "decks/polycrystal_plate/cpfe_plate_template.i"
MECH_KEYS = {
    # pipeline key -> template placeholder
    "C11": "C11", "C12": "C12", "C44": "C44", "mu": "MU", "gamma0": "GAMMA0",
    "nexp": "NEXP", "tau_c0": "TAUC0", "q_lat": "QLAT", "alpha0": "ALPHA0",
    "burgers": "BURGERS", "rho0": "RHO0", "k1": "K1", "y_c": "YC",
    "strain": "STRAIN", "t_end": "TEND", "n_steps": "NSTEPS", "nx": "NX", "ny": "NY",
    "thickness": "THICKNESS", "newton_tol": "NEWTONTOL", "max_cutbacks": "MAXCUTBACKS",
}


# --------------------------------------------------------------------------- deck
def read_pipeline_deck(path: Path) -> dict:
    """Parse the pipeline deck with PINNACLE's .i parser (block -> dict)."""
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    ir = InputParser().parse(path.read_text())
    cfg = Translator().translate(ir)
    for key in ("growth", "labeling", "orientations", "mechanics"):
        if key not in cfg:
            raise SystemExit(f"{path}: pipeline deck is missing the [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


# --------------------------------------------------------------------------- vtk
def read_legacy_vtk_points_scalars(path: Path, prefix: str = "eta_") -> tuple[np.ndarray, np.ndarray]:
    """Read POINTS and the ``SCALARS <prefix>k`` point arrays of a legacy ASCII VTK.

    Returns ``(points (n, 3), eta (n, k))``. Only what PINNACLE's writer emits
    is handled: ASCII, one value per line after ``LOOKUP_TABLE``.
    """
    lines = path.read_text().splitlines()
    i = 0
    points = None
    arrays: dict[int, np.ndarray] = {}
    n_points = 0
    while i < len(lines):
        tok = lines[i].split()
        if tok and tok[0] == "POINTS":
            n_points = int(tok[1])
            vals = []
            j = i + 1
            while len(vals) < 3 * n_points:
                vals.extend(float(v) for v in lines[j].split())
                j += 1
            points = np.asarray(vals, dtype=float).reshape(n_points, 3)
            i = j
            continue
        if tok and tok[0] == "SCALARS" and tok[1].startswith(prefix):
            k = int(tok[1][len(prefix):])
            j = i + 1
            if lines[j].startswith("LOOKUP_TABLE"):
                j += 1
            vals = []
            while len(vals) < n_points:
                vals.extend(float(v) for v in lines[j].split())
                j += 1
            arrays[k] = np.asarray(vals, dtype=float)
            i = j
            continue
        i += 1
    if points is None or not arrays:
        raise SystemExit(f"{path}: no POINTS / {prefix}* scalars found")
    eta = np.stack([arrays[k] for k in sorted(arrays)], axis=1)
    return points, eta


def _grid_from_points(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Recover the structured (nx, ny) grid of a 2-D box mesh and the node order map."""
    xs = np.unique(np.round(points[:, 0], 9))
    ys = np.unique(np.round(points[:, 1], 9))
    ix = np.searchsorted(xs, np.round(points[:, 0], 9))
    iy = np.searchsorted(ys, np.round(points[:, 1], 9))
    if len(xs) * len(ys) != points.shape[0]:
        raise SystemExit("growth mesh is not a structured 2-D grid")
    return xs, ys, ix * len(ys) + iy  # flat index in (ix, iy) C order


# --------------------------------------------------------------------------- label
def label_grains(eta_grid: np.ndarray, threshold: float, min_cells: int) -> tuple[np.ndarray, int]:
    """Grain ids on the grid: connected components within each argmax colour."""
    col = np.argmax(eta_grid, axis=-1)
    val = np.max(eta_grid, axis=-1)
    gid = np.full(col.shape, -1, dtype=int)
    nxt = 0
    for k in range(eta_grid.shape[-1]):
        mask = (col == k) & (val >= threshold)
        if not mask.any():
            continue
        lab, n = ndimage.label(mask)
        for i in range(1, n + 1):
            sel = lab == i
            if sel.sum() < min_cells:
                continue
            gid[sel] = nxt
            nxt += 1
    # Fill unlabelled cells (diffuse boundaries, specks) from the nearest labelled cell.
    if (gid < 0).any():
        idx = ndimage.distance_transform_edt(gid < 0, return_distances=False, return_indices=True)
        gid = gid[tuple(idx)]
    return gid, nxt


def stage_label(cfg: dict, deck: Path, out: Path) -> dict:
    from grain_morphology import summarize

    growth = cfg["growth"]
    run_dir = _resolve(str(growth["run_dir"]), deck)
    vtk = run_dir / "fields.vtk"
    points, eta = read_legacy_vtk_points_scalars(vtk)
    xs, ys, flat = _grid_from_points(points)
    nx, ny = len(xs), len(ys)
    eta_grid = np.zeros((nx, ny, eta.shape[1]))
    eta_grid.reshape(nx * ny, -1)[flat] = eta
    lab = cfg["labeling"]
    gid, n = label_grains(eta_grid, float(lab.get("threshold", 0.5)), int(lab.get("min_cells", 4)))
    hx, hy = float(xs[1] - xs[0]), float(ys[1] - ys[0])
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "grains.npz", gid=gid, xs=xs, ys=ys)
    morph = summarize(gid, hx, hy, label=str(run_dir))
    areas = np.bincount(gid.ravel(), minlength=n) * hx * hy
    rec = {
        "growth_run_dir": str(run_dir),
        "growth_deck": str(_resolve(str(growth.get("deck", "")), deck)) if growth.get("deck") else None,
        "grid": [nx, ny], "h_um": [hx, hy], "n_order_params": int(eta.shape[1]),
        "n_grains": int(n),
        "grain_area_um2": {"min": float(areas.min()), "median": float(np.median(areas)),
                           "max": float(areas.max())},
        "equivalent_diameter_um_median": float(np.sqrt(4 * np.median(areas) / np.pi)),
        "morphology": morph,
    }
    (out / "label.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: v for k, v in rec.items() if k != "morphology"}, indent=1))
    return rec


# --------------------------------------------------------------------------- orient
def uniform_so3_euler(n: int, seed: int) -> np.ndarray:
    """Uniform random orientations as Bunge (phi1, Phi, phi2) in degrees.

    Phi is drawn as arccos of a uniform variate so the sample is uniform in
    SO(3); drawing all three angles uniformly clusters at the poles.
    """
    rng = np.random.default_rng(seed)
    phi1 = rng.uniform(0.0, 360.0, n)
    Phi = np.degrees(np.arccos(rng.uniform(-1.0, 1.0, n)))
    phi2 = rng.uniform(0.0, 360.0, n)
    return np.column_stack([phi1, Phi, phi2])


def stage_orient(cfg: dict, deck: Path, out: Path) -> dict:
    z = np.load(out / "grains.npz")
    gid, xs, ys = z["gid"], z["xs"], z["ys"]
    n = int(gid.max()) + 1
    seed = int(cfg["orientations"]["seed"])
    euler = uniform_so3_euler(n, seed)
    thickness = float(cfg["mechanics"].get("thickness", 1.0))
    zmid = 0.5 * thickness
    with (out / "grains.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["grain_id", "phi1", "Phi", "phi2"])
        for g in range(n):
            w.writerow([g, f"{euler[g,0]:.6f}", f"{euler[g,1]:.6f}", f"{euler[g,2]:.6f}"])
    with (out / "orientations_points.csv").open("w", newline="") as fh:
        fh.write(f"# One row per growth-grid node; orientation seed {seed}; z at mid-thickness.\n")
        w = csv.writer(fh)
        w.writerow(["x", "y", "z", "phi1", "Phi", "phi2", "grain_id"])
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                g = int(gid[i, j])
                w.writerow([f"{x:.6f}", f"{y:.6f}", f"{zmid:.6f}",
                            f"{euler[g,0]:.6f}", f"{euler[g,1]:.6f}", f"{euler[g,2]:.6f}", g])
    rec = {"n_grains": n, "orientation_seed": seed, "point_file": str(out / "orientations_points.csv")}
    (out / "orient.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))
    return rec


# --------------------------------------------------------------------------- deck
def stage_deck(cfg: dict, deck: Path, out: Path) -> Path:
    mech = dict(cfg["mechanics"])
    src = TEMPLATE.read_text()
    sub = {"OUTDIR": str(out / "solver"), "POINTFILE": "orientations_points.csv"}
    for key, ph in MECH_KEYS.items():
        if key not in mech:
            raise SystemExit(f"[Mechanics] is missing '{key}'")
        v = mech[key]
        sub[ph] = f"{v:.12g}" if isinstance(v, float) else str(v)
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    leftover = sorted(set(re.findall(r"@([A-Z0-9_]+)@", re.sub(r"#.*", "", src))))
    if leftover:
        raise SystemExit(f"unsubstituted placeholders: {leftover}")
    path = out / "cpfe.i"
    header = (f"# Generated by scripts/polycrystal_plate.py from {deck.relative_to(STUDY) if deck.is_relative_to(STUDY) else deck}\n"
              f"# (stage 'deck'). Do not edit; edit the pipeline deck and regenerate.\n")
    path.write_text(header + src)
    print(f"wrote {path}")
    return path


# --------------------------------------------------------------------------- curve
def stage_curve(cfg: dict, deck: Path, out: Path) -> dict:
    """Nominal stress-strain from the CP executioner's per-step history."""
    hist = out / "solver" / "solver_history.csv"
    if not hist.exists():
        raise SystemExit(f"{hist} not found: run `pinnacle run {out/'cpfe.i'}` first")
    rows = list(csv.DictReader(hist.open()))
    mech = cfg["mechanics"]
    strain_max = float(mech["strain"])
    recs = []
    for r in rows:
        if r.get("converged", "True") not in ("True", "1", "1.0", "true"):
            continue
        lam = float(r["load_factor"])
        rec = {"step": int(float(r["step"])), "load_factor": lam,
               "nominal_strain_yy": lam * strain_max,
               "cutbacks": int(float(r.get("cutbacks", 0) or 0))}
        for k, v in r.items():
            if k.startswith("mean_stress") or k.startswith("mean_strain") or k == "time":
                rec[k] = float(v)
        recs.append(rec)
    if not recs:
        raise SystemExit("no converged steps in the history")
    cols = list(recs[0].keys())
    with (out / "results.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(recs)
    print(f"wrote {out/'results.csv'} ({len(recs)} steps; columns {cols})")
    return {"n_steps": len(recs), "columns": cols}


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    deck = Path(argv[0]).resolve()
    stages = argv[1:]
    if stages == ["all"]:
        stages = ["label", "orient", "deck"]
    cfg = read_pipeline_deck(deck)
    out = _resolve(str(cfg["mechanics"]["out_dir"]), deck)
    for st in stages:
        {"label": stage_label, "orient": stage_orient, "deck": stage_deck, "curve": stage_curve}[st](cfg, deck, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
