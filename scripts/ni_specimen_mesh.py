#!/usr/bin/env python3
"""Dog-bone meshes of the Ni single-crystal specimens (hexaform ``dogbone``).

Deck-driven (PINNACLE ``.i`` syntax):

  [Mesh]
    out_root       = ../../../results/ni_study/mesh/h020     # one directory per specimen below it
    planform       = ../../../data/extracted/ni_fig2_planform.csv
    thickness      = ../../../data/extracted/ni_table2_thickness.csv
    specimens      = all                      # or a list: 'L.1 9.1'
    station_spacing = 2.0                     # mm between t1..t5 (Fig. 2)
    grip_thickness = 0.35                     # the common thickness outside the gage (author 2026-10-15)
    h              = 0.2                      # gage cell size (mm)
    h_grip         = 0.4
    n_z            = 1
    figure         = ../../../figures/ni/mesh_h020.png

For every specimen of Table 2 the thickness stations become the gage
profile, the planform (gage length/width, grip width, fillet radius, total
length) comes from the Fig. 2 file, and hexaform writes ``mesh.e``,
``mesh.vtk`` and ``mesh.json`` under ``<out_root>/<specimen>/``; a figure
shows the planform with the mesh of one specimen and the fifteen thickness
profiles. Needs hexaform on the path (HEXAFORM_SRC, default
/mnt/extended/projects/hexaform/src).

Usage: ni_specimen_mesh.py decks/ni/mesh/h020.i
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, os.environ.get("HEXAFORM_SRC", "/mnt/extended/projects/hexaform/src"))


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "mesh" not in cfg:
        raise SystemExit(f"{path}: missing [Mesh] block")
    return cfg["mesh"]


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); m = read_deck(deck)
    from hexaform.core.validation import validate
    from hexaform.dogbone import DogboneOptions, dogbone_mesh, dogbone_summary, half_width, thickness
    from hexaform.io import exodus, vtk

    plan = {r["quantity"]: float(r["value_mm"]) for r in _rows(_resolve(str(m["planform"]), deck))}
    thick = _rows(_resolve(str(m["thickness"]), deck))
    want = m.get("specimens", "all")
    want = None if str(want) == "all" else set(str(want).split() if not isinstance(want, (list, tuple)) else [str(x) for x in want])
    spacing = float(m.get("station_spacing", 2.0)); h = float(m["h"]); h_grip = float(m.get("h_grip", 2 * h)); n_z = int(m.get("n_z", 1))
    t_grip = float(m.get("grip_thickness", plan.get("nominal_thickness", 0.35)))
    out_root = _resolve(str(m["out_root"]), deck); out_root.mkdir(parents=True, exist_ok=True)
    summary = {}
    profiles = []
    for r in thick:
        name = r["specimen"]
        if want and name not in want:
            continue
        stations = [(k * spacing, float(r[f"t{k + 1}_mm"])) for k in range(5)]
        opts = DogboneOptions(length=plan["total_length"], gage_length=plan["gage_length"], gage_width=plan["gage_width"],
                              grip_width=plan["grip_width"], fillet_radius=plan["fillet_radius"], grip_thickness=t_grip,
                              stations=stations, h=h, h_grip=h_grip, n_z=n_z)
        mesh = dogbone_mesh(opts)
        report = validate(mesh, cell_size=h)
        d = out_root / name; d.mkdir(parents=True, exist_ok=True)
        exodus.write(mesh, d / "mesh.e", title=f"Ni specimen {name}")
        vtk.write(mesh, d / "mesh.vtk")
        s = dogbone_summary(mesh); s["validation"] = report.as_dict(); s["specimen"] = name
        (d / "mesh.json").write_text(json.dumps(s, indent=2, default=str))
        summary[name] = {"elements": s["n_elements"], "nodes": s["n_nodes"], "gage_cells": s["gage_cells"],
                         "thickness_range": s["thickness_range"], "valid": bool(report.passed),
                         "min_scaled_jacobian": s["validation"].get("min_scaled_jacobian", s["validation"].get("scaled_jacobian_min"))}
        xs = np.linspace(0, plan["total_length"], 400)
        profiles.append((name, r["group"], xs, thickness(opts, xs), mesh if name == (next(iter(want)) if want else "L.1") else None, opts))
        print(name, summary[name])
    (out_root / "summary.json").write_text(json.dumps(summary, indent=2))
    if m.get("figure") and profiles:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection

        plt.rcParams.update({"font.size": 8, "axes.labelsize": 9})
        fig, axes = plt.subplots(3, 1, figsize=(8.5, 7.5), dpi=200, gridspec_kw={"height_ratios": [1.3, 1.2, 1.4]})
        name, group, xs, _, mesh, opts = next((p for p in profiles if p[4] is not None), profiles[0])
        if mesh is None:
            mesh = dogbone_mesh(opts)
        ax = axes[0]
        pts = mesh.points; cells = mesh.cells; bottom = cells[:, :4]
        ax.add_collection(PolyCollection(pts[bottom][:, :, :2], facecolors=["#d9e3f0", "#f5deb3", "#dddddd"][0], edgecolors="#555555", linewidths=0.25))
        colors = {1: "#d9e3f0", 2: "#f5deb3", 3: "#dddddd"}
        for b, c in colors.items():
            sel = mesh.block_id == b
            if sel.any():
                ax.add_collection(PolyCollection(pts[bottom[sel]][:, :, :2], facecolors=c, edgecolors="#555555", linewidths=0.25))
        x0 = 0.5 * (opts.length - opts.gage_length)
        for k in range(5):
            ax.axvline(x0 + k * spacing, color="#d1495b", lw=0.6, ls=":")
        ax.set_aspect("equal"); ax.set_xlim(0, opts.length); ax.set_ylim(-0.55 * opts.grip_width, 0.55 * opts.grip_width)
        ax.set_title(f"specimen {name}: {mesh.n_elements} HEX8, {mesh.n_nodes} nodes, gage cell {h} mm, grips {h_grip} mm, {n_z} through the thickness; blocks gage / shoulder / grip, stations dotted", fontsize=7.5)
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        ax = axes[1]
        xg = np.linspace(x0 - opts.fillet_extent - 0.5, x0 + opts.gage_length + opts.fillet_extent + 0.5, 300)
        yw = half_width(opts, xg); tz = thickness(opts, xg)
        ax.fill_between(xg, -yw, yw, color="#d9e3f0", lw=0)
        ax2 = ax.twinx(); ax2.plot(xg, tz, color="#1f4e9c", lw=1.2); ax2.set_ylabel("thickness (mm)", color="#1f4e9c")
        for k in range(5):
            ax.axvline(x0 + k * spacing, color="#d1495b", lw=0.6, ls=":")
        ax.set_xlim(xg[0], xg[-1]); ax.set_ylim(-0.55 * opts.grip_width, 0.55 * opts.grip_width); ax.set_aspect("auto")
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_title(f"planform around the gage (fillet R {opts.fillet_radius} mm, gage {opts.gage_width} x {opts.gage_length} mm) and the thickness profile of {name}", fontsize=7.5)
        ax = axes[2]
        cmap = {"9": "#1f4e9c", "10": "#2a9d8f", "H": "#6a4c93", "L": "#d1495b", "P": "#111111"}
        styles = {"x.1": "-", "x.3": "--", "x.4": ":"}
        for name_i, group_i, xs_i, t_i, _, opts_i in profiles:
            variant = ("x." + name_i.split(".")[1]) if "." in name_i else name_i   # synthetic specimens (S01...) have no variant
            ax.plot(xs_i, t_i, color=cmap.get(group_i, "#888888"), ls=styles.get(variant, "-"), lw=1.0, label=name_i)
        ax.set_xlim(0, opts.length); ax.set_xlabel("x (mm)"); ax.set_ylabel("thickness (mm)")
        ax.set_title(f"thickness profiles of the {len(profiles)} specimens (Table 2 stations, spline in the gage, linear blend to {t_grip} mm across the fillets)", fontsize=7.5)
        ax.legend(ncol=5, fontsize=6, frameon=False)
        fig.tight_layout(); fp = _resolve(str(m["figure"]), deck); fp.parent.mkdir(parents=True, exist_ok=True); fig.savefig(fp); plt.close(fig)
        print("wrote", fp)
    print("wrote", out_root / "summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
