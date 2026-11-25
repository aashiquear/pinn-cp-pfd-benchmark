#!/usr/bin/env python3
"""Part V, case K3: tensile targets from the raw NIST AM Bench 2025 records.

usage: scripts/lpbf_k3_targets.py <raw tensile data dir> <out.csv> [summary.json]

Reads `IN_718_MT*/specimen.dat` of (benzing 2025 tensiledata): time, cross-head
displacement, extensometer displacement (mm) and force (N), with the specimen width and
thickness in the file header. Engineering strain = extensometer displacement / 3 mm (the
gauge length of the extensometer of these miniature specimens), engineering stress =
force / (width x thickness).

Writes the mean engineering curve of all specimens on a common strain grid up to the
smallest strain at the ultimate stress (uniform elongation; later the specimens neck and
the curve is not a material curve), with the standard deviation between specimens, the
true stress and the logarithmic plastic strain; and per specimen: Young's modulus (fit
between 100 and 400 MPa), 0.2 % offset yield stress, ultimate stress and its strain.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

GAUGE_MM = 3.0


def read(path: Path) -> dict:
    width = thickness = None
    rows = []
    started = False
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split("\t")
        if parts[0].startswith("Specimen Width"):
            width = float(parts[1])
        elif parts[0].startswith("Specimen Thickness"):
            thickness = float(parts[1])
        elif parts[0] == "Sec":
            started = True
        elif started and len(parts) >= 4:
            try:
                rows.append([float(v) for v in parts[:4]])
            except ValueError:
                pass
    data = np.array(rows)
    strain = data[:, 2] / GAUGE_MM
    stress = data[:, 3] / (width * thickness)
    return {"strain": strain, "stress": stress, "width_mm": width, "thickness_mm": thickness}


def properties(strain: np.ndarray, stress: np.ndarray) -> dict:
    peak = int(np.argmax(stress))
    elastic = (stress > 100.0) & (stress < 400.0) & (np.arange(stress.size) < peak)
    slope, intercept = np.polyfit(strain[elastic], stress[elastic], 1)
    offset = stress[:peak] - slope * (strain[:peak] - 0.002) - intercept
    k = int(np.argmax(offset < 0.0))
    f = offset[k - 1] / (offset[k - 1] - offset[k])
    return {
        "E_GPa": slope / 1e3,
        "yield_0p2_MPa": float(stress[k - 1] + f * (stress[k] - stress[k - 1])),
        "ultimate_MPa": float(stress[peak]),
        "strain_at_ultimate": float(strain[peak]),
        "toe_strain": float(-intercept / slope),
    }


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    root = Path(sys.argv[1])
    specimens = {}
    for path in sorted(root.glob("IN_718_MT*/specimen.dat")):
        rec = read(path)
        rec.update(properties(rec["strain"], rec["stress"]))
        specimens[path.parent.name] = rec
    limit = min(r["strain_at_ultimate"] - r["toe_strain"] for r in specimens.values())
    grid = np.concatenate([np.linspace(0.0, 0.01, 41)[:-1], np.linspace(0.01, limit, 97)])
    curves = []
    for rec in specimens.values():
        peak = int(np.argmax(rec["stress"]))
        e = rec["strain"][: peak + 1] - rec["toe_strain"]  # remove the toe of the record
        s = rec["stress"][: peak + 1]
        keep = e >= 0.0
        order = np.argsort(e[keep])
        curves.append(np.interp(grid, np.concatenate([[0.0], e[keep][order]]),
                                np.concatenate([[0.0], s[keep][order]])))
    curves = np.array(curves)
    mean, sd = curves.mean(axis=0), curves.std(axis=0, ddof=1)
    modulus = float(np.mean([r["E_GPa"] for r in specimens.values()])) * 1e3
    true_stress = mean * (1.0 + grid)
    plastic = np.log(1.0 + grid) - true_stress / modulus
    out = Path(sys.argv[2])
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("# mean of %d specimens of (benzing 2025 tensiledata), build direction, 1e-3 /s; "
                 "derived by scripts/lpbf_k3_targets.py\n" % len(specimens))
        fh.write("eng_strain,eng_stress_MPa,sd_MPa,true_stress_MPa,log_plastic_strain\n")
        for row in zip(grid, mean, sd, true_stress, plastic):
            fh.write("%.6f,%.3f,%.3f,%.3f,%.6f\n" % row)
    keys = ("E_GPa", "yield_0p2_MPa", "ultimate_MPa", "strain_at_ultimate")
    summary = {
        "specimens": {n: {k: round(r[k], 4) for k in keys} | {"width_mm": r["width_mm"], "thickness_mm": r["thickness_mm"]}
                      for n, r in specimens.items()},
        "mean": {k: round(float(np.mean([r[k] for r in specimens.values()])), 4) for k in keys},
        "sd": {k: round(float(np.std([r[k] for r in specimens.values()], ddof=1)), 4) for k in keys},
        "curve_strain_limit": round(float(limit), 5),
        "gauge_mm": GAUGE_MM,
    }
    if len(sys.argv) > 3:
        Path(sys.argv[3]).write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in ("mean", "sd", "curve_strain_limit")}, indent=1))


if __name__ == "__main__":
    main()
