#!/usr/bin/env python3
"""Design of the Ni specimen surrogate's training set: synthetic single-crystal specimens.

Deck-driven (PINNACLE .i syntax):

  [Design]
    out_root   = ../../../results/ni_study/final/surrogate/data
    samples    = 30
    seed       = 0
    t_mean     = '0.325 0.400'      # mm, range of the gage mean thickness (Table 2: 0.33-0.39)
    t_taper    = '-0.04 0.06'       # mm, t5 - t1 over the 8 mm gage (Table 2: -0.02 ... +0.059)
    mesh_deck  = ../mesh/h020.i     # mesh settings copied (planform, h, grips); thickness/out_root replaced
    study_base = ../final/fem_base_h020.i   # base specimen deck of the FE runs (carried set)
    concurrent = 6
  []

Latin-hypercube sample (author 2026-10-25: 30 samples) over five variables: the three Bunge
Euler angles, uniform over the rotation group (phi1 and phi2 uniform on [0, 360), Phi as
arccos of a uniform variable on [-1, 1]), and the gage thickness as a mean and a linear
taper over the five stations (t_i = mean + taper (i - 3) / 4). Writes under out_root:
design.csv (one row per sample), orientations.csv and thickness.csv in the formats of the
Table 1 / Table 2 files (group = specimen = S01..S30, so that scripts/ni_specimen.py reads
them as it reads the real specimens), mesh.i for scripts/ni_specimen_mesh.py, and study.i
for scripts/ni_study.py (one case per sample, meshes under <out_root>/mesh). The 15 real
specimens are not in the design: they are the held-out test of the surrogate.

Usage: ni_surrogate_design.py DESIGN.i
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_figures import _resolve, _str  # noqa: E402


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "design" not in cfg:
        raise SystemExit(f"{path}: missing [Design]")
    return cfg["design"]


def _pair(v) -> tuple[float, float]:
    a, b = (float(x) for x in (v if isinstance(v, (list, tuple)) else str(v).split()))
    return a, b


def lhs(n: int, dim: int, rng: np.random.Generator) -> np.ndarray:
    """Latin hypercube on [0, 1)^dim: one point per stratum per dimension, strata shuffled."""
    u = np.empty((n, dim))
    for j in range(dim):
        u[:, j] = (rng.permutation(n) + rng.random(n)) / n
    return u


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    d = read_deck(deck)
    out = _resolve(_str(d["out_root"]), deck)
    out.mkdir(parents=True, exist_ok=True)
    n = int(d.get("samples", 30))
    rng = np.random.default_rng(int(d.get("seed", 0)))
    u = lhs(n, 5, rng)
    tm0, tm1 = _pair(d.get("t_mean", "0.325 0.400"))
    tp0, tp1 = _pair(d.get("t_taper", "-0.04 0.06"))
    rows = []
    for i in range(n):
        phi1 = 360.0 * u[i, 0]
        Phi = float(np.degrees(np.arccos(1.0 - 2.0 * u[i, 1])))
        phi2 = 360.0 * u[i, 2]
        tmean = tm0 + (tm1 - tm0) * u[i, 3]
        taper = tp0 + (tp1 - tp0) * u[i, 4]
        t = [tmean + taper * (k - 3) / 4.0 for k in range(1, 6)]
        rows.append({"specimen": f"S{i + 1:02d}", "phi1_deg": round(phi1, 2), "Phi_deg": round(Phi, 2), "phi2_deg": round(phi2, 2),
                     "t_mean_mm": round(tmean, 4), "t_taper_mm": round(taper, 4), **{f"t{k + 1}_mm": round(t[k], 4) for k in range(5)}})
    with (out / "design.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with (out / "orientations.csv").open("w", newline="") as fh:
        fh.write("# Synthetic single-crystal orientations of the surrogate design (scripts/ni_surrogate_design.py),\n# Bunge Euler angles (deg), uniform over the rotation group; format of data/extracted/ni_table1_orientations.csv.\n")
        w = csv.DictWriter(fh, fieldnames=["group", "phi1_deg", "Phi_deg", "phi2_deg", "max_schmid"]); w.writeheader()
        w.writerows({"group": r["specimen"], "phi1_deg": r["phi1_deg"], "Phi_deg": r["Phi_deg"], "phi2_deg": r["phi2_deg"], "max_schmid": ""} for r in rows)
    with (out / "thickness.csv").open("w", newline="") as fh:
        fh.write("# Synthetic gage thickness profiles of the surrogate design (mean + linear taper over the five\n# stations); format of data/extracted/ni_table2_thickness.csv.\n")
        w = csv.DictWriter(fh, fieldnames=["specimen", "group", "variant", "t1_mm", "t2_mm", "t3_mm", "t4_mm", "t5_mm"]); w.writeheader()
        w.writerows({"specimen": r["specimen"], "group": r["specimen"], "variant": "lhs", **{f"t{k}_mm": r[f"t{k}_mm"] for k in range(1, 6)}} for r in rows)
    # mesh deck: the real mesh deck with the thickness table and the output root replaced; every
    # path relative to the written deck, so the design travels to the compute host unchanged
    def rel(path: Path) -> str:
        r = os.path.relpath(path, out)
        return r if r.startswith(".") else "./" + r   # scripts/ni_study.py re-bases "./" and "../" values

    mesh_src = _resolve(_str(d["mesh_deck"]), deck)
    mt = mesh_src.read_text().split("\n")
    mesh_dir = out / "mesh"
    keep = []
    for ln in mt:
        if ln.lstrip().startswith("#"):
            continue
        k = ln.split("=")[0].strip()
        if k == "out_root":
            ln = f"  out_root        = {rel(mesh_dir)}"
        elif k == "thickness":
            ln = f"  thickness       = {rel(out / 'thickness.csv')}"
        elif k == "figure":
            ln = f"  figure          = {rel(out / 'mesh.png')}"
        elif k in ("planform",):
            ln = f"  {k:<15} = {rel(_resolve(ln.split('=', 1)[1].split('#')[0].strip(), mesh_src))}"
        keep.append(ln)
    (out / "mesh.i").write_text("# Meshes of the surrogate design's synthetic specimens (written by scripts/ni_surrogate_design.py).\n" + "\n".join(keep).strip() + "\n")
    # study deck for the FE runs
    base = _resolve(_str(d["study_base"]), deck)
    cases = "".join(f"  [./{r['specimen']}]\n    group        = {r['specimen']}\n    specimen     = {r['specimen']}\n    orientations = {rel(out / 'orientations.csv')}\n"
                    f"    thickness    = {rel(out / 'thickness.csv')}\n    mesh         = {rel(mesh_dir / r['specimen'] / 'mesh.e')}\n  []\n" for r in rows)
    (out / "study.i").write_text(f"""# FE runs of the surrogate design's {n} synthetic specimens (written by scripts/ni_surrogate_design.py;
# author 2026-10-25): the carried CP+PFD set of the base deck, one case per sample. Driver scripts/ni_study.py.
[Study]
  base         = {rel(base)}
  out_root     = {rel(out / 'runs')}
  concurrent   = {int(d.get('concurrent', 6))}
  include_base = false
[]
[Cases]
{cases}[]
""")
    print(f"wrote {n} samples under {out}: design.csv, orientations.csv, thickness.csv, mesh.i, study.i")
    t = np.array([[r[f"t{k}_mm"] for k in range(1, 6)] for r in rows])
    print("thickness range %.3f-%.3f mm; Phi range %.0f-%.0f deg" % (t.min(), t.max(), min(r["Phi_deg"] for r in rows), max(r["Phi_deg"] for r in rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
