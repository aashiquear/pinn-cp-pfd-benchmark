#!/usr/bin/env python3
"""Ship the result files behind the report: force-add them to git and write the manifest.

Deck-driven (PINNACLE ``.i`` syntax):

  [Ship]
    manifest = ../../results/SHIPPED.md
  [Series]                       # field series: which snapshots of solver/fields_step*.vtk
    [./cal_cp]  dir = ../../results/pfd_study/cp/neural_h40/base  frames = all  []
    [./shape]   dir = ../../results/pfd_study/cp/shape_h40/*       frames = every:2  []
    [./mesh80]  dir = ../../results/cp_study/mesh/h80              frames = '2 8 15 20'  []
  [Files]                        # other ignored files (meshes, networks), globs
    [./meshes]  glob = ../../results/pfd_study/mesh/*/mesh.e  []

``frames`` is ``all``, ``every:k`` (every k-th snapshot and the last) or a
list of nominal strains in % (the nearest snapshot of each, as the figure
scripts pick them, plus the last); ``fields.vtk`` and ``fields.pvd`` of the
run are always included. The snapshot-to-strain map follows
scripts/pfd_fields.py: snapshots are written every n_steps / n_snapshots
steps, the strain comes from ``results.csv``. Every file is ``git add -f``-ed
(the .gitignore keeps field files out otherwise) and listed in the manifest
with its size. Usage: ship_results.py decks/report/ship_results.i
"""
from __future__ import annotations

import csv
import glob
import subprocess
import sys
from pathlib import Path

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "ship" not in cfg:
        raise SystemExit(f"{path}: missing [Ship] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def _series(run: Path) -> list[tuple[float, Path]]:
    files = sorted((run / "solver").glob("fields_step*.vtk"))
    res = run / "results.csv"
    if not files:
        return []
    if not res.exists():
        return [(float("nan"), f) for f in files]
    rows = list(csv.DictReader(res.open()))
    key = "nominal_strain" if "nominal_strain" in rows[0] else "strain_pct"
    strains = [float(r[key]) * (1.0 if key == "nominal_strain" else 0.01) for r in rows]
    k = max(1, round(len(strains) / len(files)))
    return [(strains[min((i + 1) * k - 1, len(strains) - 1)], f) for i, f in enumerate(files)]


def select(run: Path, frames: str) -> list[Path]:
    snaps = _series(run)
    out: list[Path] = []
    if frames == "all":
        out = [f for _, f in snaps]
    elif frames.startswith("every:"):
        k = int(frames.split(":", 1)[1])
        out = [f for i, (_, f) in enumerate(snaps) if i % k == 0]
    else:
        targets = [float(x) / 100.0 for x in frames.split()]
        for t in targets:
            if snaps:
                out.append(min(snaps, key=lambda s: abs(s[0] - t))[1])
    if snaps and snaps[-1][1] not in out:
        out.append(snaps[-1][1])
    for extra in ("fields.vtk", "fields.pvd"):
        if (run / "solver" / extra).exists():
            out.append(run / "solver" / extra)
    return sorted(set(out))


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); cfg = read_deck(deck)
    manifest = _resolve(str(cfg["ship"].get("manifest", "../../results/SHIPPED.md")), deck)
    shipped: list[tuple[str, Path, list[Path]]] = []
    for name, sp in cfg.get("series", {}).items():
        frames = _str(sp.get("frames", "all"))
        for run in sorted(glob.glob(str(_resolve(str(sp["dir"]), deck)))):
            run = Path(run)
            if not (run / "solver").is_dir():
                continue
            files = select(run, frames)
            if files:
                shipped.append((f"{name}: {frames}", run, files))
    for name, fp in cfg.get("files", {}).items():
        files = [Path(p) for p in sorted(glob.glob(str(_resolve(str(fp["glob"]), deck))))]
        if files:
            shipped.append((name, STUDY_ROOT, files))
    total = 0
    lines = ["# Shipped result files", "",
             "Field series, meshes and networks behind the report's figures and tables,",
             "force-added to git by `scripts/ship_results.py decks/report/ship_results.i`",
             "(the .gitignore keeps every other field file out). Series list the",
             "snapshots the figures read; `all` series feed an animation.", "",
             "| set | run | files | size (MB) |", "|---|---|---|---|"]
    for label, run, files in shipped:
        size = sum(f.stat().st_size for f in files); total += size
        rel = run.relative_to(STUDY_ROOT) if run != STUDY_ROOT else Path(".")
        names = ", ".join(f.name for f in files) if len(files) <= 6 else f"{len(files)} files ({files[0].name} .. {files[-1].name})"
        lines.append(f"| {label} | `{rel}` | {names} | {size / 1e6:.1f} |")
        subprocess.run(["git", "add", "-f", *[str(f) for f in files]], cwd=STUDY_ROOT, check=True)
    lines += ["", f"Total: {total / 1e6:.0f} MB in {sum(len(f) for _, _, f in shipped)} files."]
    manifest.write_text("\n".join(lines) + "\n")
    subprocess.run(["git", "add", str(manifest)], cwd=STUDY_ROOT, check=True)
    print(f"shipped {sum(len(f) for _, _, f in shipped)} files, {total / 1e6:.0f} MB; manifest {manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
