#!/usr/bin/env python3
"""Per-specimen figures of a last-stage route: curve + max d, and E_xx at failure with the crack line.

Deck-driven (PINNACLE .i syntax):

  [Figures]
    root     = ../../../results/ni_study/final/fem     # <specimen>/summary.json per case (lower-case names)
    out_dir  = ../../../figures/ni/final
    prefix   = fem
    label    = 'FEM'
    parity   = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv   # experiment column only
    scales   = ../../../data/digitized/ni_fig5_panel_scales.csv
    panels   = ../../../data/digitized/ni_fig5_panels
    cp_only  = ../../../results/ni_study/cal_h020/gn2/base/results.csv
    specimens = all                                      # or a list: '10.1 10.3'
    dpi      = 200
  []

For every finished specimen it writes a curve deck and a shape deck under <out_dir>/decks/ and
runs scripts/ni_pfd_curve_figures.py and scripts/ni_pfd_shape_figures.py on them:
<out_dir>/<prefix>_curve_<spec>.png (stress-strain and max d, the experiment's failure strain
of that specimen as the dotted line; the Fig. 3 trace only for L.1) and
<out_dir>/<prefix>_shape_<spec>.png (E_xx on the experiment panel's own color scale, the d =
0.85 outline and the crack line, against the DIC panel). Experiment data only.

Usage: ni_final_specimen_figures.py FIGURES.i
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ni_final_figures import read_deck as _read  # noqa: E402
from report_figures import _resolve, _str  # noqa: E402

HERE = Path(__file__).resolve().parent


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "figures" not in cfg:
        raise SystemExit(f"{path}: missing [Figures]")
    return cfg["figures"]


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    f = read_deck(deck)
    root = _resolve(_str(f["root"]), deck)
    out_dir = _resolve(_str(f["out_dir"]), deck)
    (out_dir / "decks").mkdir(parents=True, exist_ok=True)
    prefix, label = _str(f.get("prefix", "fem")), _str(f.get("label", "FEM"))
    exp = {r["specimen"]: float(r["exp_failure_strain"]) for r in _rows(_resolve(_str(f["parity"]), deck))}
    scales = {r["specimen"]: (float(r["exx_min"]), float(r["exx_max"])) for r in _rows(_resolve(_str(f["scales"]), deck))}
    panels = _resolve(_str(f["panels"]), deck)
    want = _str(f.get("specimens", "all"))
    specs = sorted(exp) if want == "all" else want.split()
    done = 0
    for spec in specs:
        case = root / spec.lower()
        if not (case / "summary.json").exists():
            continue
        lo, hi = scales.get(spec, (0.23, 0.98))
        hide = "false" if spec == "L.1" else "true"
        cdeck = out_dir / "decks" / f"{prefix}_curve_{spec}.i"
        cdeck.write_text(f"""# {label} route, specimen {spec} (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out            = ../{prefix}_curve_{spec}.png
  markers        = {_resolve(_str(f.get('markers', '../../../data/digitized/ni_fig3_experiment_markers.csv')), deck)}
  experiment     = {_resolve(_str(f.get('experiment', '../../../data/digitized/ni_fig3_specimenL_calibration.csv')), deck)}
  cp_only        = {_resolve(_str(f['cp_only']), deck)}
  failure_target = {exp[spec]}
  hide_experiment = {hide}
  damage         = true
  strain_max     = 80
  stress_max     = 300
  stress_tick    = 50
  dpi            = {int(f.get('dpi', 200))}
[]
[Runs]
  [./run]
    dir   = {case}
    label = 'specimen {spec}, {label}'
    color = '#1f4e9c'
  []
[]
""")
        sdeck = out_dir / "decks" / f"{prefix}_shape_{spec}.i"
        sdeck.write_text(f"""# {label} route, specimen {spec} (written by scripts/ni_final_specimen_figures.py).
[Figures]
  out         = ../{prefix}_shape_{spec}.png
  panel       = {panels / (spec + '_exp.png')}
  panel_label = '{spec}: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '{lo} {hi}'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'
  color       = exx
  dpi         = {int(f.get('dpi', 200))}
[]
[Runs]
  [./run]
    dir   = {case}
    label = 'specimen {spec}, {label}'
  []
[]
""")
        for script, d in (("ni_pfd_curve_figures.py", cdeck), ("ni_pfd_shape_figures.py", sdeck)):
            subprocess.run([sys.executable, str(HERE / script), str(d)], check=True, stdout=subprocess.DEVNULL)
        done += 1
        print("drew", spec)
    print(f"{done} specimens drawn under {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
