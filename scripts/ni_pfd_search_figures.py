#!/usr/bin/env python3
"""The G_c searches of the Ni CP+PFD calibration as one figure (docs/12 §6.2), driven by one deck.

  [Figures]
    out            = ../../../figures/ni/pfd/search_L1_cal.png
    root           = ../../../results/ni_study/pfd/cal_L1_h020   # <root>/l<l0>/gc<Gc>/ runs (ni_pfd_calibrate.py)
    failure_target = 0.463
    tolerance      = 0.004
    h              = 0.2                                         # gage cell (mm), for l0/h in the legend
    dpi            = 200
  []

Left: failure gauge strain against G_c, one series per length scale (open square: a run
that stopped without separating, drawn at its last strain). Right: the same runs against
G_c / l0, which collapses the length scales. The target band is the calibration tolerance.
Void runs (directories named void_*) are skipped.

Usage: ni_pfd_search_figures.py FIGURES.i
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_figures import _resolve, _str  # noqa: E402

def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "figures" not in cfg:
        raise SystemExit(f"{path}: missing [Figures] block")
    return cfg


COLORS = {"0.4": "#d1495b", "0.6": "#edae49", "0.8": "#00798c", "1.2": "#1f4e9c"}


def runs_of(root: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for ldir in sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("l")):
        l0 = ldir.name[1:]
        for cdir in sorted(p for p in ldir.iterdir() if p.is_dir() and p.name.startswith("gc")):
            if not (cdir / "summary.json").exists():
                continue
            s = json.loads((cdir / "summary.json").read_text())
            gc = float(json.loads((cdir / "case.json").read_text())["substitutions"]["GC"])
            fs = s.get("failure_gauge_strain")
            out.setdefault(l0, []).append(
                {"Gc": gc, "fail": fs, "end": s.get("final_gauge_strain"), "l0": float(l0)}
            )
    return out


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    f = read_deck(deck)["figures"]
    root = _resolve(_str(f["root"]), deck)
    target, tol = float(f["failure_target"]), float(f.get("tolerance", 0.004))
    h = float(f.get("h", 0.2))
    data = runs_of(root)
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.4, 3.4), dpi=int(f.get("dpi", 200)))
    for ax in (a, b):
        ax.axhspan(100 * (target - tol), 100 * (target + tol), color="#dddddd", lw=0)
        ax.axhline(100 * target, color="#111111", lw=0.9, ls=":")
        ax.set_xscale("log")
        ax.set_ylim(0, 55)
        ax.set_ylabel("failure gauge strain (%)")
        for side in ("top", "right"):
            ax.spines[side].set_visible(True)
    for l0, rs in sorted(data.items(), key=lambda kv: float(kv[0])):
        c = COLORS.get(l0, "#9a9a9a")
        rs = sorted(rs, key=lambda r: r["Gc"])
        hit = [r for r in rs if r["fail"] is not None]
        miss = [r for r in rs if r["fail"] is None and r["end"] is not None]
        lab = r"$\ell_0$ = %s mm ($\ell_0/h$ = %g)" % (l0, float(l0) / h)
        for ax, xk in ((a, 1.0), (b, 1.0 / float(l0))):
            ax.plot([r["Gc"] * xk for r in hit], [100 * r["fail"] for r in hit], "o-", color=c,
                    ms=4.5, lw=1.0, label=lab if ax is a else None)
            ax.plot([r["Gc"] * xk for r in miss], [100 * r["end"] for r in miss], "s", color=c,
                    ms=5, mfc="none", mew=1.3)
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter

    for ax, ticks in ((a, [5, 10, 20, 40, 80]), (b, [10, 20, 30, 50, 80])):
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r"$G_c$ (MPa mm)")
    b.set_xlabel(r"$G_c / \ell_0$ (MPa)")
    a.legend(loc="lower right", frameon=False, fontsize=7)
    a.text(0.03, 0.95, "open square: stopped without separating", transform=a.transAxes,
           fontsize=6.5, va="top")
    fig.tight_layout()
    out = _resolve(_str(f["out"]), deck)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
