#!/usr/bin/env python3
"""Compare Gc calibrations at several plastic-work thresholds W_c (the W_c scan).

Deck-driven (PINNACLE ``.i`` syntax):

  [Figures]
    figures       = ../../../figures/pfd
    out_dir       = ../../../results/pfd_study/j2
    prefix        = j2_wc_scan_h40
    l0            = 5                         # the length scale compared across the roots
    markers       = ../../../data/digitized/al2219_fig3_experiment_markers.csv
    reference     = ../../../results/pfd_study/cp/neural_h40/base   # optional: a run drawn in grey (the CP+PFD plate)
    target_strain = 0.172
  [Roots]                                     # one calibration root (scripts/pfd_calibrate.py) per W_c
    [./wc0]    W_c = 0    dir = ../../../results/pfd_study/j2/cal_h40      []
    [./wc200]  W_c = 200  dir = ../../../results/pfd_study/j2/cal_h40_wc200 []

For every root the calibrated case of ``l0`` (``calibrated.json``) gives the
curve; the figure ``<prefix>_curves.png`` overlays them with the experiment
markers, and ``<prefix>_panels.png`` shows the calibrated Gc, the strain at
the peak stress (where the softening starts) and the peak stress against
W_c. ``<out_dir>/<prefix>_table.csv`` holds the numbers.

Usage: pfd_wc_scan_figures.py decks/pfd_study/j2/wc_scan_h40_figures.i
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
STYLE = {"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 7.5, "mathtext.fontset": "dejavusans",
         "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True}
COLOURS = ["#111111", "#1f4e9c", "#d1495b", "#00798c", "#edae49", "#6a4c93", "#2a9d8f"]


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("figures", "roots"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _curve(run: Path):
    rr = list(csv.DictReader((run / "results.csv").open()))
    e = np.array([float(r["nominal_strain"]) for r in rr])
    s = np.array([float(r["eng_stress"]) if r["eng_stress"] not in ("", "nan") else np.nan for r in rr])
    frac = json.loads((run / "fracture.json").read_text()) if (run / "fracture.json").exists() else {}
    return e, s, frac


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); cfg = read_deck(deck); f = cfg["figures"]
    figs = _resolve(str(f["figures"]), deck); figs.mkdir(parents=True, exist_ok=True)
    out_dir = _resolve(str(f.get("out_dir", f["figures"])), deck); out_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(f.get("prefix", "wc_scan")); l0 = float(f.get("l0", 5)); target = float(f.get("target_strain", 0.172))
    markers = []
    if f.get("markers"):
        for r in csv.DictReader((ln for ln in _resolve(str(f["markers"]), deck).open() if not ln.startswith("#"))):
            markers.append((float(r["strain_pct"]) / 100.0, float(r["stress_MPa"])))
    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=200)
    if f.get("reference"):
        ref = _resolve(str(f["reference"]), deck)
        if (ref / "results.csv").exists():
            e, s, _ = _curve(ref); ok = np.isfinite(s)
            ax.plot(np.concatenate([[0], 100 * e[ok]]), np.concatenate([[0], s[ok]]), color="#999999", lw=1.2, label="CP+PFD, calibrated (reference)")
    rows = []
    for k, (name, r) in enumerate(cfg["roots"].items()):
        root = _resolve(str(r["dir"]), deck); wc = float(r["W_c"])
        cal = root / "calibrated.json"
        if not cal.exists():
            print("skipping", name, "(no calibrated.json yet)"); continue
        entry = json.loads(cal.read_text())
        rec = next((v for kk, v in (entry.items() if isinstance(entry, dict) else [])
                    if abs(float(str(kk).lstrip("l")) - l0) < 1e-9), None) if isinstance(entry, dict) else None
        if rec is None and isinstance(entry, list):
            rec = next((v for v in entry if abs(float(v.get("l0", -1)) - l0) < 1e-9), None)
        if rec is None:
            print("skipping", name, "(no entry for l0", l0, ")"); continue
        gc = float(rec["Gc"]); run = root / f"l{l0:g}" / f"gc{gc:g}"
        if not (run / "results.csv").exists():   # the directory name carries the Gc as formatted by pfd_calibrate
            cands = sorted((root / f"l{l0:g}").glob("gc*"), key=lambda d: abs(float(d.name[2:]) - gc))
            run = cands[0] if cands else run
        e, s, frac = _curve(run); ok = np.isfinite(s)
        c = COLOURS[k % len(COLOURS)]
        label = f"W_c = {wc:g} MPa, Gc = {gc:.0f} MPa µm"
        ax.plot(np.concatenate([[0], 100 * e[ok]]), np.concatenate([[0], s[ok]]), color=c, lw=1.4, label=label)
        fs = frac.get("fracture_strain_node")
        if fs is not None:
            j = int(np.argmin(np.abs(e - fs))); ax.plot([100 * e[j]], [s[j]], "o", ms=5, color=c, mfc="white")
        rows.append({"W_c": wc, "Gc": gc, "fracture_strain": fs, "peak_eng_stress": frac.get("peak_eng_stress"),
                     "strain_at_peak": frac.get("strain_at_peak"), "fracture_energy_MPa": frac.get("fracture_energy_MPa"), "run": str(run)})
    if markers:
        ax.plot([100 * e for e, _ in markers], [s_ for _, s_ in markers], "x", ms=5, color="#111111", label="experiment markers (Vaughan 2024)")
    ax.axvline(100 * target, color="#888888", lw=0.8, ls=":")
    ax.set_xlim(0, 20.5); ax.set_ylim(0, 700); ax.grid(False)
    ax.set_xlabel(r"Nominal strain $\varepsilon$ (%)"); ax.set_ylabel(r"Engineering stress $\sigma$ (MPa), gross section")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(figs / f"{prefix}_curves.png"); plt.close(fig)
    if rows:
        rows.sort(key=lambda r: r["W_c"])
        with (out_dir / f"{prefix}_table.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), dpi=200)
        x = np.array([r["W_c"] for r in rows])
        for ax, key, lab in zip(axes, ("Gc", "strain_at_peak", "peak_eng_stress"), (r"calibrated $G_c$ (MPa µm)", r"strain at peak $\sigma$ (softening onset)", r"peak $\sigma$ (MPa)"), strict=True):
            y = np.array([r[key] if r[key] is not None else np.nan for r in rows], dtype=float)
            ax.plot(x, y, "o-", color="#1f4e9c", ms=5, lw=1.2)
            if key == "strain_at_peak":
                ax.axhline(target, color="#888888", lw=0.8, ls=":")
            ax.set_xlabel(r"$W_c$ (MPa)"); ax.set_ylabel(lab); ax.grid(False)
        fig.tight_layout(); fig.savefig(figs / f"{prefix}_panels.png"); plt.close(fig)
        for r in rows:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "run"})
    print("wrote", figs / f"{prefix}_curves.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
