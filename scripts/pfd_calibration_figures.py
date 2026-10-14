#!/usr/bin/env python3
"""Figures and table of a Gc(l0) calibration (scripts/pfd_calibrate.py output).

Deck-driven (PINNACLE ``.i`` syntax):

  [Figures]
    cal_root  = ../../../results/pfd_study/j2/cal_h40     # calibration.csv, calibrated.json, l*/gc*/results.csv
    out_dir   = ../../../results/pfd_study/j2/cal_h40
    figures   = ../../../figures/pfd
    prefix    = j2_cal_h40
    markers   = ../../../data/digitized/al2219_fig3_experiment_markers.csv
    target_strain = 0.172
    damage_threshold = 0.8

Writes ``<out_dir>/<prefix>_table.csv`` (one row per l0: calibrated Gc,
fracture strain, peak stress, fracture energy, steps, Newton iterations,
wall time) and three figures: the engineering stress-strain curve of every
calibrated case against the experiment markers (the fracture point marked),
the fracture strain of every trial against Gc per l0 (the bisection
record), and Gc against l0 (log-log, with the Gc l0 product).

Usage: pfd_calibration_figures.py decks/pfd_study/j2/cal_h40_figures.i
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

STYLE = {"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "mathtext.fontset": "dejavusans",
         "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True}


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "figures" not in cfg:
        raise SystemExit(f"{path}: missing [Figures] block")
    return cfg["figures"]


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _curve(d: Path):
    rows = list(csv.DictReader((d / "results.csv").open()))
    e = np.array([float(r["nominal_strain"]) for r in rows])
    s = np.array([float(r["eng_stress"]) if r["eng_stress"] not in ("", "nan") else np.nan for r in rows])
    dmax = np.array([float(r["max_damage"]) for r in rows])
    return e, s, dmax


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    f = read_deck(deck)
    root = _resolve(str(f["cal_root"]), deck)
    out = _resolve(str(f.get("out_dir", f["cal_root"])), deck)
    figs = _resolve(str(f.get("figures", f["cal_root"])), deck)
    figs.mkdir(parents=True, exist_ok=True)
    prefix = str(f.get("prefix", "cal"))
    target = float(f.get("target_strain", 0.172))
    log = list(csv.DictReader((root / "calibration.csv").open()))
    calibrated = json.loads((root / "calibrated.json").read_text())
    markers = []
    if f.get("markers"):
        for r in csv.DictReader((ln for ln in _resolve(str(f["markers"]), deck).open() if not ln.startswith("#"))):
            markers.append((float(r.get("strain_pct", 0)) / 100.0, float(r.get("stress_MPa", 0))))
    plt.rcParams.update(STYLE)
    colors = ["#1f4e9c", "#d1495b", "#00798c", "#edae49", "#6a4c93"]

    # ---- table and curves of the calibrated cases ----
    table = []
    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=200)
    for i, c in enumerate(sorted(calibrated, key=lambda r: float(r["l0"]))):
        if c.get("Gc") is None:
            continue
        d = Path(c["dir"])
        if not d.exists():
            d = root / f"l{float(c['l0']):g}" / f"gc{float(c['Gc']):.6g}"
        frac = json.loads((d / "fracture.json").read_text())
        e, s, dmax = _curve(d)
        ok = np.isfinite(s)
        ax.plot(np.concatenate([[0], 100 * e[ok]]), np.concatenate([[0], s[ok]]), color=colors[i % 5], lw=1.5,
                label=rf"$\ell_0$ = {float(c['l0']):g} µm, $G_c$ = {float(c['Gc']):.0f} MPa µm")
        fs = frac.get("fracture_strain_node")
        if fs is not None:
            j = int(np.argmin(np.abs(e - fs)))
            ax.plot([100 * e[j]], [s[j]], "o", ms=5, color=colors[i % 5], mfc="white")
        table.append({"l0_um": float(c["l0"]), "Gc_MPa_um": float(c["Gc"]), "Gc_l0": float(c["Gc"]) * float(c["l0"]),
                      "fracture_strain": fs, "fracture_strain_element": frac.get("fracture_strain_element"),
                      "peak_eng_stress_MPa": frac.get("peak_eng_stress"), "strain_at_peak": frac.get("strain_at_peak"),
                      "fracture_energy_MPa": frac.get("fracture_energy_MPa"), "final_max_damage": frac.get("final_max_damage"),
                      "steps": frac.get("steps"), "newton_iterations": frac.get("newton_iterations"), "alternations": frac.get("alternations"),
                      "trials": sum(1 for r in log if float(r["l0"]) == float(c["l0"])), "dir": str(d)})
    if markers:
        ax.plot([100 * e for e, _ in markers], [s_ for _, s_ in markers], "x", ms=5, color="#111111", label="experiment markers (Vaughan 2024)")
    ax.axvline(100 * target, color="#888888", lw=0.8, ls=":")
    ax.set_xlim(0, 20.5); ax.set_ylim(0, 700); ax.grid(False)
    ax.set_xlabel(r"Nominal strain $\varepsilon$ (%)"); ax.set_ylabel(r"Engineering stress $\sigma$ (MPa), gross section")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(figs / f"{prefix}_curves.png"); plt.close(fig)
    with (out / f"{prefix}_table.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0])); w.writeheader(); w.writerows(table)

    # ---- bisection record: fracture strain against Gc per l0 ----
    fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=200)
    for i, l0 in enumerate(sorted({float(r["l0"]) for r in log})):
        rows = [r for r in log if float(r["l0"]) == l0]
        gc = np.array([float(r["Gc"]) for r in rows])
        fs = np.array([float(r["fracture_strain"]) if r["fracture_strain"] not in ("", "None") else np.nan for r in rows])
        o = np.argsort(gc)
        ax.plot(gc[o], 100 * fs[o], "o-", color=colors[i % 5], ms=4, lw=1.0, label=rf"$\ell_0$ = {l0:g} µm")
        none = np.isnan(fs)
        if none.any():
            ax.plot(gc[none], np.full(none.sum(), 20.5), "^", color=colors[i % 5], ms=5, mfc="white")
    ax.axhline(100 * target, color="#888888", lw=0.8, ls=":")
    ax.set_xscale("log"); ax.set_ylim(0, 21); ax.grid(False)
    ax.set_xlabel(r"$G_c$ (MPa µm)"); ax.set_ylabel(r"fracture strain $\varepsilon_f$ (%)  (open triangles: none by 20%)")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(figs / f"{prefix}_bisection.png"); plt.close(fig)

    # ---- Gc against l0 ----
    if len(table) >= 2:
        fig, ax = plt.subplots(figsize=(4.6, 3.6), dpi=200)
        l0 = np.array([t["l0_um"] for t in table]); gc = np.array([t["Gc_MPa_um"] for t in table])
        ax.plot(l0, gc, "o-", color="#1f4e9c", ms=5, lw=1.2, label=r"calibrated $G_c$")
        ax.plot(l0, gc * l0 / l0[0], "s--", color="#d1495b", ms=4, lw=1.0, label=r"$G_c \ell_0$ / $\ell_{0,\min}$")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(False)
        ax.set_xlabel(r"$\ell_0$ (µm)"); ax.set_ylabel(r"$G_c$ (MPa µm)"); ax.legend(frameon=False)
        fig.tight_layout(); fig.savefig(figs / f"{prefix}_gc_l0.png"); plt.close(fig)
    for t in table:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in t.items() if k != "dir"})
    print("wrote", out / f"{prefix}_table.csv", "and figures", figs / f"{prefix}_*.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
