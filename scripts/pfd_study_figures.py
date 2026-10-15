#!/usr/bin/env python3
"""Figures of a PFD sweep study (scripts/pfd_study.py output): the Fig. 8 analogue.

Deck-driven (PINNACLE ``.i`` syntax):

  [Figures]
    study     = ../j2/mesh_l_h40.i            # the study deck (cases and out_root)
    figures   = ../../../figures/pfd
    prefix    = j2_mesh_l_h40
    param     = l0                            # the swept key (x axis of the summary panels);
                                              # 'mesh' = cell size from the mesh name hNN; 'case' = categorical
    xlabel    = '$\\ell_0$ (µm)'
    logx      = true
    markers   = ../../../data/digitized/al2219_fig3_experiment_markers.csv
    target_strain = 0.172

Reads ``<out_root>/<case>/results.csv`` and ``fracture.json`` of every case
of the study and writes two figures after (chakraborty 2016 multiscale)
Fig. 8: (a) the engineering stress-strain curves of all cases with the
fracture point marked; (b) three panels against the swept parameter: peak
stress, fracture strain and fracture energy (the integral of the curve to
the end of the path). A case whose key is not swept (the base) is drawn at
its base value.

Usage: pfd_study_figures.py decks/pfd_study/j2/mesh_l_h40_figures.i
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

STUDY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, str(STUDY_ROOT / "scripts"))

from pfd_study import _resolve, cases_of, read_study  # noqa: E402

STYLE = {"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "mathtext.fontset": "dejavusans",
         "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True}


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "figures" not in cfg:
        raise SystemExit(f"{path}: missing [Figures] block")
    cfg["figures"]["_labels"] = {k.lower(): _label(v, k) for k, v in cfg.get("labels", {}).items()}
    return cfg["figures"]


PARAM_COLOURS = ["#1f4e9c", "#d1495b", "#00798c", "#edae49", "#6a4c93", "#2a9d8f", "#8c564b", "#e377c2"]
VALUE_STYLES = ["-", "--", "-.", ":"]


def _case_style(rows: list[dict]) -> dict[str, tuple[str, str]]:
    """(colour, linestyle) per case of a categorical study: colour by the varied parameter, style by its value."""
    groups: dict[str, list[str]] = {}
    for r in rows:
        param = next(iter(r["params"]), "")
        groups.setdefault(param, []).append(r["case"])
    out = {}
    for gi, (param, names) in enumerate(groups.items()):
        if param == "":
            for n in names:
                out[n] = ("#111111", "-")
            continue
        for vi, n in enumerate(names):
            out[n] = (PARAM_COLOURS[(gi - (1 if "" in groups else 0)) % len(PARAM_COLOURS)], VALUE_STYLES[vi % len(VALUE_STYLES)])
    return out


def _label(v, default: str) -> str:
    """A deck string that the translator may have split on spaces."""
    if v is None:
        return default
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def _base_value(study_cfg: dict, study: Path, key: str) -> float | None:
    base = _resolve(str(study_cfg["study"]["base"]), study)
    for line in base.read_text().splitlines():
        s = line.split("#", 1)[0].strip()
        if s.startswith(key) and "=" in s and s.split("=", 1)[0].strip() == key:
            try:
                return float(s.split("=", 1)[1].strip())
            except ValueError:
                return None
    return None


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    f = read_deck(deck)
    study = _resolve(str(f["study"]), deck)
    cfg = read_study(study)
    out_root = _resolve(str(cfg["study"]["out_root"]), study)
    figs = _resolve(str(f.get("figures", out_root)), deck)
    figs.mkdir(parents=True, exist_ok=True)
    prefix = str(f.get("prefix", study.stem))
    key = str(f.get("param", "l0"))
    base_val = _base_value(cfg, study, key)
    target = float(f.get("target_strain", 0.172))
    markers = []
    if f.get("markers"):
        for r in csv.DictReader((ln for ln in _resolve(str(f["markers"]), deck).open() if not ln.startswith("#"))):
            markers.append((float(r.get("strain_pct", 0)) / 100.0, float(r.get("stress_MPa", 0))))
    plt.rcParams.update(STYLE)
    cmap = plt.get_cmap("viridis")
    rows = []
    for name, kv in cases_of(cfg):
        d = out_root / name
        if not (d / "fracture.json").exists():
            print("missing", d)
            continue
        frac = json.loads((d / "fracture.json").read_text())
        rr = list(csv.DictReader((d / "results.csv").open()))
        e = np.array([float(r["nominal_strain"]) for r in rr])
        s = np.array([float(r["eng_stress"]) if r["eng_stress"] not in ("", "nan") else np.nan for r in rr])
        if key == "mesh":
            # the mesh sweep: x = background cell size h (um) read from the mesh name hNN (100/NN)
            import re as _re

            src = str(kv.get("mesh", "")) or str(_resolve(str(cfg["study"]["base"]), study).read_text())
            m = _re.search(r"h(\d+)_hole", src)
            val = 100.0 / float(m.group(1)) if m else None
        elif key == "case":
            val = None  # categorical: one point per case, in deck order
        else:
            val = float(kv[key]) if key in kv else base_val
        rows.append({"case": name, "value": val, "params": kv, "e": e, "s": s, "frac": frac})
    if key != "case":
        rows.sort(key=lambda r: (r["value"] is None, r["value"]))
    # (a) curves
    labels = f.get("_labels", {})
    styles = _case_style(rows) if key == "case" else {}
    many = len(rows) > 7
    fig, ax = plt.subplots(figsize=(8.2 if many else 6.2, 4.4), dpi=200)
    vals = [r["value"] for r in rows if r["value"] is not None]
    lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
    for r in rows:
        t = 0.5 if r["value"] is None or hi == lo else (np.log(r["value"] / lo) / np.log(hi / lo) if f.get("logx", True) in (True, "true") and lo > 0 else (r["value"] - lo) / (hi - lo))
        c, ls = styles.get(r["case"], (cmap(0.1 + 0.8 * t), "-"))
        ok = np.isfinite(r["s"])
        label = labels.get(r["case"].lower()) or (
            f"h = {r['value']:g} µm" if key == "mesh" and r["value"] is not None
            else f"{key} = {r['value']:g}" if key not in ("mesh", "case") and r["value"] is not None
            else (", ".join(f"{k} = {v}" for k, v in r["params"].items()) or "base")
        )
        r["label"] = label
        ax.plot(np.concatenate([[0], 100 * r["e"][ok]]), np.concatenate([[0], r["s"][ok]]), color=c, lw=1.4, ls=ls, label=label)
        fs = r["frac"].get("fracture_strain_node")
        if fs is not None:
            j = int(np.argmin(np.abs(r["e"] - fs)))
            ax.plot([100 * r["e"][j]], [r["s"][j]], "o", ms=5, color=c if isinstance(c, str) else c, mfc="white")
    if markers:
        ax.plot([100 * e for e, _ in markers], [s_ for _, s_ in markers], "x", ms=5, color="#111111", label="experiment markers")
    ax.axvline(100 * target, color="#888888", lw=0.8, ls=":")
    ax.set_xlim(0, 20.5); ax.set_ylim(0, 700); ax.grid(False)
    ax.set_xlabel(r"Nominal strain $\varepsilon$ (%)"); ax.set_ylabel(r"Engineering stress $\sigma$ (MPa), gross section")
    if many:
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=7)
    else:
        ax.legend(loc="lower right", frameon=False, ncol=2 if len(rows) > 4 else 1, fontsize=7.5)
    fig.tight_layout(); fig.savefig(figs / f"{prefix}_curves.png"); plt.close(fig)
    # (b) panels
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), dpi=200)
    categorical = key == "case"
    x = np.arange(len(rows), dtype=float) if categorical else np.array([r["value"] if r["value"] is not None else np.nan for r in rows])
    for ax, keyf, lab in zip(axes, ("peak_eng_stress", "fracture_strain_node", "fracture_energy_MPa"),
                             (r"peak $\sigma$ (MPa)", r"fracture strain $\varepsilon_f$", r"$\int \sigma\,d\varepsilon$ (MPa)"), strict=True):
        y = np.array([r["frac"].get(keyf) if r["frac"].get(keyf) is not None else np.nan for r in rows], dtype=float)
        ax.plot(x, y, "o-", color="#1f4e9c", ms=5, lw=1.2)
        if keyf == "fracture_strain_node":
            ax.axhline(target, color="#888888", lw=0.8, ls=":")
            none = np.isnan(y)
            if none.any():
                ax.plot(x[none], np.full(none.sum(), 0.2), "^", color="#1f4e9c", ms=5, mfc="white")
        if categorical:
            ax.set_xticks(x)
            ax.set_xticklabels([r["label"] for r in rows], rotation=45, ha="right", fontsize=6.5)
        elif f.get("logx", True) in (True, "true"):
            ax.set_xscale("log")
        ax.set_xlabel("" if categorical else _label(f.get("xlabel"), key)); ax.set_ylabel(lab); ax.grid(False)
    fig.tight_layout(); fig.savefig(figs / f"{prefix}_panels.png"); plt.close(fig)
    for r in rows:
        print(r["case"], r["value"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r["frac"].items() if k in ("fracture_strain_node", "peak_eng_stress", "fracture_energy_MPa", "steps_converged", "newton_iterations")})
    print("wrote", figs / f"{prefix}_curves.png", figs / f"{prefix}_panels.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
