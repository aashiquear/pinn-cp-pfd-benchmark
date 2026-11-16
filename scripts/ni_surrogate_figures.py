#!/usr/bin/env python3
"""Figures of the Ni specimen surrogate: held-out synthetic specimens and the 15 real specimens.

Deck-driven (PINNACLE .i syntax):

  [Figures]
    out_dir  = ../../../figures/ni/final
    prefix   = surrogate
    model    = ../../../results/ni_study/tl/final/surrogate/model_h020     # parametric_surrogate output
    targets  = ../../../results/ni_study/tl/final/surrogate/targets_h020   # samples.csv, real_specimens.csv
    fem      = ../../../results/ni_study/tl/final/fem                      # FE route, <specimen>/summary.json
    parity   = ../../../data/digitized/ni_fig4b_failure_strain_parity.csv  # experiment column only
    fracture = ../../../data/digitized/ni_fig5_fracture_descriptors.csv    # experiment, Fig. 5 panels
    dpi      = 200
  []

Writes <out_dir>/<prefix>_holdout_curves.png (engineering stress-strain curve of every held-out
synthetic specimen: FE run, surrogate, material-point baseline, with both failure strains),
<prefix>_holdout_parity.png (failure strain, fracture location and peak stress: FE against
surrogate, training and held-out samples), <prefix>_real.png (the 15 real specimens: surrogate
against the FE route and against the experiment) and <prefix>_table.csv. Only the experiment's
columns of the data files are read. No solver run.

Usage: ni_surrogate_figures.py FIGURES.i
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ni_final_figures import GROUP_COLORS, VARIANT_MARKERS, _marker, experiment, predicted  # noqa: E402
from report_figures import _resolve, _str  # noqa: E402

FE, SUR, BASE = "#111111", "#d1495b", "#999999"
SCALARS = (("failure_gauge_strain", "strain to failure"), ("location_fraction", "fracture location (of the gage)"),
           ("peak_eng_stress", "peak engineering stress (MPa)"))


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "figures" not in cfg:
        raise SystemExit(f"{path}: missing [Figures]")
    return cfg["figures"]


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(path.open()))


def _rms(v) -> float:
    v = np.asarray(v, dtype=float)
    return float(np.sqrt(np.mean(v**2))) if v.size else float("nan")


def holdout_curves(model: Path, names: dict[int, str], rows: list[dict], out: Path, dpi: int) -> None:
    hold = [r for r in rows if r["split"] == "holdout"]
    n = len(hold)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.9 * nrow), squeeze=False, dpi=dpi, sharex=True, sharey=True)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    for ax, r in zip(axes.ravel(), hold):
        i = int(r["sample"])
        c = np.genfromtxt(model / "holdout" / f"{i:04d}.csv", delimiter=",", names=True)
        ax.plot(c["gauge_strain"], c["baseline"], "-", color=BASE, lw=0.9, label="material point (network input)")
        ax.plot(c["gauge_strain"], c["eng_stress_true"], "-", color=FE, lw=1.4, label="FE run")
        ax.plot(c["gauge_strain"], c["eng_stress_pred"], "--", color=SUR, lw=1.4, label="surrogate")
        ft, fp = float(r["s_failure_gauge_strain_true"]), float(r["s_failure_gauge_strain_pred"])
        ax.axvline(ft, color=FE, lw=0.8, ls=":")
        ax.axvline(fp, color=SUR, lw=0.8, ls=":")
        ax.set_title("%s: failure %.3f (FE), %.3f (surrogate); curve error %.1f%%"
                     % (names.get(i, str(i)), ft, fp, 100 * float(r["rel_l2"])), fontsize=7)
        ax.set_xlim(0, 0.65); ax.set_ylim(0, 260)
    for ax in axes[-1]:
        ax.set_xlabel("gauge engineering strain")
    for ax in axes[:, 0]:
        ax.set_ylabel("engineering stress (MPa)")
    axes[0, 0].legend(loc="lower right", frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def holdout_parity(rows: list[dict], out: Path, dpi: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8), dpi=dpi)
    for ax, (key, label) in zip(axes, SCALARS):
        for split, kw, lab in (("train", dict(marker="o", color="#888888", mfc="none", ms=5, ls="none"), "training (24)"),
                               ("holdout", dict(marker="o", color=SUR, ms=6, ls="none"), "held out (6)")):
            sel = [r for r in rows if r["split"] == split]
            t = np.array([float(r[f"s_{key}_true"]) for r in sel]); p = np.array([float(r[f"s_{key}_pred"]) for r in sel])
            err = _rms((p - t) / t) if key != "location_fraction" else _rms(p - t)
            unit = "%.1f%%" % (100 * err) if key != "location_fraction" else "%.3f of the gage" % err
            ax.plot(t, p, label=f"{lab}: RMS error {unit}", **kw)
        allv = np.array([float(r[f"s_{key}_true"]) for r in rows])
        lo, hi = (0.0, 1.0) if key == "location_fraction" else (0.92 * allv.min(), 1.06 * allv.max())
        ax.plot([lo, hi], [lo, hi], "-", color="#111111", lw=0.9, zorder=0)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.set_xlabel(f"FE run: {label}", fontsize=8); ax.set_ylabel("surrogate", fontsize=8)
        ax.legend(loc="upper left", frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def real(table: list[dict], out: Path, dpi: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.9), dpi=dpi)
    panels = (("fem_failure_strain", "sur_failure_strain", "FE route: strain to failure", (0.3, 0.7), True),
              ("exp_failure_strain", "sur_failure_strain", "experimental strain to failure (Fig. 4(b))", (0.3, 0.7), True),
              ("exp_location", "sur_location", "experimental fracture location (of the gage, Fig. 5)", (0.0, 1.0), False))
    for ax, (kx, ky, xlabel, lim, rel) in zip(axes, panels):
        ax.plot(lim, lim, "-", color="#111111", lw=0.9)
        if rel:
            ax.fill_between(lim, [0.9 * x for x in lim], [1.1 * x for x in lim], color="#dddddd", lw=0, label=r"$\pm$10%")
        else:
            ax.fill_between(lim, [lim[0] - 0.08, lim[1] - 0.08], [lim[0] + 0.08, lim[1] + 0.08], color="#dddddd", lw=0,
                            label=r"$\pm$0.08 (reading of Fig. 5)")
        errs = []
        for r in table:
            c, m = _marker(r["specimen"])
            ax.plot(r[kx], r[ky], color=c, marker=m, ms=8 if m == "*" else 6, mfc="none" if m == "o" else c, mew=1.4, ls="none")
            errs.append((r[ky] - r[kx]) / r[kx] if rel else r[ky] - r[kx])
        ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("surrogate: strain to failure" if rel else "surrogate: fracture location", fontsize=8)
        ax.set_title("RMS difference %.1f%%" % (100 * _rms(errs)) if rel else "RMS difference %.2f of the gage" % _rms(errs), fontsize=9)
    for g, c in GROUP_COLORS.items():
        axes[0].plot([], [], "s", color=c, label=f"group {g}")
    for v, m in VARIANT_MARKERS.items():
        axes[0].plot([], [], m, color="#555555", mfc="none" if m == "o" else "#555555", label=f"variant x.{v}")
    for ax in axes:
        ax.legend(loc="upper left", frameon=False, fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    f = read_deck(deck)
    out_dir = _resolve(_str(f["out_dir"]), deck)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix, dpi = _str(f.get("prefix", "surrogate")), int(f.get("dpi", 200))
    model, targets = _resolve(_str(f["model"]), deck), _resolve(_str(f["targets"]), deck)
    rows = [r for r in _rows(model / "surrogate_holdout.csv") if r["split"] in ("train", "holdout")]
    names = {int(r["sample"]): r["specimen"] for r in _rows(targets / "samples.csv")}
    holdout_curves(model, names, rows, out_dir / f"{prefix}_holdout_curves.png", dpi)
    holdout_parity(rows, out_dir / f"{prefix}_holdout_parity.png", dpi)
    exp = experiment(f, deck)
    fem = predicted(_resolve(_str(f["fem"]), deck))
    specs = [r["specimen"] for r in _rows(targets / "real_specimens.csv")]
    preds = _rows(model / "surrogate_predictions.csv")
    if len(specs) != len(preds):
        raise SystemExit("real_specimens.csv and surrogate_predictions.csv differ in length")
    table = []
    for spec, p in zip(specs, preds):
        e, fe = exp[spec], fem[spec]
        table.append({"specimen": spec, "exp_failure_strain": e["strain"], "fem_failure_strain": fe["strain"],
                      "sur_failure_strain": float(p["s_failure_gauge_strain"]), "exp_location": e["location"],
                      "fem_location": fe["location"], "sur_location": float(p["s_location_fraction"]),
                      "fem_peak_MPa": fe["peak"], "sur_peak_MPa": float(p["s_peak_eng_stress"]),
                      "ensemble_spread": float(p["spread"]), "flagged": int(p["flagged"])})
    real(table, out_dir / f"{prefix}_real.png", dpi)
    with (out_dir / f"{prefix}_table.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0])); w.writeheader(); w.writerows(table)
    a = lambda k: np.array([r[k] for r in table])  # noqa: E731
    sel = np.array([r["specimen"].split(".")[0] in ("L", "P", "10") for r in table])
    for lab, m in (("all 15", np.ones(len(table), bool)), ("L, P, 10", sel)):
        print("%s: surrogate vs FE %.3f, surrogate vs experiment %.3f, FE vs experiment %.3f (RMS relative, failure strain); "
              "peak stress vs FE %.3f; location vs FE %.3f, vs experiment %.3f (RMS, of the gage)"
              % (lab, _rms((a("sur_failure_strain")[m] - a("fem_failure_strain")[m]) / a("fem_failure_strain")[m]),
                 _rms((a("sur_failure_strain")[m] - a("exp_failure_strain")[m]) / a("exp_failure_strain")[m]),
                 _rms((a("fem_failure_strain")[m] - a("exp_failure_strain")[m]) / a("exp_failure_strain")[m]),
                 _rms((a("sur_peak_MPa")[m] - a("fem_peak_MPa")[m]) / a("fem_peak_MPa")[m]),
                 _rms(a("sur_location")[m] - a("fem_location")[m]), _rms(a("sur_location")[m] - a("exp_location")[m])))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
