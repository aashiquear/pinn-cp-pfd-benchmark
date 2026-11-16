#!/usr/bin/env python3
"""Figures of the Ni CP calibration (Part III, docs/12 §2), driven by one deck.

Reads ``decks/ni/figures/cal_figures.i`` and writes, for the calibrated
iteration (``[Figures] full_set``), the Al2219 calibration set (same style,
scripts/report_figures.py):

  cal_curve.png        calibrated FE curve against the Fig. 3 experiment markers
  cal_oat_panels.png   2 x 2: one parameter perturbed per panel
  cal_oat_overlay.png  all perturbations on one axis
  cal_error.png        per-marker absolute and relative error
  cal_sensitivity.png  dS/d(ln p) of each parameter at the markers

and ``cal_convergence.png`` (the center of every iteration listed in
``[Iterations]`` on one axis, and the RMSE against the iteration).

Curves are gauge strain (``gauge_strain``) against engineering stress
(``eng_stress``) from each case's ``results.csv``, prepended with (0, 0) and
cut at the last step before the first load cutback (docs/12 §2.4).

Usage: ni_calibration_figures.py decks/ni/figures/cal_figures.i
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_figures import COLORS, _resolve, _str  # noqa: E402  (shared style)

PARAMS = ("rho0", "k1", "y_c", "tau_c0")
IT_COLORS = ("#9a9a9a", "#d1495b", "#edae49", "#1f4e9c")


def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("figures", "iterations", "parameters"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


CUT_AT_CUTBACK = True  # set from the deck: runs made before PINNACLE ddc76c3 only


def read_curve(path: Path) -> tuple[np.ndarray, np.ndarray, float | None]:
    """Gauge strain (%) and engineering stress (MPa), with (0, 0) prepended.

    results.csv holds converged steps only. Before PINNACLE ddc76c3 a cutback
    left the time increment halved (twice the strain rate from there on), so for
    such runs the curve is cut at the first cutback (``cut_at_first_cutback``);
    the third value is the gauge strain (%) of the cut, or None."""
    rows = list(csv.DictReader(open(path)))
    cut = None
    keep = []
    for r in rows:
        if CUT_AT_CUTBACK and float(r["cutbacks"]) > 0:
            cut = 100.0 * float(r["gauge_strain"])
            break
        keep.append(r)
    e = np.array([100.0 * float(r["gauge_strain"]) for r in keep])
    s = np.array([float(r["eng_stress"]) for r in keep])
    return np.concatenate([[0.0], e]), np.concatenate([[0.0], s]), cut


def read_markers(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = list(csv.DictReader([ln for ln in open(path) if not ln.startswith("#")]))
    return (np.array([100.0 * float(r["strain"]) for r in rows]),
            np.array([float(r["stress_MPa"]) for r in rows]))


def box(ax, f: dict, ymax: float, xlabel=True, ylabel=True):
    tick = float(f["stress_tick"])
    for k in (1, 2, 4, 10, 20):  # coarser tick blocks when the curves run high
        if ymax / (k * tick) <= 8:
            break
    tick *= k
    top = max(float(f["stress_max"]), tick * (np.floor(ymax / tick) + 1))
    ax.set_xlim(0.0, float(f["strain_max"]))
    ax.set_ylim(0.0, top)
    ax.set_yticks(np.arange(0.0, top + 1e-9, tick))
    ax.grid(False)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
    if xlabel:
        ax.set_xlabel(r"Engineering (gauge) strain $\varepsilon$ (%)")
    if ylabel:
        ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)")


def exp_marker(ax, ex, ey, label=True):
    ax.plot(ex, ey, "o", ms=5.5, color=COLORS["exp"], mfc="white", mew=1.3,
            label="Experiment, specimen L (Rezwan 2026, Fig. 3)" if label else None, zorder=5)


def fmt(v: float) -> str:
    """A parameter value in mathtext, scientific below 1e-3 or above 1e4."""
    if v == 0 or 1e-3 <= abs(v) < 1e4:
        return "%.4g" % v
    m, ex = ("%.3e" % v).split("e")
    return r"%s$\times10^{%d}$" % (m, int(ex))


def set_label(center: dict, prm: dict, head: str) -> str:
    def one(n):
        u = _str(prm[n]["unit"])
        return "%s = %s%s" % (_str(prm[n]["symbol"]), fmt(center[n]), (" " + u) if u else "")
    return "%s\n%s, %s\n%s, %s" % (head, one("rho0"), one("k1"), one("y_c"), one("tau_c0"))


def cut_note(ax, cut):
    if cut is not None:
        ax.axvline(cut, color="#999999", lw=0.9, ls=":", label="first load cutback (curve ends)")


def early_end(e: np.ndarray, e0: np.ndarray, f: dict) -> str:
    """Legend suffix for a perturbation curve that ends before the center curve
    (author, 2026-10-18: say so on the figure rather than rerun), else ''."""
    if e[-1] < e0[-1] - 1.0:
        return ",\nends at %.1f%% (%s)" % (e[-1], _str(f.get("early_end_reason", "run stopped")))
    return ""


def figs_iteration(name: str, it: dict, cfg: dict, deck: Path, out: Path, draw: bool) -> dict:
    f, prm = cfg["figures"], cfg["parameters"]
    root = _resolve(_str(it["root"]), deck)
    gn = json.loads((root / "gauss_newton.json").read_text())
    center = gn["center"]
    ex, ey = read_markers(_resolve(_str(f["markers"]), deck))
    dpi = int(f["dpi"])
    e0, s0, cut0 = read_curve(root / "base" / "results.csv")
    head = "CP FE, %s: RMSE %.1f MPa" % (_str(it["label"]), gn["rmse_center_MPa"])
    summary = {"label": _str(it["label"]), "curve": (e0, s0), "rmse": gn["rmse_center_MPa"],
               "pred": gn.get("rmse_predicted_linear_MPa")}
    if not draw:
        return summary
    name = "cal"

    # center curve
    fig, ax = plt.subplots(figsize=(5.6, 4.0), dpi=dpi)
    ax.plot(e0, s0, color=COLORS["fe"], lw=1.8, label=set_label(center, prm, head))
    exp_marker(ax, ex, ey)
    cut_note(ax, cut0)
    box(ax, f, max(s0.max(), ey.max()))
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout(); fig.savefig(out / f"{name}_curve.png"); plt.close(fig)

    # one at a time
    # only the parameters the iteration perturbed (2026-11-15: tau_c0 is held in the restart)
    params = [p for p in PARAMS if (root / f"{p}_p" / "results.csv").exists()]
    curves = {p: read_curve(root / f"{p}_p" / "results.csv") for p in params}
    ymax = max([s0.max(), ey.max()] + [c[1].max() for c in curves.values()])
    if len(params) == 4:
        fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), dpi=dpi, sharex=True, sharey=True)
    else:
        fig, axes = plt.subplots(1, len(params), figsize=(4.0 * len(params), 3.8), dpi=dpi, sharey=True)
    for ax, p in zip(np.atleast_1d(axes).ravel(), params):
        e, s, _cut = curves[p]
        factor = gn["perturbations"][p]["value"] / center[p]
        u = _str(prm[p]["unit"])
        ax.plot(e0, s0, color=COLORS["fe"], lw=1.6, label="center")
        ax.plot(e, s, color=COLORS[p], lw=1.6,
                label=r"%s $\times$ %.1f = %s%s" % (_str(prm[p]["symbol"]), factor,
                                                     fmt(gn["perturbations"][p]["value"]), (" " + u) if u else "")
                + early_end(e, e0, f))
        if early_end(e, e0, f):
            ax.plot(e[-1], s[-1], "x", ms=7, mew=1.6, color=COLORS[p])
        exp_marker(ax, ex, ey, label=False)
        box(ax, f, ymax, xlabel=False, ylabel=False)
        ax.legend(loc="upper left", frameon=False)
    grid = np.atleast_2d(axes)
    for ax in grid[-1]:
        ax.set_xlabel(r"Engineering (gauge) strain $\varepsilon$ (%)")
    for ax in grid[:, 0]:
        ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)")
    fig.suptitle(_str(it["label"]), fontsize=10)
    fig.tight_layout(); fig.savefig(out / f"{name}_oat_panels.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=dpi)
    ax.plot(e0, s0, color=COLORS["fe"], lw=2.0, label="center, " + _str(it["label"]))
    for p in params:
        e, s, _ = curves[p]
        factor = gn["perturbations"][p]["value"] / center[p]
        ax.plot(e, s, color=COLORS[p], lw=1.5,
                label=r"%s $\times$ %.1f" % (_str(prm[p]["symbol"]), factor) + early_end(e, e0, f))
        if early_end(e, e0, f):
            ax.plot(e[-1], s[-1], "x", ms=7, mew=1.6, color=COLORS[p])
    exp_marker(ax, ex, ey)
    box(ax, f, ymax)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(out / f"{name}_oat_overlay.png"); plt.close(fig)

    # per-marker error
    rows = gn["per_marker_center"]
    me = np.array([100.0 * r["strain"] for r in rows])
    rel = np.array([100.0 * r["rel_err"] for r in rows])
    ab = np.array([r["fe"] - r["experiment"] for r in rows])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3), dpi=dpi)
    for ax, y, lab in ((a1, ab, r"$\sigma_{FE} - \sigma_{exp}$ (MPa)"), (a2, rel, "relative error (%)")):
        ax.axhline(0.0, color="#777", lw=0.8)
        ax.vlines(me, 0.0, y, color=COLORS["fe"], lw=1.2)
        ax.plot(me, y, "o", ms=5, color=COLORS["fe"])
        lim = 1.15 * max(abs(y).max(), 1e-9)
        ax.set_ylim(-lim, lim)
        ax.set_xlim(0.0, float(f["strain_max"]))
        ax.set_xlabel(r"Engineering (gauge) strain $\varepsilon$ (%)")
        ax.set_ylabel(lab)
        ax.grid(False)
    fig.suptitle("%s: RMSE %.1f MPa" % (_str(it["label"]), gn["rmse_center_MPa"]), fontsize=10)
    fig.tight_layout(); fig.savefig(out / f"{name}_error.png"); plt.close(fig)

    # sensitivity
    J = np.array(gn["jacobian_dS_dlnp"])
    cols = list(center)
    # a held parameter has a zero column and no perturbation curve
    fig, ax = plt.subplots(figsize=(5.6, 3.8), dpi=dpi)
    ax.axhline(0.0, color="#777", lw=0.8)
    for j, p in enumerate(cols):
        if p not in curves:
            continue
        # markers beyond the end of a perturbation curve carry a clamped value: not drawn
        ok = me <= curves[p][0][-1] + 1e-9
        ax.plot(me[ok], J[ok, j], "-o", ms=4, lw=1.5, color=COLORS[p],
                label=_str(prm[p]["symbol"]) + ("" if ok.all() else " (to %.0f%%)" % me[ok][-1]))
    ax.set_xlim(0.0, float(f["strain_max"]))
    ax.set_xlabel(r"Engineering (gauge) strain $\varepsilon$ (%)")
    ax.set_ylabel(r"$\partial\sigma / \partial \ln p$ (MPa)")
    ax.grid(False)
    ax.legend(loc="upper left", frameon=False, ncol=4)
    ax.set_title(_str(it["label"]), fontsize=10)
    fig.tight_layout(); fig.savefig(out / f"{name}_sensitivity.png"); plt.close(fig)
    return summary


def fig_convergence(its: list[dict], cfg: dict, deck: Path, out: Path) -> None:
    f = cfg["figures"]
    ex, ey = read_markers(_resolve(_str(f["markers"]), deck))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.0), dpi=int(f["dpi"]),
                                 gridspec_kw={"width_ratios": [1.6, 1.0]})
    for k, it in enumerate(its):
        e, s = it["curve"]
        a1.plot(e, s, color=IT_COLORS[k % len(IT_COLORS)], lw=1.7, label=it["label"])
    exp_marker(a1, ex, ey)
    box(a1, f, max([ey.max()] + [it["curve"][1].max() for it in its]))
    a1.legend(loc="upper left", frameon=False)
    n = np.arange(1, len(its) + 1)
    a2.semilogy(n, [it["rmse"] for it in its], "-o", color=COLORS["fe"], lw=1.5, label="FE center")
    a2.semilogy(n[1:], [it["pred"] for it in its[:-1]], "s", mfc="white", color="#777",
                label="linear prediction of the previous step")
    a2.set_xticks(n)
    # labeled decades and halves, so two iterations still show a readable scale
    vals = [it["rmse"] for it in its] + [it["pred"] for it in its[:-1] if it["pred"] is not None]
    lo, hi = 0.7 * min(vals), 1.3 * max(vals)
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200) if lo <= t <= hi]
    a2.set_ylim(lo, hi)
    a2.set_yticks(ticks, [str(t) for t in ticks])
    a2.minorticks_off()
    a2.set_xlabel("Gauss-Newton iteration")
    a2.set_ylabel("RMSE against the markers (MPa)")
    a2.grid(False)
    a2.legend(loc="upper right", frameon=False)
    fig.tight_layout(); fig.savefig(out / "cal_convergence.png"); plt.close(fig)


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    global CUT_AT_CUTBACK
    CUT_AT_CUTBACK = str(cfg["figures"].get("cut_at_first_cutback", "true")).lower() in ("true", "1", "yes")
    out = _resolve(_str(cfg["figures"]["out_dir"]), deck)
    out.mkdir(parents=True, exist_ok=True)
    full = _str(cfg["figures"]["full_set"])
    if full not in cfg["iterations"]:
        raise SystemExit(f"{deck}: full_set {full!r} is not an iteration")
    its = [figs_iteration(n, it, cfg, deck, out, n == full) for n, it in cfg["iterations"].items()]
    fig_convergence(its, cfg, deck, out)
    print("wrote", sorted(p.name for p in out.glob("*.png")))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
