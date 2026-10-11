#!/usr/bin/env python3
"""Report figures for the CP FE calibration, driven by one figure deck.

Reads ``decks/report/cp_calibration_figures.i`` (PINNACLE ``.i`` syntax) and
writes the figures of ``docs/reports/cp_fe_calibration_report.md``. Nothing
about which run, which markers or which axis limits is hard-wired here; the
deck records all of it so a figure can be rebuilt from versioned inputs.

Style rules (author's decision 2026-10-10, recorded in the deck header):
  * every stress-strain curve is prepended with the initial condition (0, 0);
  * parameters are written as mathematical symbols (mathtext);
  * no grid lines; all four axes spines drawn (closed box);
  * the box runs from (0, 0) to (strain_max, stress_max), one tick block above
    the largest computed stress;
  * engineering stress on both sides (FE column ``eng_stress``).

Figures written to ``[Figures] out_dir``:
  cal_final_curve.png     calibrated FE curve against the experiment markers
  cal_oat_panels.png      2 x 2: one parameter perturbed per panel
  cal_oat_overlay.png     all perturbations on one axis
  cal_error.png           per-marker relative error of the calibrated set
  cal_sensitivity.png     dS/d(ln p) of each parameter along the curve
  grain_structure.png     the 30-grain plate used by every run

Usage: report_figures.py decks/report/cp_calibration_figures.i
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

STUDY = Path(__file__).resolve().parent.parent

COLORS = {"fe": "#1f4e9c", "exp": "#111111", "rho0": "#d1495b", "k1": "#edae49",
          "y_c": "#00798c", "tau_c0": "#66a182"}

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8,
    "mathtext.fontset": "dejavusans", "axes.linewidth": 0.9,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
})


# --------------------------------------------------------------------------- deck
def read_deck(path: Path) -> dict:
    sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("figures", "calibrated", "perturbations"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    """Deck values quoted as a single token come back as a one-element list."""
    if isinstance(v, list):
        v = v[0] if v else ""
    return "" if str(v).strip() == "-" else str(v)


# --------------------------------------------------------------------------- data
def read_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Nominal strain (%) and engineering stress (MPa), with (0, 0) prepended."""
    rows = list(csv.DictReader(open(path)))
    e = np.array([100.0 * float(r["nominal_strain"]) for r in rows])
    s = np.array([float(r["eng_stress"]) for r in rows])
    return np.concatenate([[0.0], e]), np.concatenate([[0.0], s])


def read_markers(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = [ln for ln in open(path) if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    return (np.array([float(r["strain_pct"]) for r in rows]),
            np.array([float(r["stress_MPa"]) for r in rows]))


# --------------------------------------------------------------------------- style
def box(ax, fig_cfg: dict, xlabel=True, ylabel=True):
    ax.set_xlim(0.0, float(fig_cfg["strain_max"]))
    ax.set_ylim(0.0, float(fig_cfg["stress_max"]))
    ax.set_yticks(np.arange(0.0, float(fig_cfg["stress_max"]) + 1e-9, float(fig_cfg["stress_tick"])))
    ax.grid(False)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
    if xlabel:
        ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)")
    if ylabel:
        ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)")


def exp_marker(ax, ex, ey, label=True):
    ax.plot(ex, ey, "o", ms=5.5, color=COLORS["exp"], mfc="white", mew=1.3,
            label="Experiment (Vaughan 2024, Fig. 3 of Rezwan 2025)" if label else None, zorder=5)


def calibrated_label(cal: dict) -> str:
    return ("CP FE, calibrated set\n"
            r"$\rho_0$ = %.2f $\mu$m$^{-2}$, $k_1$ = %.4f" "\n"
            r"$y_c$ = %.4f $\mu$m, $\tau_{c0}$ = %.2f MPa"
            % (cal["rho0"], cal["k1"], cal["y_c"], cal["tau_c0"]))


# --------------------------------------------------------------------------- figures
def fig_final(cfg, deck, out):
    f = cfg["figures"]; cal = cfg["calibrated"]
    e, s = read_curve(_resolve(cal["results"], deck))
    ex, ey = read_markers(_resolve(f["markers"], deck))
    fig, ax = plt.subplots(figsize=(5.6, 4.0), dpi=int(f["dpi"]))
    ax.plot(e, s, color=COLORS["fe"], lw=1.8, label=calibrated_label(cal))
    exp_marker(ax, ex, ey)
    box(ax, f)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(out / "cal_final_curve.png"); plt.close(fig)


def fig_oat(cfg, deck, out):
    f = cfg["figures"]; cal = cfg["calibrated"]; pert = cfg["perturbations"]
    e0, s0 = read_curve(_resolve(cal["results"], deck))
    ex, ey = read_markers(_resolve(f["markers"], deck))
    names = list(pert)
    # 2 x 2 panels
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), dpi=int(f["dpi"]), sharex=True, sharey=True)
    for ax, name in zip(axes.ravel(), names):
        p = pert[name]
        e, s = read_curve(_resolve(p["results"], deck))
        unit = _str(p["unit"])
        sym = _str(p["symbol"])
        ax.plot(e0, s0, color=COLORS["fe"], lw=1.6, label="calibrated")
        ax.plot(e, s, color=COLORS[name], lw=1.6,
                label=r"%s $\times$ %.1f = %g %s" % (sym, float(p["factor"]), float(p["value"]), unit))
        exp_marker(ax, ex, ey, label=False)
        box(ax, f, xlabel=False, ylabel=False)
        ax.legend(loc="lower right", frameon=False)
    for ax in axes[1]:
        ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)")
    fig.tight_layout(); fig.savefig(out / "cal_oat_panels.png"); plt.close(fig)
    # overlay
    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=int(f["dpi"]))
    ax.plot(e0, s0, color=COLORS["fe"], lw=2.0, label="calibrated set")
    for name in names:
        p = pert[name]
        e, s = read_curve(_resolve(p["results"], deck))
        ax.plot(e, s, color=COLORS[name], lw=1.5,
                label=r"%s $\times$ %.1f" % (_str(p["symbol"]), float(p["factor"])))
    exp_marker(ax, ex, ey)
    box(ax, f)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(); fig.savefig(out / "cal_oat_overlay.png"); plt.close(fig)


def fig_error(cfg, deck, out):
    f = cfg["figures"]; cal = cfg["calibrated"]
    gn = json.loads(_resolve(cal["gauss_newton"], deck).read_text())
    rows = gn["per_marker_center"]
    e = np.array([100.0 * r["strain"] for r in rows])
    rel = np.array([100.0 * r["rel_err"] for r in rows])
    ab = np.array([r["fe"] - r["experiment"] for r in rows])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3), dpi=int(f["dpi"]))
    for ax, y, lab in ((a1, ab, r"$\sigma_{FE} - \sigma_{exp}$ (MPa)"), (a2, rel, "relative error (%)")):
        ax.axhline(0.0, color="#777", lw=0.8)
        ax.vlines(e, 0.0, y, color=COLORS["fe"], lw=1.2)
        ax.plot(e, y, "o", ms=5, color=COLORS["fe"])
        ax.set_xlim(0.0, float(f["strain_max"]))
        ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)")
        ax.set_ylabel(lab)
        ax.grid(False)
    a1.set_ylim(-50, 20); a2.set_ylim(-32, 8)
    fig.tight_layout(); fig.savefig(out / "cal_error.png"); plt.close(fig)


def fig_sensitivity(cfg, deck, out):
    f = cfg["figures"]; cal = cfg["calibrated"]; pert = cfg["perturbations"]
    gn = json.loads(_resolve(cal["gauss_newton"], deck).read_text())
    e = np.array([100.0 * r["strain"] for r in gn["per_marker_center"]])
    J = np.array(gn["jacobian_dS_dlnp"])  # (n_markers, 4) columns rho0, k1, y_c, tau_c0
    cols = list(gn["center"])
    fig, ax = plt.subplots(figsize=(5.6, 3.8), dpi=int(f["dpi"]))
    ax.axhline(0.0, color="#777", lw=0.8)
    for j, name in enumerate(cols):
        ax.plot(e, J[:, j], "-o", ms=4, lw=1.5, color=COLORS[name], label=_str(pert[name]["symbol"]))
    ax.set_xlim(0.0, float(f["strain_max"]))
    ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)")
    ax.set_ylabel(r"$\partial\sigma / \partial \ln p$ (MPa)")
    ax.grid(False)
    ax.legend(loc="upper left", frameon=False, ncol=4)
    fig.tight_layout(); fig.savefig(out / "cal_sensitivity.png"); plt.close(fig)


def fig_grains(cfg, deck, out):
    f = cfg["figures"]; cal = cfg["calibrated"]
    z = np.load(_resolve(cal["grains"], deck))
    gid, xs, ys = z["gid"], z["xs"], z["ys"]
    rng = np.random.default_rng(3)
    n = int(gid.max()) + 1
    lut = rng.permutation(n)
    cmap = plt.get_cmap("tab20", 20)
    img = cmap(lut[gid] % 20)
    edge = np.zeros_like(gid, dtype=bool)
    edge[1:, :] |= gid[1:, :] != gid[:-1, :]
    edge[:, 1:] |= gid[:, 1:] != gid[:, :-1]
    img[edge] = (0.12, 0.12, 0.12, 1.0)
    fig, ax = plt.subplots(figsize=(6.4, 3.6), dpi=int(f["dpi"]))
    ax.imshow(np.transpose(img, (1, 0, 2)), origin="lower",
              extent=[xs[0], xs[-1], ys[0], ys[-1]], interpolation="nearest")
    ax.set_xlabel(r"$x$ ($\mu$m), loading axis"); ax.set_ylabel(r"$y$ ($\mu$m)")
    ax.set_xticks([0, 25, 50, 75, 100]); ax.set_yticks([0, 25, 50])
    ax.annotate("", xy=(112, 25), xytext=(102, 25), annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#111"))
    ax.text(107, 29, r"$u_x$", ha="center", fontsize=10)
    ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ys[0], ys[-1])
    fig.tight_layout(); fig.savefig(out / "grain_structure.png"); plt.close(fig)


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    out = _resolve(cfg["figures"]["out_dir"], deck)
    out.mkdir(parents=True, exist_ok=True)
    for fn in (fig_final, fig_oat, fig_error, fig_sensitivity, fig_grains):
        fn(cfg, deck, out)
    print("wrote", sorted(p.name for p in out.glob("*.png")))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
