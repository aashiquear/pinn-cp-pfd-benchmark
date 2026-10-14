#!/usr/bin/env python3
"""Fit the J2 isotropic-hardening law to the digitized experimental curve (0-D).

Deck-driven (PINNACLE ``.i`` syntax):

  [Fit]
    curve     = ../../data/digitized/al2219_fig3_calibration_stress_strain.csv
    series    = experiment_Vaughan2024          # rows of the digitized file to fit
    markers   = ../../data/digitized/al2219_fig3_experiment_markers.csv
    out_dir   = ../../results/pfd_study/j2_fit
    figure    = ../../figures/pfd/j2_fit.png
    C11 = 108.0e3   C12 = 61.3e3   C44 = 28.5e3  # MPa: E and nu from the Voigt average
    measure   = true                            # fit in true stress / log strain (finite-strain J2)
    plastic_min = 5.0e-4                        # plastic strain below which points are not fitted
    hardening = 'power piecewise'               # laws to fit: power (sigma_y0 + K alpha^n) and/or
                                                # piecewise (table of alpha, sigma_y for PINNACLE's
                                                # hardening_mode = piecewise)
    table_points = 24                           # piecewise: bins of alpha (log-spaced)

The elastic constants are NOT fitted: E and nu are the Voigt average of the
cubic constants the CP study used, so the J2 and CP plates share the same
elastic response. The digitized engineering curve is converted to true
stress and logarithmic strain (``measure = true``), the plastic strain
alpha = eps_true - sigma_true / E is formed, and sigma_y0, K, n of
sigma_y(alpha) = sigma_y0 + K alpha^n are fitted by least squares on the
points with alpha >= plastic_min; the piecewise law is the median true
stress per log-spaced bin of alpha, made non-decreasing. Writes
``<out_dir>/j2_fit.json`` (the parameters of every law in the deck's units
and the residuals at the markers), ``j2_fit_<law>.csv`` (the fitted curve
back in engineering measures) and the figure.

Usage: j2_fit.py decks/pfd_study/j2_fit.i
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    if "fit" not in cfg:
        raise SystemExit(f"{path}: missing [Fit] block")
    return cfg["fit"]


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def voigt_isotropic(C11: float, C12: float, C44: float) -> tuple[float, float]:
    """E and nu of the Voigt average of a cubic crystal."""
    K = (C11 + 2 * C12) / 3.0
    G = (C11 - C12 + 3 * C44) / 5.0
    E = 9 * K * G / (3 * K + G)
    nu = (3 * K - 2 * G) / (2 * (3 * K + G))
    return float(E), float(nu)


def fit_power(alpha: np.ndarray, sig: np.ndarray) -> tuple[float, float, float]:
    """Least squares for sigma_y0 + K alpha^n: n by scan, (sigma_y0, K) linear for each n."""
    best = None
    for n in np.linspace(0.05, 1.0, 381):
        A = np.column_stack([np.ones_like(alpha), alpha**n])
        coef, res, *_ = np.linalg.lstsq(A, sig, rcond=None)
        r = float(np.sum((A @ coef - sig) ** 2))
        if coef[1] > 0 and (best is None or r < best[0]):
            best = (r, float(coef[0]), float(coef[1]), float(n))
    assert best is not None
    return best[1], best[2], best[3]


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    deck = Path(argv[1]).resolve()
    f = read_deck(deck)
    out = _resolve(str(f["out_dir"]), deck)
    out.mkdir(parents=True, exist_ok=True)
    rows = [r for r in csv.DictReader((ln for ln in _resolve(str(f["curve"]), deck).open() if not ln.startswith("#"))) if r["series"] == str(f.get("series", "experiment_Vaughan2024"))]
    eps_e = np.array([float(r["strain_pct"]) for r in rows]) / 100.0
    sig_e = np.array([float(r["stress_MPa"]) for r in rows])
    order = np.argsort(eps_e)
    eps_e, sig_e = eps_e[order], sig_e[order]
    E, nu = voigt_isotropic(float(f["C11"]), float(f["C12"]), float(f["C44"]))
    true = bool(f.get("measure", "true") in (True, "true", "True"))
    if true:
        sig_t, eps_t = sig_e * (1.0 + eps_e), np.log1p(eps_e)
    else:
        sig_t, eps_t = sig_e, eps_e
    alpha = eps_t - sig_t / E
    keep = alpha >= float(f.get("plastic_min", 5.0e-4))
    laws = [str(x) for x in f.get("hardening", "power")] if isinstance(f.get("hardening"), (list, tuple)) else str(f.get("hardening", "power")).split()

    def make_sigma_of_eps(yield_fn):
        def sigma_of_eps(eps_eng: np.ndarray) -> np.ndarray:
            e_t = np.log1p(eps_eng) if true else eps_eng
            s = np.empty_like(e_t)
            for i, et in enumerate(e_t):
                # solve sigma = E (et - alpha) = sigma_y(alpha) for alpha >= 0
                if E * et <= yield_fn(0.0):
                    s[i] = E * et
                    continue
                lo, hi = 0.0, et
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    if E * (et - mid) > yield_fn(mid):
                        lo = mid
                    else:
                        hi = mid
                s[i] = yield_fn(0.5 * (lo + hi))
            return s / (1.0 + eps_eng) if true else s
        return sigma_of_eps

    mrows = []
    mpath = f.get("markers")
    if mpath:
        for r in csv.DictReader((ln for ln in _resolve(str(mpath), deck).open() if not ln.startswith("#"))):
            mrows.append((float(r.get("strain_pct", r.get("strain", 0))) / 100.0, float(r.get("stress_MPa", r.get("stress", 0)))))
    grid = np.linspace(0.0, float(eps_e.max()), 400)
    plastic = alpha >= float(f.get("plastic_min", 5.0e-4))
    rec: dict = {"E_MPa": E, "nu": nu, "measure": "true" if true else "engineering", "fit_points": int(keep.sum()), "laws": {},
                 "source": {"curve": str(f["curve"]), "series": str(f.get("series", "")), "elastic": "Voigt average of C11 C12 C44 (CP study)"}}
    curves = {}
    for law in laws:
        if law == "power":
            sy0, K, n = fit_power(alpha[keep], sig_t[keep])
            yield_fn = lambda a, sy0=sy0, K=K, n=n: sy0 + K * max(a, 0.0) ** n  # noqa: E731
            entry = {"sigma_y0_MPa": sy0, "K_MPa": K, "n": n, "label": rf"$\sigma_y = {sy0:.0f} + {K:.0f}\,\alpha^{{{n:.2f}}}$ MPa"}
        elif law == "piecewise":
            nb = int(f.get("table_points", 24))
            a_min = float(f.get("plastic_min", 5.0e-4)); a_max = float(alpha[keep].max())
            edges = np.geomspace(a_min, a_max * 1.0001, nb + 1)
            pts_a, pts_s = [], []
            for lo_, hi_ in zip(edges[:-1], edges[1:], strict=True):
                m = (alpha >= lo_) & (alpha < hi_)
                if m.sum() >= 2:
                    pts_a.append(float(np.median(alpha[m]))); pts_s.append(float(np.median(sig_t[m])))
            # the last digitized point anchors the end of the table (a bin median
            # would stop short of the final hardening), then non-decreasing
            pts_a.append(float(alpha[keep].max())); pts_s.append(float(sig_t[keep][np.argmax(alpha[keep])]))
            pts_s = np.maximum.accumulate(np.array(pts_s)).tolist()
            pts_a = [0.0] + pts_a; pts_s = [pts_s[0]] + pts_s
            yield_fn = lambda a, A=np.array(pts_a), S=np.array(pts_s): float(np.interp(a, A, S))  # noqa: E731
            entry = {"plastic_strain_points": pts_a, "yield_stress_points": pts_s, "label": f"J2 piecewise table ({len(pts_a)} points)"}
        else:
            raise SystemExit(f"unknown hardening law {law!r}")
        sig_of = make_sigma_of_eps(yield_fn)
        fit_e = sig_of(grid); curves[law] = fit_e
        pred = sig_of(eps_e)
        entry["rms_plastic_MPa"] = float(np.sqrt(np.mean((pred[plastic] - sig_e[plastic]) ** 2)))
        entry["max_abs_dev_plastic_MPa"] = float(np.max(np.abs(pred[plastic] - sig_e[plastic])))
        entry["markers"] = [{"strain_pct": 100 * e, "stress_exp_MPa": s_, "stress_fit_MPa": float(sig_of(np.array([e]))[0]),
                             "rel_error_pct": 100 * (float(sig_of(np.array([e]))[0]) - s_) / s_} for e, s_ in mrows]
        rec["laws"][law] = entry
        with (out / f"j2_fit_{law}.csv").open("w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["nominal_strain", "eng_stress_fit"]); w.writerows(zip(grid, fit_e, strict=True))
    (out / "j2_fit.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("laws", "source")}))
    for law, e in rec["laws"].items():
        print(law, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in e.items() if k not in ("markers", "plastic_strain_points", "yield_stress_points", "label")})
        for m in e["markers"]:
            print("  ", {k: round(v, 2) for k, v in m.items()})
    if f.get("figure"):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "mathtext.fontset": "dejavusans",
                             "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
        fig, ax = plt.subplots(figsize=(6.0, 4.2), dpi=200)
        ax.plot(np.concatenate([[0], 100 * eps_e]), np.concatenate([[0], sig_e]), color="#d1495b", lw=1.6, label="experiment (digitized)")
        for (law, fit_e), ls, col in zip(curves.items(), ("--", "-."), ("#1f4e9c", "#00798c"), strict=False):
            ax.plot(100 * grid, fit_e, color=col, lw=1.4, ls=ls, label=rec["laws"][law]["label"])
        if mrows:
            ax.plot([100 * e for e, _ in mrows], [s_ for _, s_ in mrows], "o", ms=4, color="#111111", label="experiment markers")
        ax.set_xlim(0, 20.5); ax.set_ylim(0, 700); ax.grid(False)
        ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)"); ax.set_ylabel(r"Engineering stress $\sigma$ (MPa)"); ax.legend(loc="lower right", frameon=False)
        fp = _resolve(str(f["figure"]), deck); fp.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout(); fig.savefig(fp); plt.close(fig)
        print("wrote", fp)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
