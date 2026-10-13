#!/usr/bin/env python3
"""Compare load-step predictor runs (PINN+FEM hybrid) with the cold FE run.

Deck-driven (PINNACLE ``.i`` syntax):

  [Compare]
    out_dir   = ../../results/cp_study/hybrid          # where the tables go
    figures   = ../../figures/report                    # where the figures go
    prefix    = hybrid_h40                              # file stem
  [Runs]
    [./cold]     dir = ../../results/cp_study/hybrid/h40_cold     label = 'cold FE' []
    [./linear]   dir = ../../results/cp_study/hybrid/h40_linear   label = 'linear predictor' []
    [./network]  dir = ../../results/cp_study/hybrid/h40_network  label = 'network predictor' []
  []

The first run is the reference. For every run the script reads
``solver/solver_history.csv`` (accepted steps) and ``results.csv`` (the curve,
from ``polycrystal_plate.py ... curve``) and reports:

* that the curves coincide (max |stress difference| against the reference,
  which must be at solver tolerance: the corrector decides the answer);
* Newton iterations per step and in total, Newton wall time, predictor wall
  time, training time, and their sums: the cost side;
* the predictor's relative error per step (0 = exact proposal, 1 = no
  better than the previous step): the quality side;
* for a network trained on the assembled residual (``mode = fem_residual``),
  the residual-to-reaction ratio of the proposal before and after training
  (``predictor_residual_ratio_initial`` / ``predictor_residual_ratio``), the
  assemblies spent and the steps the executioner accepted without Newton
  (``corrector_used = none``).

Writes ``<prefix>_summary.csv`` (one row per run), ``<prefix>_steps.csv`` (per
step, every run) and three figures: iterations/time per step, the predictor
error along the load path, and the residual ratio before/after training.

Usage: hybrid_compare.py decks/report/hybrid_compare.i
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("compare", "runs"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _str(v) -> str:
    return " ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def read_history(d: Path) -> list[dict]:
    rows = list(csv.DictReader((d / "solver" / "solver_history.csv").open()))
    out = []
    for r in rows:
        conv = str(r.get("converged", "1")).strip().lower()
        if conv not in ("1", "1.0", "true") and not conv.startswith("1.0"):
            continue
        out.append({k: (float(v) if v not in ("", None) else np.nan) for k, v in r.items()})
    return out


def read_curve(d: Path):
    p = d / "results.csv"
    if not p.exists():
        return None
    rows = list(csv.DictReader(p.open()))
    return np.array([float(r["nominal_strain"]) for r in rows]), np.array([float(r["eng_stress"]) for r in rows])


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve()
    cfg = read_deck(deck)
    comp = cfg["compare"]
    out_dir = _resolve(str(comp["out_dir"]), deck); out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = _resolve(str(comp.get("figures", comp["out_dir"])), deck); fig_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(comp.get("prefix", "hybrid"))
    runs = {name: (_resolve(str(r["dir"]), deck), _str(r.get("label", name))) for name, r in cfg["runs"].items()}
    missing = [n for n, (d, _) in runs.items() if not (d / "solver" / "solver_history.csv").exists()]
    for n in missing:
        print(f"skipping {n}: no solver_history.csv yet")
        runs.pop(n)
    names = list(runs)
    ref = names[0]
    hist = {n: read_history(d) for n, (d, _) in runs.items()}
    curves = {n: read_curve(d) for n, (d, _) in runs.items()}

    summary = []
    step_rows = []
    for n in names:
        h = hist[n]
        iters = np.array([r.get("iterations", np.nan) for r in h])
        newton = np.array([r.get("newton_time_s", np.nan) for r in h])
        pred_t = np.array([r.get("predictor_time_s", 0.0) for r in h])
        pred_e = np.array([r.get("predictor_rel_error", np.nan) for r in h])
        train_t = np.array([r.get("predictor_train_time_s", 0.0) for r in h])
        asm_t = np.nan_to_num(np.array([r.get("predictor_assembly_time_s", 0.0) for r in h]))
        asm_n = np.nan_to_num(np.array([r.get("predictor_assemblies", 0.0) for r in h]))
        rr0 = np.array([r.get("predictor_residual_ratio_initial", np.nan) for r in h])
        rr1 = np.array([r.get("predictor_residual_ratio", np.nan) for r in h])
        used = [str(r.get("corrector_used", "")) for r in h]
        pred_t = np.nan_to_num(pred_t); train_t = np.nan_to_num(train_t)
        curve_dev = np.nan
        if curves[n] is not None and curves[ref] is not None:
            e0, s0 = curves[ref]; e, s = curves[n]
            curve_dev = float(np.max(np.abs(np.interp(e0, e, s) - s0)))
        summary.append({
            "run": n, "label": runs[n][1], "steps": len(h),
            "newton_iterations": int(np.nansum(iters)),
            "newton_time_s": float(np.nansum(newton)),
            "predictor_time_s": float(np.sum(pred_t)),
            "train_time_s": float(np.sum(train_t)),
            "total_solve_time_s": float(np.nansum(newton) + np.sum(pred_t)),
            "predictor_rel_error_median": float(np.nanmedian(pred_e)) if np.isfinite(pred_e).any() else np.nan,
            "predictor_rel_error_max": float(np.nanmax(pred_e)) if np.isfinite(pred_e).any() else np.nan,
            "residual_ratio_initial_median": float(np.nanmedian(rr0)) if np.isfinite(rr0).any() else np.nan,
            "residual_ratio_median": float(np.nanmedian(rr1)) if np.isfinite(rr1).any() else np.nan,
            "assemblies": int(np.sum(asm_n)),
            "assembly_time_s": float(np.sum(asm_t)),
            "steps_accepted_without_newton": int(sum(1 for u in used if u == "none")),
            "max_curve_dev_vs_ref_MPa": curve_dev,
        })
        for r, it, nt, pt, pe, a0, a1, u in zip(h, iters, newton, pred_t, pred_e, rr0, rr1, used, strict=True):
            step_rows.append({"run": n, "step": int(r["step"]), "load_factor": r["load_factor"],
                              "iterations": it, "newton_time_s": nt, "predictor_time_s": pt,
                              "predictor_rel_error": pe, "residual_ratio_initial": a0,
                              "residual_ratio": a1, "corrector_used": u})
    with (out_dir / f"{prefix}_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0])); w.writeheader(); w.writerows(summary)
    with (out_dir / f"{prefix}_steps.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(step_rows[0])); w.writeheader(); w.writerows(step_rows)
    for s in summary:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()})

    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8,
                         "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    colors = ["#111111", "#1f4e9c", "#d1495b", "#00798c", "#edae49", "#6a4c93", "#2a9d8f"]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), dpi=200)
    for i, n in enumerate(names):
        h = hist[n]
        lam = 100.0 * 0.20 * np.array([r["load_factor"] for r in h])
        it = np.array([r.get("iterations", np.nan) for r in h])
        t = np.array([r.get("newton_time_s", np.nan) for r in h]) + np.nan_to_num(np.array([r.get("predictor_time_s", 0.0) for r in h]))
        axes[0].plot(lam, it, lw=1.2, color=colors[i % len(colors)], label=runs[n][1])
        axes[1].plot(lam, t, lw=1.2, color=colors[i % len(colors)], label=runs[n][1])
    axes[0].set_ylabel("Newton iterations per step"); axes[1].set_ylabel("solve + predictor time per step (s)")
    for ax in axes:
        ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)"); ax.set_xlim(0, 20.5); ax.grid(False); ax.legend(frameon=False)
    axes[0].set_ylim(0, None)
    fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_cost.png"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.6, 3.4), dpi=200)
    for i, n in enumerate(names[1:], start=1):
        h = hist[n]
        lam = 100.0 * 0.20 * np.array([r["load_factor"] for r in h])
        pe = np.array([r.get("predictor_rel_error", np.nan) for r in h])
        ax.plot(lam, pe, lw=1.2, color=colors[i % len(colors)], label=runs[n][1])
    ax.axhline(1.0, color="#777", lw=0.8)
    ax.set_yscale("log"); ax.set_xlim(0, 20.5)
    ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)"); ax.set_ylabel(r"$\|u_{pred}-u\| / \|u-u_{prev}\|$")
    ax.grid(False); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_prediction_error.png"); plt.close(fig)
    # Residual ratio before/after training, for the runs that report it.
    with_ratio = [n for n in names if any(np.isfinite(r.get("predictor_residual_ratio", np.nan)) for r in hist[n])]
    if with_ratio:
        fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=200)
        for i, n in enumerate(with_ratio):
            h = hist[n]
            lam = 100.0 * 0.20 * np.array([r["load_factor"] for r in h])
            a0 = np.array([r.get("predictor_residual_ratio_initial", np.nan) for r in h])
            a1 = np.array([r.get("predictor_residual_ratio", np.nan) for r in h])
            c = colors[names.index(n) % len(colors)]
            ax.plot(lam, a0, lw=1.0, ls="--", color=c, label=f"{runs[n][1]}: before training (baseline)")
            ax.plot(lam, a1, lw=1.4, color=c, label=f"{runs[n][1]}: after training")
        tol = float(comp.get("corrector_tol", 0.0))
        if tol > 0:
            ax.axhline(tol, color="#777", lw=0.8, ls=":")
        ax.set_yscale("log"); ax.set_xlim(0, 20.5)
        ax.set_xlabel(r"Engineering strain $\varepsilon$ (%)"); ax.set_ylabel(r"$\|R\|_{free} / \|R\|_{Dirichlet}$")
        ax.grid(False); ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
        fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_residual_ratio.png", bbox_inches="tight"); plt.close(fig)
    print("wrote", out_dir / f"{prefix}_summary.csv", fig_dir / f"{prefix}_cost.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
