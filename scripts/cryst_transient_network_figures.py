#!/usr/bin/env python3
"""Cases A2 and A5, network routes against the FE march on the same window: table and figures.

Usage: python3 scripts/cryst_transient_network_figures.py
Reads results/crystallization/single_field/<case>_{window_fem,pinn_seed*,hybrid_seed0}/ for
the cases of scripts/cryst_transient_network_decks.py; writes
results/crystallization/single_field/network_transient_table.csv and
figures/crystallization/{a2_network,a5_network}.png.

Measures. `amount` is int phi dx (A2, nm: the front position up to a constant) or
int phi 4 pi r^2 dr (A5, nm^3). PINN: at the slab ends (summary.json) and, for seed 0, at
every step through the hybrid run's `predictor_mass`. Front velocity (A2): slope of the
amount over t >= 2. Newton iterations: sum over the steps, FE march against the hybrid.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/single_field"
FIG = ROOT / "figures/crystallization"
V_REF = 0.644490  # travelling-wave reference of case A2 (scripts/cryst_reference.py)
SEED_COLORS = ("#1f4e79", "#c8553d", "#3a7d44")


def history(run: str):
    path = RES / run / "solver_history.csv"
    return np.genfromtxt(path, delimiter=",", names=True) if path.exists() else None


def slab_amounts(run: str):
    path = RES / run / "summary.json"
    if not path.exists():
        return None
    inv = json.load(path.open())["phase_field"]["invariants"]
    return np.array([[e["t"], e["mass"]] for e in inv])


def final_profile(run: str, radial: bool):
    m = read_vtk(RES / run / "fields.vtk")
    x = m["nodes"][:, 0]
    key = "phi" if "phi" in m["point"] else "displacement"
    u = np.asarray(m["point"][key], dtype=float)
    u = u[:, 0] if u.ndim == 2 else u
    order = np.argsort(x)
    return x[order], u[order]


def amount_of(x: np.ndarray, u: np.ndarray, radial: bool) -> float:
    w = 4.0 * np.pi * x**2 if radial else np.ones_like(x)
    return float(np.trapezoid(u * w, x))


def case_figure(cases: list[tuple[str, str]], radial: bool, out: str, rows: list[dict]) -> None:
    fig, axes = plt.subplots(len(cases), 3, figsize=(12.0, 3.3 * len(cases)), squeeze=False)
    unit = r"solid volume $\int\phi\,4\pi r^2dr$ (nm$^3$)" if radial else r"solid length $\int\phi\,dx$ (nm)"
    for i, (case, label) in enumerate(cases):
        fe = history(f"{case}_window_fem")
        hyb = history(f"{case}_hybrid_seed0")
        ax = axes[i, 0]
        if fe is not None:
            ax.plot(fe["t"], fe["mass"], color="black", lw=1.4, label="FE")
        if hyb is not None and "predictor_mass" in hyb.dtype.names:
            ax.plot(hyb["t"], hyb["predictor_mass"], color=SEED_COLORS[0], lw=1.0, ls="--",
                    label="PINN, seed 0 (every step)")
        for seed, color in enumerate(SEED_COLORS):
            pts = slab_amounts(f"{case}_pinn_seed{seed}")
            if pts is None:
                continue
            ax.plot(pts[:, 0], pts[:, 1], "o", ms=4, mfc="none", color=color, label=f"PINN, seed {seed} (slab ends)")
            entry = {"case": case, "route": "PINN", "seed": seed, "amount_end": pts[-1, 1],
                     "amount_end_fe": float(fe["mass"][-1]) if fe is not None else "",
                     "velocity": "", "velocity_rel_error": "", "newton_total": "", "pred_rms_mean": "", "pred_rms_max": ""}
            if not radial:
                sel = pts[:, 0] >= 2.0
                v = float(np.polyfit(pts[sel, 0], pts[sel, 1], 1)[0])
                entry.update(velocity=v, velocity_rel_error=v / V_REF - 1.0)
            rows.append(entry)
        ax.set_xlabel(r"$t\,M_\phi$ (GPa$^{-1}$)"); ax.set_ylabel(unit); ax.set_title(label, fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        if fe is not None:
            entry = {"case": case, "route": "FE", "seed": "", "amount_end": float(fe["mass"][-1]), "amount_end_fe": "",
                     "velocity": "", "velocity_rel_error": "", "newton_total": int(fe["newton_iters"].sum()),
                     "pred_rms_mean": "", "pred_rms_max": ""}
            if not radial:
                sel = fe["t"] >= 2.0
                v = float(np.polyfit(fe["t"][sel], fe["mass"][sel], 1)[0])
                entry.update(velocity=v, velocity_rel_error=v / V_REF - 1.0)
            rows.append(entry)
        ax = axes[i, 1]
        if fe is not None:
            ax.plot(fe["t"], fe["newton_iters"], color="black", lw=2.2, label=f"FE ({int(fe['newton_iters'].sum())} in total)")
        if hyb is not None:
            ax.plot(hyb["t"], hyb["newton_iters"], color=SEED_COLORS[0], lw=1.0,
                    label=f"PINN+FEM ({int(hyb['newton_iters'].sum())} in total)")
            twin = ax.twinx()
            twin.semilogy(hyb["t"], hyb["predictor_rms"], color="0.55", lw=0.8)
            twin.set_ylabel("rms distance of the prediction (grey)")
            entry = {"case": case, "route": "PINN+FEM", "seed": 0, "amount_end": float(hyb["mass"][-1]),
                     "amount_end_fe": float(fe["mass"][-1]) if fe is not None else "", "velocity": "", "velocity_rel_error": "",
                     "newton_total": int(hyb["newton_iters"].sum()), "pred_rms_mean": float(hyb["predictor_rms"].mean()),
                     "pred_rms_max": float(hyb["predictor_rms"].max())}
            if not radial:
                sel = hyb["t"] >= 2.0
                v = float(np.polyfit(hyb["t"][sel], hyb["mass"][sel], 1)[0])
                entry.update(velocity=v, velocity_rel_error=v / V_REF - 1.0)
                vp = float(np.polyfit(hyb["t"][sel], hyb["predictor_mass"][sel], 1)[0])
                rows.append({**entry, "route": "PINN (seed 0, every step)", "velocity": vp,
                             "velocity_rel_error": vp / V_REF - 1.0, "newton_total": "",
                             "amount_end": float(hyb["predictor_mass"][-1])})
            rows.append(entry)
        ax.set_xlabel(r"$t\,M_\phi$ (GPa$^{-1}$)"); ax.set_ylabel("Newton iterations per step")
        ax.set_ylim(0, None); ax.legend(fontsize=7, frameon=False, loc="lower left")
        ax = axes[i, 2]
        if (RES / f"{case}_window_fem/fields.vtk").exists():
            x, u = final_profile(f"{case}_window_fem", radial)
            ax.plot(x, u, color="black", lw=1.4, label="FE")
        # the network's field at the end of the horizon, on the FE nodes (written by the hybrid run)
        hyb_file = RES / f"{case}_hybrid_seed0/fields.vtk"
        if hyb_file.exists():
            m = read_vtk(hyb_file)
            if "field_phi_pinn" in m["point"]:
                order = np.argsort(m["nodes"][:, 0])
                ax.plot(m["nodes"][order, 0], np.asarray(m["point"]["field_phi_pinn"])[order],
                        color=SEED_COLORS[0], lw=1.0, ls="--", label="PINN, seed 0")
        ax.set_xlabel(r"$r$ (nm)" if radial else r"$x$ (nm)"); ax.set_ylabel(r"$\phi$ at the end of the horizon")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / out, bbox_inches="tight")
    plt.close(fig)
    print("wrote figures/crystallization/" + out)


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    rows: list[dict] = []
    case_figure([("a2_front_velocity_standard", "A2: front, $\\Delta T/T_f = 0.20$")], False, "a2_network.png", rows)
    case_figure([("a5_seed_subcritical", "A5: seed of 0.9 critical radii"),
                 ("a5_seed_supercritical", "A5: seed of 1.1 critical radii")], True, "a5_network.png", rows)
    if rows:
        with (RES / "network_transient_table.csv").open("w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader(); wr.writerows(rows)
        for r in rows:
            print(" ".join(f"{k}={v:.5g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))


if __name__ == "__main__":
    main()
