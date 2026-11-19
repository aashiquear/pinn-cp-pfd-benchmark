#!/usr/bin/env python3
"""Part IV, stage D (phase-field crystal): figures and measures.

Usage: python3 scripts/cryst_pfc_figures.py [d1 | d2 | series | poly | all]
Reads results/crystallization/pfc/<run>/ (VTK series of PINNACLE, cached as frames.npz) and
the reference files of scripts/cryst_pfc_reference.py; writes figures/crystallization/d*.png,
d*.gif and results/crystallization/pfc/*_summary.csv.

Measures
  crystalline fraction   share of the box where the Gaussian average (width a/1.5) of
                         (psi - psi0)^2 exceeds half of its maximum in the frame
  equivalent radius      of the full crystal (the decks solve a quarter): sqrt(4 X A_box / pi)
  free-energy density    `energy` of solver_history.csv divided by the box area
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from scipy import ndimage  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cryst_pfc_collect import load  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/pfc"
FIG = ROOT / "figures/crystallization"
A = 4.0 * np.pi / np.sqrt(3.0)


def fraction(psi: np.ndarray, psi0: float, h: float) -> float:
    local = ndimage.gaussian_filter((psi - psi0) ** 2, sigma=(A / 1.5) / h, mode="mirror")
    return float((local > 0.5 * local.max()).mean()) if local.max() > 1e-4 else 0.0


def full(psi: np.ndarray) -> np.ndarray:
    """The quarter mirrored to the full crystal."""
    top = np.concatenate([psi[:, ::-1], psi[:, 1:]], axis=1)
    return np.concatenate([top[::-1], top[1:]], axis=0)


def gif(frames: list[np.ndarray], path: Path, vmin: float, vmax: float, size: float = 3.4) -> None:
    fig, ax = plt.subplots(figsize=(size, size * frames[0].shape[0] / frames[0].shape[1]))
    ax.set_xticks([]); ax.set_yticks([])
    im = ax.imshow(frames[0], origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
    fig.subplots_adjust(0, 0, 1, 1)
    FuncAnimation(fig, lambda k: (im.set_data(frames[k]),), frames=len(frames)).save(path, writer=PillowWriter(fps=6), dpi=80)
    plt.close(fig)


def params_of(run: str) -> tuple[float, float]:
    m = re.search(r"eps([0-9.]+)_psi(-?[0-9.]+)", run)
    return float(m.group(1)), float(m.group(2).rstrip("."))


def growth_first() -> None:
    run = "d2_growth_eps0.3_psi-0.34_fem"
    eps, psi0 = params_of(run)
    d = load(run)
    h = float(d["x"][1] - d["x"][0])
    frac = np.array([fraction(p, psi0, h) for p in d["psi"]])
    ref = np.load(RES / "reference_grow_eps0.3_psi-0.34_r12.npz")
    fig, axes = plt.subplots(1, 6, figsize=(15.5, 2.9), gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.25, 1.25]})
    picks = [0, len(d["t"]) // 5, len(d["t"]) // 2, len(d["t"]) - 1]
    for ax, k in zip(axes[:4], picks):
        ax.imshow(full(d["psi"][k]), origin="lower", cmap="gray", vmin=-0.75, vmax=0.5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"$t$ = {d['t'][k]:.0f}, $X$ = {frac[k]:.2f}", fontsize=9)
    axes[4].plot(ref["t"], ref["fraction"], "k-", lw=1.2, label="spectral reference")
    axes[4].plot(d["t"], frac, "o", ms=3.5, mfc="none", color="#c8553d", label="PINNACLE FE")
    axes[4].set_xlabel("$t$"); axes[4].set_ylabel("crystalline fraction $X$"); axes[4].legend(frameon=False, fontsize=8)
    axes[5].plot(ref["t"], ref["f"], "k-", lw=1.2, label="spectral reference")
    axes[5].plot(d["t_hist"], d["energy"] / d["area"], "-", color="#c8553d", label="PINNACLE FE")
    axes[5].set_xlabel("$t$"); axes[5].set_ylabel("free-energy density"); axes[5].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "d2_growth.png", bbox_inches="tight")
    plt.close(fig)
    gif([full(p) for p in d["psi"]], FIG / "d2_growth.gif", -0.75, 0.5)
    for tq in (300, 600, 900, 1200, 1500):
        print(tq, "X fe", round(float(np.interp(tq, d["t"], frac)), 4), "ref", round(float(np.interp(tq, ref["t"], ref["fraction"])), 4),
              "f fe", round(float(np.interp(tq, d["t_hist"], d["energy"] / d["area"])), 6), "ref", round(float(np.interp(tq, ref["t"], ref["f"])), 6))
    print("newton per step", d["newton"].mean(), "wrote figures/crystallization/d2_growth.png")


def phase_points() -> None:
    """D1: mesh series of the triangular crystal, and the phase reached from a random state."""
    area = 2.0 * A * 4.0 * np.pi
    print("triangular crystal, eps 0.3, psi0 -0.25: spectral f = 0.016382, peak - mean = 0.9157")
    rows = []
    for per_a in (8, 12, 16):
        d = load(f"d1_triangular_eps0.3_psi-0.25_h{per_a}_fem")
        f, peak = d["energy"][-1] / area, float(d["psi"][-1].max()) + 0.25
        rows.append({"elements_per_a": per_a, "f": f, "f_error": f / 0.016382 - 1, "peak": peak, "peak_error": peak / 0.9157 - 1})
        print(rows[-1])
    table = list(csv.DictReader((ROOT / "data/reference/cryst_pfc_phase_diagram.csv").open()))
    runs = sorted(RES.glob("d1_phase_eps0.3_psi*_fem"), key=lambda r: -params_of(r.name)[1])
    runs = [r for r in runs if (r / "fields.pvd").exists()]
    fig = plt.figure(figsize=(15.5, 5.4))
    grid = fig.add_gridspec(2, len(runs) + 3)
    ax = fig.add_subplot(grid[:, :3])
    eps = np.array([float(r["eps"]) for r in table])
    for key, color, label in (("liquidus", "#1b4965", "liquid | coexistence"), ("solidus", "#1b4965", None),
                              ("triangular_side", "#c8553d", "triangular | stripes"), ("stripes_side", "#c8553d", None)):
        ax.plot([float(r[f"relaxed_{key}"]) for r in table], eps, "-", color=color, lw=1.3, label=label)
        ax.plot([float(r[f"one_mode_{key}"]) for r in table], eps, ":", color=color, lw=1.0)
    ax.plot([float(r["stability_limit"]) for r in table], eps, "--", color="0.5", lw=0.9, label="liquid stability limit")
    ax.plot([], [], "k:", lw=1.0, label="one-mode approximation")
    summary = []
    for k, run in enumerate(runs):
        e, psi0 = params_of(run.name)
        d = load(run.name)
        last = d["psi"][-1]
        spread = float(last.max() - last.min())
        skew = float(((last - last.mean()) ** 3).mean() / max(last.std(), 1e-12) ** 3)
        kind = "liquid" if spread < 0.05 else ("stripes" if abs(skew) < 0.25 else ("triangular" if skew > 0.6 else "mixed"))
        summary.append({"run": run.name, "eps": e, "psi0": psi0, "spread": spread, "skewness": skew, "state": kind,
                        "f": d["energy"][-1] / d["area"]})
        ax.plot(psi0, e, {"liquid": "o", "stripes": "s", "triangular": "^", "mixed": "D"}[kind], color="k", mfc="w", ms=6)
        for row, frame in ((0, len(d["psi"]) // 4), (1, len(d["psi"]) - 1)):
            a = fig.add_subplot(grid[row, 3 + k])
            a.imshow(d["psi"][frame], origin="lower", cmap="gray", vmin=psi0 - 0.75, vmax=psi0 + 1.0)
            a.set_xticks([]); a.set_yticks([])
            a.set_title(f"$\\psi_0$ = {psi0:g}, $t$ = {d['t'][frame]:.0f}" + (f"\n{kind}" if row else ""), fontsize=8)
    ax.set_xlabel("mean density $\\psi_0$"); ax.set_ylabel("reduced temperature $\\epsilon$")
    ax.set_xlim(-0.6, 0.02); ax.set_ylim(0.0, 0.52); ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "d1_phase_diagram.png", bbox_inches="tight")
    plt.close(fig)
    for name, data in (("d1_mesh_summary.csv", rows), ("d1_phase_summary.csv", summary)):
        with (RES / name).open("w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
            wr.writeheader(); wr.writerows(data)
    for r in summary:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote figures/crystallization/d1_phase_diagram.png")


def envelope(psi: np.ndarray, psi0: float, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Radius of the crystal outline (half-maximum of the local order) against the polar angle."""
    h = float(x[1] - x[0])
    local = ndimage.gaussian_filter((psi - psi0) ** 2, sigma=(A / 1.5) / h, mode="mirror")
    inside = local > 0.5 * local.max()
    angles = np.deg2rad(np.arange(1.0, 90.0, 2.0))
    radii = []
    rr = np.arange(0.0, min(x[-1], y[-1]), 0.25 * h)
    for a in angles:
        cols = np.clip(np.round(rr * np.cos(a) / h).astype(int), 0, len(x) - 1)
        rows = np.clip(np.round(rr * np.sin(a) / (y[1] - y[0])).astype(int), 0, len(y) - 1)
        line = inside[rows, cols]
        radii.append(rr[np.argmin(line)] if not line.all() else rr[-1])
    return angles, np.array(radii)


def growth_series() -> None:
    """D2: growth against mean density (eps = 0.3) and against eps (equal supersaturation)."""
    table = {float(r["eps"]): r for r in csv.DictReader((ROOT / "data/reference/cryst_pfc_phase_diagram.csv").open())}
    density = ["d2_growth_eps0.3_psi-0.36_fem", "d2_growth_eps0.3_psi-0.36_r24_fem", "d2_growth_eps0.3_psi-0.35_fem",
               "d2_growth_eps0.3_psi-0.35_r24_fem", "d2_growth_eps0.3_psi-0.34_fem", "d2_growth_eps0.3_psi-0.33_fem",
               "d2_growth_eps0.3_psi-0.32_fem"]
    eps_runs = ["d2_growth_eps0.1_psi-0.2001_r16_fem", "d2_growth_eps0.1_psi-0.2001_r32_fem", "d2_growth_eps0.2_psi-0.28_r16_fem",
                "d2_growth_eps0.3_psi-0.34_r16_fem", "d2_growth_eps0.4_psi-0.3894_r16_fem", "d2_growth_eps0.5_psi-0.4328_r16_fem"]
    checks = ["d2_growth_eps0.3_psi-0.34_dt1_fem", "d2_growth_eps0.3_psi-0.34_h14_fem"]
    rows = []
    for title, runs, out in (("density", density, "d2_density_series.png"), ("eps", eps_runs, "d2_eps_series.png")):
        runs = [r for r in runs if (RES / r / "fields.pvd").exists()]
        fig, axes = plt.subplots(2, len(runs) + 2, figsize=(2.55 * (len(runs) + 2), 5.3), squeeze=False)
        big = fig.add_subplot(axes[0, 0].get_gridspec()[:, len(runs):])
        for a in axes[:, len(runs):].ravel():
            a.remove()
        colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(runs)))
        for k, run in enumerate(runs):
            eps, psi0 = params_of(run)
            d = load(run)
            h = float(d["x"][1] - d["x"][0])
            frac = np.array([fraction(p, psi0, h) for p in d["psi"]])
            radius = np.sqrt(4.0 * frac * d["area"] / np.pi)
            row = table[eps]
            liq, sol = float(row["relaxed_liquidus"]), float(row["relaxed_solidus"])
            sat = (psi0 - liq) / (sol - liq)
            grew = frac[-1] > 2.0 * frac[0]
            # speed of the equivalent radius between X = 0.05 and 0.25 (before the box is felt)
            window = (frac > 0.05) & (frac < 0.25)
            speed = float(np.polyfit(d["t"][window], radius[window], 1)[0]) if window.sum() > 3 else float("nan")
            k_mid = int(np.argmin(np.abs(frac - 0.25))) if grew else len(frac) - 1
            ang, rad = envelope(d["psi"][k_mid], psi0, d["x"], d["y"])
            rows.append({"series": title, "run": run, "eps": eps, "psi0": psi0, "supersaturation": sat,
                         "final_fraction": float(frac[-1]), "lever_rule": min(max(sat, 0.0), 1.0), "radius_speed": speed,
                         "outline_rmax_over_rmin": float(rad.max() / max(rad.min(), 1e-9)) if grew else float("nan"),
                         "newton_per_step": float(d["newton"].mean())})
            seed_r = float(m.group(1)) if (m := re.search(r"_r(\d+)_", run)) else 12.0
            rows[-1]["seed_radius"] = seed_r
            label = (f"$\\psi_0$ = {psi0:g}" if title == "density" else f"$\\epsilon$ = {eps:g}") + (
                f", $R_0$ = {seed_r:g}" if seed_r not in (12.0, 16.0) or (title == "eps" and eps == 0.1) else "")
            for r_i, frame in ((0, k_mid), (1, len(frac) - 1)):
                lo, hi = float(d["psi"][0].min()), float(d["psi"][0].max())   # the seed sets the gray scale
                axes[r_i, k].imshow(full(d["psi"][frame]), origin="lower", cmap="gray", vmin=lo - 0.1, vmax=hi + 0.2)
                axes[r_i, k].set_xticks([]); axes[r_i, k].set_yticks([])
                axes[r_i, k].set_title(f"{label}, $t$ = {d['t'][frame]:.0f}\n$X$ = {frac[frame]:.2f}", fontsize=8)
            big.plot(d["t"] if title == "density" else d["t"] / d["t"][-1], frac, "-", color=colors[k], label=label + f" ($S$ = {sat:.2f})")
            ref = RES / f"reference_grow_eps{eps:g}_psi{psi0:g}_r{seed_r:g}.npz"
            if ref.exists():
                rd = np.load(ref)
                big.plot(rd["t"] if title == "density" else rd["t"] / d["t"][-1], rd["fraction"], ":", color=colors[k], lw=1.0)
                rows[-1]["reference_fraction_at_end"] = float(np.interp(d["t"][-1], rd["t"], rd["fraction"]))
            gif([full(p) for p in d["psi"]], FIG / f"{run.replace('_fem', '')}.gif", lo - 0.1, hi + 0.2)
        big.set_xlabel("$t$" if title == "density" else "$t / t_{end}$"); big.set_ylabel("crystalline fraction $X$")
        big.plot([], [], "k:", lw=1.0, label="spectral reference")
        big.legend(frameon=False, fontsize=7.5)
        fig.tight_layout()
        fig.savefig(FIG / out, bbox_inches="tight")
        plt.close(fig)
        print("wrote", out)
    base = load("d2_growth_eps0.3_psi-0.34_fem")
    hb = float(base["x"][1] - base["x"][0])
    fb = np.array([fraction(p, -0.34, hb) for p in base["psi"]])
    ref = np.load(RES / "reference_grow_eps0.3_psi-0.34_r12.npz")
    print("X at t = 600, 1200: reference", np.interp([600, 1200], ref["t"], ref["fraction"]).round(4),
          "h10 dt2", np.interp([600, 1200], base["t"], fb).round(4))
    for run in checks:
        if (RES / run / "fields.pvd").exists():
            d = load(run)
            h = float(d["x"][1] - d["x"][0])
            f = np.array([fraction(p, -0.34, h) for p in d["psi"]])
            print(run, np.interp([600, 1200], d["t"], f).round(4))
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0].keys()).index(k) if k in rows[0] else 99)
    with (RES / "d2_series_summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader(); wr.writerows(rows)
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})


def polycrystal() -> None:
    """D3: grains, and the particles with five or seven neighbours between them."""
    from scipy.spatial import Delaunay

    run = "d3_polycrystal_eps0.3_psi-0.32_fem"
    d = load(run)
    x, y = d["x"], d["y"]
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.0))
    picks = [0, len(d["t"]) // 6, len(d["t"]) // 3, len(d["t"]) - 1]
    counts = []
    for ax, k in zip(axes, picks):
        psi = d["psi"][k]
        ax.imshow(psi, origin="lower", cmap="gray", extent=(x[0], x[-1], y[0], y[-1]))
        peaks = (psi == ndimage.maximum_filter(psi, size=5)) & (psi > -0.32 + 0.3 * (psi.max() + 0.32))
        iy, ix = np.nonzero(peaks)
        pts = np.column_stack([x[ix], y[iy]])
        five = seven = 0
        if len(pts) > 10:
            tri = Delaunay(pts)
            # neighbours within 1.4 lattice constants (drops the long edges of the hull)
            nb = [set() for _ in pts]
            for simplex in tri.simplices:
                for a in range(3):
                    p, q = simplex[a], simplex[(a + 1) % 3]
                    if np.hypot(*(pts[p] - pts[q])) < 1.4 * A:
                        nb[p].add(q); nb[q].add(p)
            coord = np.array([len(n) for n in nb])
            margin = 1.2 * A
            interior = (pts[:, 0] > margin) & (pts[:, 0] < x[-1] - margin) & (pts[:, 1] > margin) & (pts[:, 1] < y[-1] - margin)
            # only particles surrounded by crystal: every neighbour is itself well coordinated
            bulk = interior & np.array([all(len(nb[q]) >= 5 for q in n) and len(n) >= 5 for n in nb])
            five, seven = int((coord[bulk] == 5).sum()), int((coord[bulk] == 7).sum())
            ax.plot(*pts[bulk & (coord == 5)].T, "o", ms=5, mfc="none", color="#e63946")
            ax.plot(*pts[bulk & (coord == 7)].T, "s", ms=5, mfc="none", color="#2a9d8f")
        counts.append({"t": float(d["t"][k]), "particles": int(len(pts)), "five_fold": five, "seven_fold": seven})
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"$t$ = {d['t'][k]:.0f}: {len(pts)} particles, {five} five-, {seven} seven-fold", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIG / "d3_polycrystal.png", bbox_inches="tight")
    plt.close(fig)
    gif(list(d["psi"]), FIG / "d3_polycrystal.gif", float(d["psi"][-1].min()), float(d["psi"][-1].max()), size=4.2)
    with (RES / "d3_summary.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(counts[0].keys()))
        wr.writeheader(); wr.writerows(counts)
    print(counts, "newton per step", d["newton"].mean())


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("d1", "all"):
        phase_points()
    if what in ("d2", "all"):
        growth_first()
    if what in ("series", "all"):
        growth_series()
    if what in ("poly", "all"):
        polycrystal()


if __name__ == "__main__":
    main()
