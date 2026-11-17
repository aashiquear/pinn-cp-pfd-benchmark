#!/usr/bin/env python3
"""Figures and animations of Part IV, stage A (single-field model, cases A1-A5).

Usage: python3 scripts/cryst_single_field_figures.py decks/crystallization/figures/single_field.i

Reads the PINNACLE run directories (solver_history.csv, fields.vtk, the fields_step*.vtk
series), the reference tables of scripts/cryst_reference.py and the digitized curve of the
review's Fig. 3(c); writes PNG figures and GIF animations to figures/crystallization/.
Figure rules of the study: no grid, closed axes box, symbols for quantities.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from cryst_reference import read_deck  # noqa: E402
from vtk_legacy import read_vtk, series_files  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.grid": False, "figure.dpi": 150, "savefig.dpi": 200,
                     "axes.linewidth": 0.8, "legend.frameon": False, "mathtext.fontset": "cm"})
K_B = 1.380649e-23
C = {"fem": "#1f4e79", "pinn": "#c0392b", "hybrid": "#2e8b57", "ref": "black", "pub": "#888888",
     "cnt": "#b8860b"}


def history(run: Path) -> dict[str, np.ndarray]:
    rows = list(csv.DictReader((run / "solver_history.csv").open()))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if rows[0][k] not in ("", None)
            and re.match(r"^-?[0-9.eE+-]+$|^nan$", rows[0][k])}


def profile(path: Path, name: str = "phi") -> tuple[np.ndarray, np.ndarray]:
    m = read_vtk(path)
    x = m["nodes"][:, 0]
    u = m["point"][name]
    o = np.argsort(x)
    return x[o], u[o]


def series(run: Path) -> list[tuple[float, np.ndarray, np.ndarray]]:
    """(time, x, phi) of every snapshot; times from the .pvd index."""
    times = [float(t) for t in re.findall(r'timestep="([^"]+)"', (run / "fields.pvd").read_text())]
    files = series_files(run)
    out = []
    for (step, path), t in zip(files, times):
        x, u = profile(path)
        out.append((t, x, u))
    return out


def save(fig, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(w_pad=1.6)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out.relative_to(ROOT))


def main() -> None:
    deck = read_deck(Path(sys.argv[1]).resolve())
    res = ROOT / deck["Paths"]["results"]
    ref_dir = ROOT / deck["Paths"]["reference"]
    figs = ROOT / deck["Paths"]["figures"]
    mat = {k: float(v) for k, v in deck["Material"].items()}
    delta = mat["delta_nm"]

    # ------------------------------------------------------------------ A1
    run = res / "a1_planar_interface_fem"
    h = history(run)
    x, u = profile(run / "fields.vtk")
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.2))
    xs = np.linspace(-2, 2, 400)
    ax[0].plot(xs, 1 / (1 + np.exp(xs / delta)), color=C["ref"], lw=1.0, label="closed form")
    k = (np.abs(x) <= 2) & (np.arange(len(x)) % 8 == 0)
    ax[0].plot(x[k], u[k], "o", ms=3.2, mfc="none", color=C["fem"], label="PINNACLE, FE")
    s0 = series(run)[0]
    ax[0].plot(s0[1], s0[2], ":", color=C["pub"], lw=1.0, label="initial profile")
    ax[0].set_xlim(-2, 2); ax[0].set_ylim(-0.03, 1.03)
    ax[0].set_xlabel(r"$x$ (nm)"); ax[0].set_ylabel(r"$\phi$"); ax[0].legend(loc="upper right")
    ax[1].plot(h["t"], h["energy"], color=C["fem"], lw=1.3, label="PINNACLE, FE")
    ax[1].axhline(mat["gamma"], color=C["ref"], lw=0.8, ls="--", label=r"$\gamma_{SL}=0.302$ J m$^{-2}$")
    ax[1].set_xlim(0, h["t"][-1]); ax[1].set_ylim(0.29, 0.38)
    ax[1].set_xlabel(r"$t\,M_\phi$ (GPa$^{-1}$)"); ax[1].set_ylabel(r"$F/A$ (J m$^{-2}$)"); ax[1].legend()
    save(fig, figs / "a1_planar_interface.png")

    # ------------------------------------------------------------------ A2
    front = {r["interpolation"]: float(r["v_nm_per_time"]) for r in
             csv.DictReader((ref_dir / "ni_front/front_table.csv").open())
             if abs(float(r["reduced_undercooling"]) - 0.2) < 1e-9}
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.2))
    labels = {"tanh": r"$p=\phi^2(3-2\phi)$", "standard": r"$p=\phi^3(10-15\phi+6\phi^2)$"}
    styles = {"tanh": "-", "standard": "--"}
    table = []
    for tag in ("tanh", "standard"):
        hh = history(res / f"a2_front_velocity_{tag}_fem")
        pos = hh["mass"]  # int phi dx = solid length
        sel = hh["t"] > 5
        v = np.polyfit(hh["t"][sel], pos[sel], 1)[0]
        table.append((tag, v, front[tag]))
        ax[0].plot(hh["t"], pos, styles[tag], color=C["fem"], lw=1.3, label=f"FE, {labels[tag]}")
        ax[0].plot(hh["t"], pos[sel][0] + front[tag] * (hh["t"] - hh["t"][sel][0]), color=C["ref"], lw=0.6)
        vel = np.gradient(pos, hh["t"])
        ax[1].plot(hh["t"], vel, styles[tag], color=C["fem"], lw=1.3, label=f"FE, {labels[tag]}")
        ax[1].axhline(front[tag], color=C["ref"], lw=0.6)
    ax[0].plot([], [], color=C["ref"], lw=0.6, label="reference slope")
    ax[0].set_xlim(0, 35); ax[0].set_ylim(0, 32)
    ax[0].set_xlabel(r"$t\,M_\phi$ (GPa$^{-1}$)"); ax[0].set_ylabel(r"solid length $\int\phi\,dx$ (nm)")
    ax[0].legend(loc="upper left")
    ax[1].set_xlim(0, 35); ax[1].set_ylim(0.60, 0.70)
    ax[1].set_xlabel(r"$t\,M_\phi$ (GPa$^{-1}$)"); ax[1].set_ylabel(r"$v/M_\phi$ (nm GPa)")
    ax[1].plot([], [], color=C["ref"], lw=0.6, label="reference"); ax[1].legend(loc="lower right")
    save(fig, figs / "a2_front_velocity.png")
    with (res / "front_fem_table.csv").open("w") as fh:
        fh.write("interpolation,v_fem,v_reference,rel_error\n")
        for tag, v, r in table:
            fh.write(f"{tag},{v:.8f},{r:.8f},{v / r - 1:.3e}\n")

    # ------------------------------------------------------------------ A3
    ref_rows = list(csv.DictReader((ref_dir / "ni_nucleus_3d/nucleus_table.csv").open()))
    profs = json.loads((ref_dir / "ni_nucleus_3d/nucleus_profiles.json").read_text())
    std = [r for r in ref_rows if r["well"] == "quartic" and r["interpolation"] == "standard"
           and r["temperature_scaling"] == "proportional"]
    fem = list(csv.DictReader((res / "nucleus_fem_table.csv").open()))
    fem_std = [r for r in fem if r["well"] == "quartic" and r["interpolation"] == "standard"]
    dig = np.loadtxt(ROOT / deck["Paths"]["digitized"], delimiter=",", comments="#", skiprows=3)
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.3))
    # (a) radial profiles, the seven undercoolings of the published figure
    cmap = plt.get_cmap("viridis")
    shown = [0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]
    for i, xu in enumerate(shown):
        col = cmap(i / (len(shown) - 1) * 0.9)
        p = profs[f"quartic:standard:proportional:{xu:g}"]
        ax[0].plot(np.array(p["r_m"]) * 1e9, p["phi"], color=col, lw=1.0)
        r, u = profile(res / f"a3_ni_nucleus_x{xu:.3f}_fem/fields.vtk")
        kk = (r <= 3.0) & (np.arange(len(r)) % 24 == 0)
        ax[0].plot(r[kk], u[kk], "o", ms=2.8, mfc="none", color=col)
    ax[0].plot([], [], color="k", lw=1.0, label="reference (shooting)")
    ax[0].plot([], [], "o", ms=2.8, mfc="none", color="k", label="PINNACLE, FE")
    ax[0].set_xlim(0, 3); ax[0].set_ylim(0, 1.0)
    ax[0].set_xlabel(r"$r$ (nm)"); ax[0].set_ylabel(r"$\phi$"); ax[0].legend(loc="upper right")
    ax[0].text(0.52, 0.62, r"$\Delta T/T_f = 0.15$ (outer)" "\n" r"to $0.30$ (inner)", transform=ax[0].transAxes, fontsize=8)
    # (b) barrier
    xr = np.array([float(r["reduced_undercooling"]) for r in std])
    ax[1].plot(dig[:, 0], dig[:, 1], color=C["pub"], lw=4.0, alpha=0.55,
               label="published (Fig. 3c, digitized)")
    ax[1].plot(xr, [float(r["W_star_over_kT"]) for r in std], color=C["ref"], lw=0.9, label="reference (shooting)")
    ax[1].plot(xr, [float(r["W_cnt_over_kT"]) for r in std], "--", color=C["cnt"], lw=0.9,
               label=r"droplet model, $\gamma\propto T$")
    ax[1].plot([float(r["reduced_undercooling"]) for r in fem_std], [float(r["W_star_over_kT"]) for r in fem_std],
               "o", ms=4.5, mfc="none", color=C["fem"], label="PINNACLE, FE")
    for route, mk in (("pinn", "s"), ("hybrid", "^")):
        post = res / f"a3_ni_nucleus_x0.200_{route}/postprocessors.csv"
        if post.exists():
            w = float(post.read_text().splitlines()[-1])
            ax[1].plot([0.2], [w / (K_B * 1382.64 * 1e18)], mk, ms=7, mfc="none", color=C[route],
                       label={"pinn": "PINNACLE, PINN", "hybrid": "PINNACLE, PINN+FEM"}[route])
    ax[1].set_xlim(0.1, 0.3); ax[1].set_ylim(0, 120)
    ax[1].set_xlabel(r"$\Delta T/T_f$"); ax[1].set_ylabel(r"$W^*/kT$"); ax[1].legend(fontsize=7.5)
    # (c) radius
    ax[2].plot(xr, [float(r["R_equimolar_m"]) * 1e10 for r in std], color=C["ref"], lw=0.9,
                   label=r"$R_e$, reference")
    ax[2].plot(xr, [float(r["R_cnt_m"]) * 1e10 for r in std], "--", color=C["cnt"], lw=0.9,
                   label=r"$R^*=2\gamma/|\Delta f|$")
    ax[2].plot([float(r["reduced_undercooling"]) for r in fem_std],
                   [float(r["R_equimolar_nm"]) * 10 for r in fem_std], "o", ms=4.5, mfc="none",
                   color=C["fem"], label=r"$R_e$, PINNACLE FE")
    ax[2].set_xlim(0.1, 0.3); ax[2].set_ylim(0, 25)
    ax[2].set_xlabel(r"$\Delta T/T_f$"); ax[2].set_ylabel(r"radius ($\mathrm{\AA}$)"); ax[2].legend()
    save(fig, figs / "a3_ni_nucleus.png")

    # ------------------------------------------------------------------ A3 routes at 0.2
    p = profs["quartic:standard:proportional:0.2"]
    rr, uu = np.array(p["r_m"]) * 1e9, np.array(p["phi"])
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.2))
    ax[0].plot(rr, uu, color=C["ref"], lw=1.0, label="reference (shooting)")
    names = {"fem": "FE", "pinn": "PINN", "hybrid": "PINN+FEM"}
    marks = {"fem": "o", "pinn": "s", "hybrid": "^"}
    for i, route in enumerate(("fem", "pinn", "hybrid")):
        f = res / f"a3_ni_nucleus_x0.200_{route}/fields.vtk"
        if not f.exists():
            continue
        r, u = profile(f)
        kk = (r <= 3.0) & ((np.arange(len(r)) + 13 * i) % 40 == 0)
        ax[0].plot(r[kk], u[kk], marks[route], ms=3.5, mfc="none", color=C[route], label=names[route])
        sel = r <= rr[-1]
        ax[1].semilogy(r[sel], np.abs(u[sel] - np.interp(r[sel], rr, uu)) + 1e-12, color=C[route],
                       lw=1.1, label=names[route])
    ax[0].set_xlim(0, 3); ax[0].set_ylim(0, 1); ax[0].set_xlabel(r"$r$ (nm)"); ax[0].set_ylabel(r"$\phi$")
    ax[0].legend()
    ax[1].set_xlim(0, 3); ax[1].set_ylim(1e-8, 1e-2)
    ax[1].set_xlabel(r"$r$ (nm)"); ax[1].set_ylabel(r"$|\phi-\phi_{\mathrm{ref}}|$"); ax[1].legend()
    save(fig, figs / "a3_routes_x0.200.png")

    # ------------------------------------------------------------------ A4
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    pairs = [("quartic", "standard", "standard PF", "-", "black"),
             ("gl_fcc", "standard", r"$g_{fcc}$ & $p_{PF}$", "-", "#7f7f7f"),
             ("quartic", "gl_bcc", r"GL: $g_{bcc}$ & $p_{bcc}$", "--", "#1f4e79"),
             ("gl_fcc", "gl_fcc", r"GL: $g_{fcc}$ & $p_{fcc}$", "--", "#c0392b")]
    for well, interp, lab, ls, col in pairs:
        rows = [r for r in ref_rows if r["well"] == well and r["interpolation"] == interp
                and r["temperature_scaling"] == "proportional"]
        ax.plot([float(r["reduced_undercooling"]) for r in rows], [float(r["W_star_over_kT"]) for r in rows],
                ls, color=col, lw=1.1, label=lab)
        pts = [r for r in fem if r["well"] == well and r["interpolation"] == interp]
        ax.plot([float(r["reduced_undercooling"]) for r in pts], [float(r["W_star_over_kT"]) for r in pts],
                "o", ms=4.5, mfc="none", color=col)
    ax.plot([], [], "o", ms=4.5, mfc="none", color="k", label="PINNACLE, FE")
    ax.set_xlim(0.1, 0.3); ax.set_ylim(0, 150)
    ax.set_xlabel(r"$\Delta T/T_f$"); ax.set_ylabel(r"$W^*/kT$"); ax.legend(fontsize=8)
    save(fig, figs / "a4_functions.png")

    # ------------------------------------------------------------------ A5
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.2))
    for tag, ls, lab in (("sub", "--", r"seed $0.9\,R$"), ("super", "-", r"seed $1.1\,R$")):
        hh = history(res / f"a5_seed_{tag}critical_fem")
        ax[0].plot(hh["t"], hh["mass"], ls, color=C["fem"], lw=1.3, label=lab)
        ax[1].plot(hh["t"], hh["energy"], ls, color=C["fem"], lw=1.3, label=lab)
        ax[2].plot(hh["t"], hh["phi_max"], ls, color=C["fem"], lw=1.3, label=lab)
    ax[0].axhline(5.2901, color=C["ref"], lw=0.6); ax[0].text(0.2, 5.6, "nucleus", fontsize=8)
    ax[1].axhline(mat["w_star_aJ"], color=C["ref"], lw=0.6); ax[1].text(4.2, 0.97, r"$W^*$", fontsize=8)
    ax[0].set_ylim(0, 30); ax[1].set_ylim(-1.5, 1.2); ax[2].set_ylim(0, 1.03)
    for a, yl in zip(ax, (r"solid volume $\int\phi\,dV$ (nm$^3$)", r"$F$ (aJ)", r"$\phi(r=0)$")):
        a.set_xlim(0, 8); a.set_xlabel(r"$t\,M_\phi$ (GPa$^{-1}$)"); a.set_ylabel(yl); a.legend()
    save(fig, figs / "a5_seeds.png")

    # ------------------------------------------------------------------ animations
    for name, ms in deck.get("Animations", {}).items():
        run = res / name
        frames = series(run)
        fig, a = plt.subplots(figsize=(5.0, 3.2))
        (line,) = a.plot([], [], color=C["fem"], lw=1.5)
        xmax = {"a1": (-2, 2), "a2": (0, 40), "a5": (0, 5)}[name[:2]]
        a.set_xlim(*xmax); a.set_ylim(-0.03, 1.03)
        a.set_xlabel(r"$x$ (nm)" if name[:2] != "a5" else r"$r$ (nm)"); a.set_ylabel(r"$\phi$")
        title = a.set_title("")

        def draw(i, frames=frames, line=line, title=title):
            t, xx, uu2 = frames[i]
            line.set_data(xx, uu2)
            title.set_text(rf"$t\,M_\phi = {t:.2f}$ GPa$^{{-1}}$")
            return line, title

        anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
        out = figs / f"{name.replace('_fem', '')}.gif"
        anim.save(out, writer=PillowWriter(fps=max(1, int(1000 / float(ms)))), dpi=110)
        plt.close(fig)
        print("wrote", out.relative_to(ROOT), len(frames), "frames")


if __name__ == "__main__":
    main()
