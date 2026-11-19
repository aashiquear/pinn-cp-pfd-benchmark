#!/usr/bin/env python3
"""Part IV, stage D: independent reference for the phase-field crystal model (no PINNACLE).

Usage:
    python3 scripts/cryst_pfc_reference.py cell EPS PSI0      free energy of the three phases
    python3 scripts/cryst_pfc_reference.py diagram            coexistence densities, eps list
    python3 scripts/cryst_pfc_reference.py grow EPS PSI0 [N_A N_ROWS RADIUS T_END DT AMP PER_A]

Method: Fourier pseudo-spectral in a periodic box, first-order semi-implicit in time (the
linear operator implicit, psi^3 explicit):
    psi_k(t+dt) = [psi_k - dt k^2 (psi^3)_k] / [1 + dt k^2 (-eps + (1 - k^2)^2)].
Free energy density f = < psi/2 (-eps + (1 + lap)^2) psi + psi^4/4 >.

`cell`: a unit cell (a x a sqrt(3), 32 x 56 points) relaxed at fixed mean density; the cell
size is scaled to the minimum of f. `diagram` (cells at scale 1): f(psi0) of liquid, triangular and striped
phase on a psi0 grid; coexistence by the common tangent (equal mu = df/dpsi0 and equal
f - mu psi0), from the relaxed cells and from the one-mode approximation
(elder 2004 plastic). Output: data/reference/cryst_pfc_phase_diagram.csv.
`grow`: the seed of scripts/cryst_pfc_decks.py at the centre of the periodic box
(n_a a) x (n_rows 4 pi), i.e. the mirror image of the quarter box the decks solve; writes
results/crystallization/pfc/reference_grow_eps<eps>_psi<psi0>_r<radius>.npz (t, crystal fraction, f).
AMP is the seed amplitude (the decks print the command with it).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
A = 4.0 * np.pi / np.sqrt(3.0)
Q = np.sqrt(3.0) / 2.0


def wavenumbers(nx: int, ny: int, lx: float, ly: float) -> np.ndarray:
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, lx / nx)
    ky = 2.0 * np.pi * np.fft.rfftfreq(ny, ly / ny)
    return kx[:, None] ** 2 + ky[None, :] ** 2


def energy_density(psi: np.ndarray, k2: np.ndarray, eps: float) -> float:
    pk = np.fft.rfft2(psi)
    lin = np.fft.irfft2((-eps + (1.0 - k2) ** 2) * pk, s=psi.shape)
    return float(np.mean(0.5 * psi * lin + 0.25 * psi**4))


def relax(psi: np.ndarray, k2: np.ndarray, eps: float, dt: float, steps: int, tol: float = 0.0) -> np.ndarray:
    denom = 1.0 + dt * k2 * (-eps + (1.0 - k2) ** 2)
    pk = np.fft.rfft2(psi)
    for _ in range(steps):
        nk = np.fft.rfft2(psi**3)
        new = (pk - dt * k2 * nk) / denom
        if tol and np.max(np.abs(new - pk)) / psi.size < tol:
            pk = new
            break
        pk = new
        psi = np.fft.irfft2(pk, s=psi.shape)
    return np.fft.irfft2(pk, s=psi.shape)


def tri_amplitude(eps: float, psi0: float) -> float:
    root = 15.0 * eps - 36.0 * psi0**2
    if root < 0:
        return 0.0
    return 0.8 * (psi0 - np.sqrt(root) / 3.0) if psi0 <= 0 else 0.8 * (psi0 + np.sqrt(root) / 3.0)


def tri_pattern(x: np.ndarray, y: np.ndarray, amp: float, scale: float = 1.0) -> np.ndarray:
    return amp * (np.cos(Q * x / scale) * np.cos(0.5 * y / scale) - 0.5 * np.cos(y / scale))


def cell_energy(eps: float, psi0: float, phase: str, scale: float, n: tuple[int, int] = (32, 56)) -> tuple[float, float]:
    """(f, amplitude max - mean) of the relaxed phase in a unit cell scaled by `scale`."""
    lx, ly = A * scale, 4.0 * np.pi * scale
    x = (np.arange(n[0]) * lx / n[0])[:, None]
    y = (np.arange(n[1]) * ly / n[1])[None, :]
    if phase == "triangular":
        psi = psi0 + tri_pattern(x, y, tri_amplitude(eps, psi0) or -0.1, scale)
    else:
        amp2 = 4.0 * (eps - 3.0 * psi0**2) / 3.0
        psi = psi0 + np.sqrt(max(amp2, 1e-4)) * np.cos(y / scale) + 0.0 * x
    k2 = wavenumbers(n[0], n[1], lx, ly)
    psi = relax(psi, k2, eps, 0.5, 20000, tol=1e-15)
    return energy_density(psi, k2, eps), float(psi.max() - psi0)


def phase_energy(eps: float, psi0: float, phase: str) -> tuple[float, float, float]:
    """(f, amplitude, optimal scale) with the cell size at the minimum of f."""
    if phase == "liquid":
        return 0.5 * (1.0 - eps) * psi0**2 + 0.25 * psi0**4, 0.0, 1.0
    res = minimize_scalar(lambda s: cell_energy(eps, psi0, phase, s)[0], bounds=(0.9, 1.1), method="bounded",
                          options={"xatol": 1e-5})
    f, amp = cell_energy(eps, psi0, phase, res.x)
    return f, amp, float(res.x)


def one_mode(eps: float, psi0: float, phase: str) -> float:
    liquid = 0.5 * (1.0 - eps) * psi0**2 + 0.25 * psi0**4
    if phase == "liquid":
        return liquid
    if phase == "stripes":
        a2 = 4.0 * (eps - 3.0 * psi0**2) / 3.0
        return liquid if a2 <= 0 else liquid + 0.25 * (3.0 * psi0**2 - eps) * a2 + 3.0 * a2**2 / 32.0
    # psi0 + A P, P = cos(Qx) cos(y/2) - cos(y)/2: every wave has |k| = 1, so (1 + lap)^2 P = 0 and
    # f = f_liquid + (3 psi0^2 - eps)/2 A^2 <P^2> + psi0 A^3 <P^3> + A^4 <P^4>/4,
    # <P^2> = 3/8, <P^3> = -3/16, <P^4> = 45/128; minimal at tri_amplitude (elder 2004 plastic)
    amp = tri_amplitude(eps, psi0)
    return liquid + 3.0 / 16.0 * (3.0 * psi0**2 - eps) * amp**2 - 3.0 / 16.0 * psi0 * amp**3 + 45.0 / 512.0 * amp**4


def main() -> None:
    cmd = sys.argv[1]
    if cmd == "cell":
        eps, psi0 = float(sys.argv[2]), float(sys.argv[3])
        for phase in ("liquid", "triangular", "stripes"):
            f, amp, scale = phase_energy(eps, psi0, phase)
            print(f"{phase:11s} f = {f:.8f}  one-mode {one_mode(eps, psi0, phase):.8f}  peak - mean = {amp:.5f}  cell scale {scale:.5f}")
    elif cmd == "diagram":
        import csv

        out = ROOT / "data/reference/cryst_pfc_phase_diagram.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        rows = []
        for eps in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
            lim = np.sqrt(eps / 3.0)
            grid = np.linspace(-1.9 * lim, 0.0, 58)
            # cells at the one-mode size (scale 1): the optimum differs by about 1e-3 and changes f
            # in the seventh digit (checked with `cell`)
            liquid = lambda p, e=eps: 0.5 * (1.0 - e) * p**2 + 0.25 * p**4  # noqa: E731
            relaxed = {}
            for ph in ("triangular", "stripes"):
                vals = np.array([cell_energy(eps, p, ph, 1.0) for p in grid])
                ok = vals[:, 1] > 1e-3          # the periodic state survived the relaxation
                relaxed[ph] = (grid[ok], vals[ok, 0])
            approx = {
                "triangular": (grid[15.0 * eps - 36.0 * grid**2 > 0],
                               np.array([one_mode(eps, p, "triangular") for p in grid[15.0 * eps - 36.0 * grid**2 > 0]])),
                "stripes": (grid[eps - 3.0 * grid**2 > 0],
                            np.array([one_mode(eps, p, "stripes") for p in grid[eps - 3.0 * grid**2 > 0]])),
            }
            row = {"eps": eps, "stability_limit": -lim}
            for label, tab in (("relaxed", relaxed), ("one_mode", approx)):
                tri, stripes = (curve(*tab[ph]) for ph in ("triangular", "stripes"))
                pl, pt = common_tangent((liquid, lambda p, e=eps: (1.0 - e) * p + p**3), tri, -1.2 * lim, -lim)
                p1, p2 = common_tangent(tri, stripes, -0.6 * lim, -0.4 * lim)
                row[f"{label}_liquidus"], row[f"{label}_solidus"] = pl, pt
                row[f"{label}_triangular_side"], row[f"{label}_stripes_side"] = p1, p2
            rows.append(row)
            print({k: round(float(v), 5) for k, v in row.items()}, flush=True)
        with out.open("w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
        print("wrote", out)
    elif cmd == "grow":
        grow(*[float(v) for v in sys.argv[2:]])


def curve(x: np.ndarray, f: np.ndarray):
    """(f, f') callables from samples of one phase (cubic spline on its range of existence)."""
    from scipy.interpolate import CubicSpline

    spline = CubicSpline(x, f)
    return (lambda p: spline(np.clip(p, x[0], x[-1])) + spline(np.clip(p, x[0], x[-1]), 1) * (p - np.clip(p, x[0], x[-1])),
            lambda p: spline(np.clip(p, x[0], x[-1]), 1))


def common_tangent(a, b, guess_a: float, guess_b: float) -> tuple[float, float]:
    """Densities of two coexisting phases: equal mu = f' and equal f - mu psi0."""
    from scipy.optimize import fsolve

    (fa, da), (fb, db) = a, b

    def eqs(v):
        x, y = v
        return [da(x) - db(y), (fa(x) - da(x) * x) - (fb(y) - db(y) * y)]

    sol, _, ok, _ = fsolve(eqs, [guess_a, guess_b], full_output=True, xtol=1e-12)
    return (float(sol[0]), float(sol[1])) if ok == 1 else (float("nan"), float("nan"))


def grow(eps: float, psi0: float, n_a: float = 24, n_rows: float = 10, radius: float = 12.0, t_end: float = 2000.0,
         dt: float = 0.125, amp: float = 0.0, per_a: float = 16) -> None:
    lx, ly = n_a * A, n_rows * 4.0 * np.pi
    nx = int(round(n_a * per_a))
    ny = int(round(ly / (A / per_a) / 2)) * 2
    x = (np.arange(nx) * lx / nx - lx / 2)[:, None]
    y = (np.arange(ny) * ly / ny - ly / 2)[None, :]
    env = 0.5 * (1.0 - np.tanh((np.sqrt(x**2 + y**2) - radius) / 3.0))
    psi = psi0 + env * tri_pattern(x, y, amp or tri_amplitude(eps, psi0))
    k2 = wavenumbers(nx, ny, lx, ly)
    smooth = np.exp(-0.5 * k2 * (A / 1.5) ** 2)
    ts, frac, fs, frames = [], [], [], []
    steps = int(round(t_end / dt))
    every = max(steps // 100, 1)
    for step in range(steps + 1):
        if step % every == 0:
            local = np.fft.irfft2(smooth * np.fft.rfft2((psi - psi0) ** 2), s=psi.shape)
            ts.append(step * dt)
            frac.append(float((local > 0.5 * local.max()).mean()) if local.max() > 1e-4 else 0.0)
            fs.append(energy_density(psi, k2, eps))
            if step % (every * 10) == 0:
                frames.append(psi.astype(np.float32))
        psi = relax(psi, k2, eps, dt, 1)
    tag = ("" if dt == 0.125 else f"_dt{dt:g}") + ("" if per_a == 16 else f"_n{per_a:g}")
    out = ROOT / f"results/crystallization/pfc/reference_grow_eps{eps:g}_psi{psi0:g}_r{radius:g}{tag}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, t=np.array(ts), fraction=np.array(frac), f=np.array(fs), frames=np.array(frames),
                        lx=lx, ly=ly)
    print("wrote", out, "final fraction", frac[-1], "f", fs[0], "->", fs[-1])


if __name__ == "__main__":
    main()
