#!/usr/bin/env python3
"""Independent reference for case F1: crystal nucleation in the hard-sphere liquid.

Phase-field theory of (granasy 2003 hardsphere): structural order parameter m (1 liquid,
0 crystal) and volume fraction phi = pi sigma^3 rho / 6. Units: kT = 1, sigma = 1.

    F = int { (b/2) |grad m|^2 + w g(m) + [1 - p(m)] f_S(phi) + p(m) f_L(phi) } dr
    g = m^2 (1 - m)^2 / 4,  p = m^3 (10 - 15 m + 6 m^2)

The nucleus is an extremum of the grand potential Omega = F - mu_inf N:
    b (m'' + 2 m'/r) = w g'(m) + p'(m) [f_L(phi) - f_S(phi)]
    [1 - p(m)] mu_S(phi) + p(m) mu_L(phi) = mu_inf            (phi is slaved to m)
    W* = int 4 pi r^2 { (b/2) m'^2 + Delta omega(m) } dr

Bulk free energies: integrated from the equations of state of Hall (J. Chem. Phys. 57, 2252,
1972): the rational form for the fluid (coefficients as reproduced in Eckert's thesis, Univ.
Bayreuth, Eq. 1.11, and Schinagl et al., Condens. Matter Phys. 2, 313, 1999, Eq. 11) and the
polynomial fit for the fcc solid (Hall's Eq. 13, as reproduced in arXiv cond-mat/0703434
Eq. 27 and cond-mat/0510445 Eq. 5.17). The paper cites Hall without printing the forms; this
pair is `derived`: it reproduces the paper's solid densities 0.5796, 0.5875, 0.5948 (equal
pressure with the liquid at 0.5207, 0.5277, 0.5343) to four digits, the rational solid form
does not. b and w follow from the orientation-averaged interface free energy 0.6133 kT/sigma^2
and 10-90 % thickness 5.410 sigma (`reported_text`), Eq. 10-11 of the paper.

Usage: python3 scripts/cryst_hard_sphere_reference.py            (prints and writes the table)
Writes results/crystallization/hard_sphere/reference.csv and reference_profiles.npz.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, fsolve

ROOT = Path(__file__).resolve().parents[1]
ETA_CP = np.pi * np.sqrt(2.0) / 6.0
PHI_L_COEX = 0.494
GAMMA, D1090 = 0.6133, 5.410


def z_fluid(y):
    return (1 + y + y**2 - 0.67825 * y**3 - y**4 - 0.5 * y**5 - 1.7 * y**6) / (
        1 - 3 * y + 3 * y**2 - 1.04305 * y**3)


def z_solid(y):
    """Hall's polynomial fit for the fcc solid (his Eq. 13), beta = 4 (1 - phi/phi_cp)."""
    beta = 4.0 * (1.0 - y / ETA_CP)
    return ((12.0 - 3.0 * beta) / beta + 2.557696 + 0.1253077 * beta + 0.1762393 * beta**2
            - 1.053308 * beta**3 + 2.818621 * beta**4 - 2.921934 * beta**5 + 1.118413 * beta**6)


def z_solid_rational(y):
    """Hall's rational solid form (not the one the paper used; kept for the comparison)."""
    xi = ETA_CP - y
    return (1 + y + y**2 - 0.67825 * y**3 - y**4 - 0.5 * y**5
            - 6.028 * np.exp(xi * (7.9 - 3.9 * xi)) * y**6) / (1 - 3 * y + 3 * y**2 - 1.04305 * y**3)


class Bulk:
    """f(phi), mu(phi) = df/drho, pressure, by integrating d(a)/d(rho) = Z/rho (a per particle)."""

    def __init__(self, z, lo: float, hi: float, n: int = 20001) -> None:
        self.phi = np.linspace(lo, hi, n)
        zz = z(self.phi)
        a = cumulative_trapezoid(zz / self.phi, self.phi, initial=0.0)  # d rho/rho = d phi/phi
        self._a = CubicSpline(self.phi, a)
        self._z = z
        self.shift = 0.0

    def a(self, phi):
        return self._a(phi) + self.shift

    def mu(self, phi):  # chemical potential per particle
        return self.a(phi) + self._z(phi)

    def dmu(self, phi):  # d mu / d phi
        return self._a(phi, 1) + (self._z(phi + 1e-6) - self._z(phi - 1e-6)) / 2e-6

    def f(self, phi):  # per volume
        return 6.0 * phi / np.pi * self.a(phi)

    def pressure(self, phi):
        return 6.0 * phi / np.pi * self._z(phi)


LIQ = Bulk(z_fluid, 0.30, 0.70)
SOL = Bulk(z_solid, 0.30, 0.70)
# Coexistence: the liquid at the paper's 0.494; the solid at equal pressure (0.5466 with these
# two equations of state; the paper prints 0.546), and the solid constant from equal chemical
# potential there. Fixing mu at 0.546 instead leaves a pressure mismatch and gives W* = 22.
PHI_S_COEX = float(brentq(lambda y: SOL.pressure(y) - LIQ.pressure(PHI_L_COEX), 0.50, 0.60))
SOL.shift = float(LIQ.mu(PHI_L_COEX) - SOL.mu(PHI_S_COEX))


def p_fun(m):
    return m**3 * (10 - 15 * m + 6 * m * m)


def dp_fun(m):
    return 30 * m * m * (1 - m) ** 2


def g_fun(m):
    return 0.25 * m * m * (1 - m) ** 2


def dg_fun(m):
    return 0.5 * m * (1 - m) * (1 - 2 * m)


def slaved_phi(m, mu_inf: float):
    """Solve [1 - p] mu_S(phi) + p mu_L(phi) = mu_inf for phi (monotonic in phi)."""
    m = np.atleast_1d(np.asarray(m, dtype=float))
    lo, hi = np.full(m.shape, 0.35), np.full(m.shape, 0.69)
    pm = p_fun(np.clip(m, 0.0, 1.0))
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        val = (1 - pm) * SOL.mu(mid) + pm * LIQ.mu(mid) - mu_inf
        lo = np.where(val < 0, mid, lo)
        hi = np.where(val < 0, hi, mid)
    return 0.5 * (lo + hi)


class Slaved:
    """Tables on a grid of m for one liquid state: phi(m), the bulk part of Delta omega and
    its derivative p'(m) [f_L(phi) - f_S(phi)] (the density drops out by the envelope theorem)."""

    def __init__(self, phi_inf: float, n: int = 4001) -> None:
        self.mu_inf, self.p_inf = float(LIQ.mu(phi_inf)), float(LIQ.pressure(phi_inf))
        m = np.linspace(0.0, 1.0, n)
        phi = slaved_phi(m, self.mu_inf)
        pm = p_fun(m)
        excess = (1 - pm) * SOL.f(phi) + pm * LIQ.f(phi) - self.mu_inf * 6.0 * phi / np.pi + self.p_inf
        self.phi = CubicSpline(m, phi)
        self.excess = CubicSpline(m, excess)
        self.force = CubicSpline(m, dp_fun(m) * (LIQ.f(phi) - SOL.f(phi)))


def calibrate() -> tuple[float, float]:
    """b and w from the planar interface at coexistence (Eq. 10-11 of the paper)."""
    table = Slaved(PHI_L_COEX)
    m = np.linspace(0.0, 1.0, 200001)
    inner = (m >= 0.1) & (m <= 0.9)

    def equations(x):
        b, w = np.exp(x)
        excess = np.maximum(w * g_fun(m) + table.excess(m), 0.0)
        gam = np.sqrt(2 * b) * np.trapezoid(np.sqrt(excess), m)
        d = np.sqrt(b / 2) * np.trapezoid(excess[inner] ** -0.5, m[inner])
        return [gam / GAMMA - 1, d / D1090 - 1]

    b, w = np.exp(fsolve(equations, np.log([3.0, 6.0]), xtol=1e-12))
    return float(b), float(w)


def nucleus(phi_inf: float, b: float, w: float, r_max: float = 80.0):
    """Shooting on m(0); returns W*, the profile and the slaved density."""
    table = Slaved(phi_inf)

    def force(m):
        return (w * dg_fun(m) + float(table.force(min(max(m, 0.0), 1.0)))) / b

    def rhs(r, y):
        return [y[1], force(y[0]) - 2 * y[1] / r]

    def over(r, y):  # passed the liquid value
        return y[0] - 1.0

    def back(r, y):  # turned back
        return y[1]

    over.terminal = back.terminal = True
    back.direction = -1

    def shoot(m0):
        f0, r0 = force(m0), 1e-3
        return solve_ivp(rhs, [r0, r_max], [m0 + f0 * r0**2 / 6, f0 * r0 / 3], events=[over, back],
                         rtol=1e-10, atol=1e-13, dense_output=True)

    lo, hi = 0.0, 0.995  # m(0): too small overshoots m = 1, too large turns back
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if shoot(mid).t_events[0].size:
            lo = mid
        else:
            hi = mid
    sol = shoot(0.5 * (lo + hi))
    r = np.linspace(1e-3, sol.t[-1], 6000)
    m, dm = sol.sol(r)
    last = int(np.argmax(m))
    r, m, dm = r[: last + 1], np.clip(m[: last + 1], 0.0, 1.0), dm[: last + 1]
    dens = 0.5 * b * dm**2 + w * g_fun(m) + table.excess(m)
    w_star = float(np.trapezoid(4 * np.pi * r**2 * dens, r))
    return w_star, r, m, table.phi(m)


# Interface data of the paper (p. 3; `reported_text`): free energy in kT/sigma^2 and 10-90 %
# thickness in sigma for three orientations, and their average, which the paper uses for the
# headline barrier.
INTERFACES = {"average": (0.6133, 5.410), "111": (0.58, 5.94), "110": (0.64, 4.75), "100": (0.62, 5.54)}
VOLUME_FRACTIONS = (0.5207, 0.5240, 0.5277, 0.5310, 0.5343, 0.5400, 0.5500)


def main() -> None:
    global GAMMA, D1090
    out = ROOT / "results/crystallization/hard_sphere"
    out.mkdir(parents=True, exist_ok=True)
    print(f"coexistence: phi_L {PHI_L_COEX}, phi_S {PHI_S_COEX:.4f}, pressure {LIQ.pressure(PHI_L_COEX):.3f} "
          "(paper: 0.494, 0.546, 11.81)")
    rows, profiles = [], {}
    for interface, (GAMMA, D1090) in INTERFACES.items():
        b, w = calibrate()
        print(f"interface {interface}: gamma {GAMMA}, d10-90 {D1090} -> b = {b:.5f} kT/sigma, w = {w:.5f} kT/sigma^3")
        profiles[f"b_{interface}"], profiles[f"w_{interface}"] = b, w
        for phi_inf in VOLUME_FRACTIONS:
            mu_inf, p_inf = float(LIQ.mu(phi_inf)), float(LIQ.pressure(phi_inf))
            grid = np.linspace(0.45, 0.69, 4001)
            gain = LIQ.f(phi_inf) + mu_inf * 6 / np.pi * (grid - phi_inf) - SOL.f(grid)
            d_omega = float(gain.max())  # driving force per volume of the classical droplet
            w_star, r, m, phi = nucleus(phi_inf, b, w)
            rows.append({"interface": interface, "gamma": GAMMA, "d1090": D1090, "b": b, "w": w,
                         "phi_inf": phi_inf, "mu_inf": mu_inf, "pressure_inf": p_inf, "driving_force": d_omega,
                         "W_cnt": 16 * np.pi / 3 * GAMMA**3 / d_omega**2, "W_pft": w_star,
                         "m_center": float(m[0]), "phi_center": float(phi[0]),
                         "phi_solid_equal_mu": float(slaved_phi(0.0, mu_inf)[0]),
                         "phi_solid_equal_pressure": float(brentq(lambda y: SOL.pressure(y) - p_inf, 0.5, 0.69))})
            profiles[f"{interface}_{phi_inf:.4f}"] = np.vstack([r, m, phi])
            print("  " + " ".join(f"{k}={v:.5g}" for k, v in rows[-1].items() if k in
                                  ("phi_inf", "driving_force", "W_cnt", "W_pft", "m_center", "phi_center",
                                   "phi_solid_equal_pressure")))
    with (out / "reference.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    np.savez(out / "reference_profiles.npz", **profiles)
    print("paper: W*/kT = 40.63 at phi_inf = 0.5207 (average interface); "
          "phi_S = 0.5796, 0.5875, 0.5948 at 0.5207, 0.5277, 0.5343")


if __name__ == "__main__":
    main()
