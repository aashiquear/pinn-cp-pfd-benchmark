#!/usr/bin/env python3
"""Thermodynamics and dimensionless groups of the Ag-Cu regular-solution eutectic model (case F2).

Model: (lewis 2004 eutectic), (ratkai 2017 nanoscale). Free energies per mole, liquid reference:

    G_L = Omega_L c (1 - c) + R T [c ln c + (1 - c) ln(1 - c)]
    G_S = (1 - c) L_A v (T - T_A)/T_A + c L_B v (T - T_B)/T_B + Omega_S c (1 - c) + R T [...]
    Omega_L = 16,000 - 4 T,  Omega_S = 26,000 - 4 T  J/mol            (c = mole fraction of Cu)

Data (`reported_text`): L_A = 1097 J/cm^3, L_B = 1728 J/cm^3, v = 8.7 cm^3/mol, T_A = 1234 K,
T_B = 1358 K (lewis 2004 eutectic); gamma_A = 172, gamma_B = 227 mJ/m^2, 10-90 % interface
thickness 1 nm, eps_c^2 = eps_phi^2 / 4, D_L = 1.4e-9 m^2/s, M_phi = 0.25 m^3/(J s)
(ratkai 2017 nanoscale).

Dimensionless form used by the decks: length xi = 1 nm, time xi^2/D_L, energy density R T/v.

Usage: python3 scripts/cryst_eutectic_thermo.py [T ...]    (prints the phase diagram points and
the dimensionless groups at the given temperatures; default 900 K)
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.optimize import fsolve

R = 8.314462
V_M = 8.7e-6  # m^3/mol
T_A, T_B = 1234.0, 1358.0
H_A, H_B = 1097.0e6 * V_M, 1728.0e6 * V_M  # J/mol: 9543.9, 15033.6
GAMMA_A, GAMMA_B = 0.172, 0.227  # J/m^2
D1090 = 1.0e-9  # m
D_L = 1.4e-9  # m^2/s
M_PHI = 0.25  # m^3/(J s)
XI = 1.0e-9  # m


def omega(t: float, phase: str) -> float:
    return (26000.0 if phase == "S" else 16000.0) - 4.0 * t


def g(c, t: float, phase: str):
    base = 0.0 if phase == "L" else (1 - c) * H_A * (t - T_A) / T_A + c * H_B * (t - T_B) / T_B
    return base + omega(t, phase) * c * (1 - c) + R * t * (c * np.log(c) + (1 - c) * np.log(1 - c))


def dg(c, t: float, phase: str):
    base = 0.0 if phase == "L" else -H_A * (t - T_A) / T_A + H_B * (t - T_B) / T_B
    return base + omega(t, phase) * (1 - 2 * c) + R * t * np.log(c / (1 - c))


def eutectic() -> tuple[float, float, float, float]:
    """(c_alpha, c_beta, c_E, T_E) of the three-phase equilibrium."""

    def eqs(x):
        ca, cb, cl, t = x
        mu = dg(cl, t, "L")
        return [dg(ca, t, "S") - mu, dg(cb, t, "S") - mu,
                g(ca, t, "S") - ca * mu - (g(cl, t, "L") - cl * mu),
                g(cb, t, "S") - cb * mu - (g(cl, t, "L") - cl * mu)]

    return tuple(fsolve(eqs, [0.15, 0.85, 0.35, 1055.0], xtol=1e-13))  # type: ignore[return-value]


def solid_solvus(t: float) -> tuple[float, float]:
    """Compositions of the two solid solutions in equilibrium with each other at T."""

    def eqs(x):
        ca, cb = x
        return [dg(ca, t, "S") - dg(cb, t, "S"),
                g(ca, t, "S") - ca * dg(ca, t, "S") - (g(cb, t, "S") - cb * dg(cb, t, "S"))]

    ca, cb = fsolve(eqs, [0.05, 0.95], xtol=1e-13)
    return float(ca), float(cb)


def groups(t: float) -> dict[str, float]:
    """Dimensionless groups of the decks at temperature T."""
    scale = R * t / V_M  # J/m^3
    delta_a = D1090 / (2.0 * np.log(9.0) * np.sqrt(2.0))  # 10-90 % width = 2 ln 9 sqrt(2) delta
    eps2_t = 6.0 * np.sqrt(2.0) * GAMMA_A * delta_a * t / T_A  # J/m
    delta_b = delta_a * (GAMMA_A / T_A) / (GAMMA_B / T_B)  # the same eps^2 for both components
    ca, cb = solid_solvus(t)
    kappa_phi = eps2_t / scale / XI**2
    return {
        "T": t,
        "solid_a": H_A * (t - T_A) / (T_A * R * t),
        "solid_b": H_B * (t - T_B) / (T_B * R * t),
        "omega_solid": omega(t, "S") / (R * t),
        "omega_liquid": omega(t, "L") / (R * t),
        "kappa_phi": kappa_phi,
        "kappa_c": kappa_phi / 4.0,
        "well_a": 6.0 * np.sqrt(2.0) * GAMMA_A / delta_a * t / T_A / scale,
        "well_b": 6.0 * np.sqrt(2.0) * GAMMA_B / delta_b * t / T_B / scale,
        "m_phi": M_PHI * scale * XI**2 / D_L,
        "c_alpha": ca,
        "c_beta": cb,
        "time_unit_s": XI**2 / D_L,
        "energy_unit_J_m3": scale,
    }


def main() -> None:
    ca, cb, ce, te = eutectic()
    print(f"eutectic: T_E = {te:.2f} K, c_E = {ce:.4f}, c_alpha = {ca:.4f}, c_beta = {cb:.4f} "
          "(papers: 1055.5 K, 0.35)")
    for t in [float(a) for a in sys.argv[1:]] or [900.0]:
        print(" ".join(f"{k}={v:.6g}" for k, v in groups(t).items()))


if __name__ == "__main__":
    main()
