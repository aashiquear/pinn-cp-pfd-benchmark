#!/usr/bin/env python3
"""Independent reference solutions for the single-field phase-field model of crystallization.

Part IV of the study (phase-field crystallization benchmark). Nothing in this file calls
PINNACLE: it is the independently verified numerical reference (accuracy tier 2 of
PROJECT_GUIDELINES.md) that the PINNACLE decks are compared with.

Model (granasy 2019 nucleation, Eq. 12, isotropic, s = 1), phi = 0 liquid, phi = 1 crystal:

    F = int { eps^2/2 |grad phi|^2 + w g(phi) + p(phi) Delta_f } dV

Closed forms for the quartic well g = phi^2 (1 - phi)^2 / 4 (derived in the report, section 2):

    planar profile      phi(x) = 1 / (1 + exp(x / delta)),   delta = sqrt(2 eps^2 / w)
    interface energy    gamma  = eps^2 / (6 delta) = (eps / 6) sqrt(w / 2)
    10-90 % thickness   d      = 2 ln(9) delta
    front velocity      v      = -6 M Delta_f delta      (p = phi^2 (3 - 2 phi) only)

For a general well the same quantities follow from the first integral
eps^2/2 phi'^2 = w g(phi):  gamma = eps sqrt(2 w) I_g with I_g = int_0^1 sqrt(g) dphi and
d = eps / sqrt(2 w) int_0.1^0.9 dphi / sqrt(g).

Critical nucleus (Eq. 15), spherical symmetry (dim = 3) or circular (dim = 2):

    phi'' + (dim - 1)/r phi' = (w g' + p' Delta_f) / eps^2,   phi'(0) = 0,  phi(inf) = 0

solved by shooting on phi(0) in the scaled variables r/delta_s, f/w, with
delta_s = sqrt(2 eps^2 / w), where it reads phi'' + (dim-1)/r phi' = 2 (g' + D p'),
D = Delta_f / w, and W* = S_dim delta_s^dim w int r^(dim-1) [phi'^2/4 + g + D p] dr.

Usage:
    python3 scripts/cryst_reference.py <deck.i>

The deck (decks/crystallization/reference/*.i) names the material, the functions and the
undercoolings; the script writes CSV tables (and the radial profiles) to its out_dir.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

K_B = 1.380649e-23  # J/K (exact, SI 2019)

# ---------------------------------------------------------------------------------------
# double wells g(phi) and interpolation functions p(phi) of the review's Table 1
# each entry: (function, derivative)
# ---------------------------------------------------------------------------------------
WELLS = {
    # standard PF, GL bcc, dynamic tanh: quartic
    "quartic": (lambda u: 0.25 * u**2 * (1 - u) ** 2, lambda u: 0.5 * u * (1 - u) * (1 - 2 * u)),
    # GL fcc: g = phi^2 (1 - phi^2)^2 / 6
    "gl_fcc": (
        lambda u: u**2 * (1 - u**2) ** 2 / 6.0,
        lambda u: (2 * u * (1 - u**2) ** 2 - 4 * u**3 * (1 - u**2)) / 6.0,
    ),
}
INTERPOLATIONS = {
    "standard": (lambda u: u**3 * (10 - 15 * u + 6 * u**2), lambda u: 30 * u**2 * (1 - u) ** 2),
    "gl_bcc": (lambda u: u**3 * (4 - 3 * u), lambda u: 12 * u**2 * (1 - u)),
    "gl_fcc": (lambda u: u**4 * (3 - 2 * u**2), lambda u: 12 * u**3 * (1 - u**2)),
    "tanh": (lambda u: u**2 * (3 - 2 * u), lambda u: 6 * u * (1 - u)),
    "quadratic": (lambda u: u**2, lambda u: 2 * u),
}


@dataclass
class Model:
    """eps^2 [J/m] and w [J/m^3] at the melting point, fixed by gamma and the 10-90 % width."""

    eps2: float
    w: float
    well: str
    interp: str

    @property
    def delta_s(self) -> float:
        return math.sqrt(2.0 * self.eps2 / self.w)


def well_integrals(well: str) -> tuple[float, float]:
    """I_g = int_0^1 sqrt(g), J_g = int_0.1^0.9 dphi / sqrt(g)."""
    g = WELLS[well][0]
    i_g = quad(lambda u: math.sqrt(max(g(u), 0.0)), 0.0, 1.0, epsabs=1e-13, epsrel=1e-13)[0]
    j_g = quad(lambda u: 1.0 / math.sqrt(g(u)), 0.1, 0.9, epsabs=1e-13, epsrel=1e-13)[0]
    return i_g, j_g


def calibrate(gamma: float, d1090: float, well: str, interp: str) -> Model:
    """Model parameters from the interface free energy and the 10-90 % thickness.

    gamma = eps sqrt(2 w) I_g and d = eps / sqrt(2 w) J_g give
    eps^2 = gamma d / (I_g J_g),  w = gamma J_g / (2 d I_g).
    """
    i_g, j_g = well_integrals(well)
    return Model(eps2=gamma * d1090 / (i_g * j_g), w=gamma * j_g / (2.0 * d1090 * i_g),
                 well=well, interp=interp)


# ---------------------------------------------------------------------------------------
# critical nucleus by shooting
# ---------------------------------------------------------------------------------------
def _shoot(phi0: float, big_d: float, well: str, interp: str, dim: int, r_max: float):
    dg = WELLS[well][1]
    dp = INTERPOLATIONS[interp][1]

    def rhs(r, y):
        u, v = y
        src = 2.0 * (dg(u) + big_d * dp(u))
        if r < 1e-12:
            return [v, src / dim]  # limit of phi'' + (dim-1)/r phi' at r = 0
        return [v, src - (dim - 1) / r * v]

    def under(r, y):  # the profile turns back before reaching the liquid
        return y[1]

    under.terminal = True
    under.direction = 1.0

    def over(r, y):  # the profile crosses phi = 0
        return y[0]

    over.terminal = True
    over.direction = -1.0
    return solve_ivp(rhs, (0.0, r_max), [phi0, 0.0], events=(under, over), method="DOP853",
                     rtol=1e-11, atol=1e-13, dense_output=True, first_step=1e-6)


def critical_nucleus(big_d: float, well: str = "quartic", interp: str = "standard", dim: int = 3,
                     n_out: int = 4001):
    """Scaled critical nucleus for D = Delta_f / w < 0.

    Returns r (in delta_s), phi(r), the scaled barrier W* / (S_dim delta_s^dim w) and the
    center value phi(0).
    """
    g = WELLS[well][0]
    p = INTERPOLATIONS[interp][0]
    # lower bracket: the nonzero root of g + D p = 0 (energy of the start must be below 0)
    us = np.linspace(1e-6, 1.0, 200001)
    e = g(us) + big_d * p(us)
    neg = np.where(e < 0.0)[0]
    if len(neg) == 0:
        raise ValueError("no driving force: g + D p >= 0 everywhere")
    lo = us[neg[0]]
    # upper bracket: the solid minimum of g + D p (may sit slightly off 1 for some p)
    hi = us[np.argmin(e)]
    hi = hi - 1e-15 if hi < 1.0 else 1.0 - 1e-16
    r_cnt = (dim - 1) * well_integrals(well)[0] * 2.0 / abs(big_d)  # 2 gamma / |Delta f| scaled
    r_max = max(60.0, 6.0 * r_cnt)
    sol = None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        s = _shoot(mid, big_d, well, interp, dim, r_max)
        if len(s.t_events[1]) > 0:  # overshoot: started too high
            hi = mid
        else:  # undershoot (or ran to r_max without crossing): started too low
            lo = mid
        sol = s
        if hi - lo < 1e-15:
            break
    # take the undershooting side, cut where it turns back, extend with zeros
    s = _shoot(lo, big_d, well, interp, dim, r_max)
    r_end = s.t_events[0][0] if len(s.t_events[0]) > 0 else s.t[-1]
    r = np.linspace(0.0, r_end, n_out)
    y = s.sol(r)
    u, v = y[0], y[1]
    integrand = r ** (dim - 1) * (0.25 * v**2 + g(u) + big_d * p(u))
    w_scaled = np.trapezoid(integrand, r)
    del sol
    return r, u, w_scaled, lo


def equimolar_radius(r: np.ndarray, u: np.ndarray, dim: int) -> float:
    """Radius of the step profile of amplitude phi(0) with the same radial integral."""
    integral = np.trapezoid(r ** (dim - 1) * u, r)
    return (dim * integral / u[0]) ** (1.0 / dim)


def sphere_factor(dim: int) -> float:
    return {1: 2.0, 2: 2.0 * math.pi, 3: 4.0 * math.pi}[dim]


# ---------------------------------------------------------------------------------------
# deck-driven tables
# ---------------------------------------------------------------------------------------
def read_deck(path: Path) -> dict:
    """Minimal reader of the PINNACLE-style `[Block] key = value` study decks."""
    out: dict = {}
    block = None
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            name = line[1:-1].strip()
            block = None if name == "" else out.setdefault(name, {})
            continue
        if "=" in line and block is not None:
            k, v = (s.strip() for s in line.split("=", 1))
            block[k] = v.strip("'\"")
    return out


def nickel_tables(deck: dict, out_dir: Path) -> None:
    mat = deck["Material"]
    t_f = float(mat["T_f"])
    h_f = float(mat["H_f"])  # J/mol
    v_m = float(mat["molar_mass"]) / float(mat["density"])  # m^3/mol
    gamma_eq = float(mat["gamma_eq"])
    d1090 = float(mat["d_10_90"])
    study = deck["Study"]
    dim = int(study.get("dim", 3))
    reduced = [float(s) for s in study["reduced_undercoolings"].split()]
    pairs = [tuple(s.split(":")) for s in study["functions"].split()]
    scalings = study.get("temperature_scaling", "proportional").split()
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    profiles = {}
    for well, interp in pairs:
        model = calibrate(gamma_eq, d1090, well, interp)
        i_g, _ = well_integrals(well)
        for scaling in scalings:
            for x in reduced:
                temp = t_f * (1.0 - x)
                ratio = temp / t_f if scaling == "proportional" else 1.0
                w_t = model.w * ratio
                eps2_t = model.eps2 * ratio
                gamma_t = gamma_eq * ratio
                delta_f = -h_f * x / v_m  # Turnbull: H_f (T - T_f) / (T_f v_m)
                big_d = delta_f / w_t
                r, u, w_scaled, phi0 = critical_nucleus(big_d, well, interp, dim)
                w_star = sphere_factor(dim) * model.delta_s**dim * w_t * w_scaled
                r_e = equimolar_radius(r, u, dim) * model.delta_s
                # classical droplet with the same gamma(T)
                r_cnt = (dim - 1) * gamma_t / abs(delta_f)
                if dim == 3:
                    w_cnt = 16.0 * math.pi / 3.0 * gamma_t**3 / delta_f**2
                else:
                    w_cnt = math.pi * gamma_t**2 / abs(delta_f)  # per unit length
                kt = K_B * temp
                rows.append({
                    "well": well, "interpolation": interp, "temperature_scaling": scaling,
                    "dim": dim, "reduced_undercooling": x, "T_K": temp,
                    "delta_f_J_m3": delta_f, "w_J_m3": w_t, "eps2_J_m": eps2_t,
                    "delta_s_m": model.delta_s, "D": big_d, "phi_center": phi0,
                    "W_star_J": w_star, "W_star_over_kT": w_star / kt,
                    "R_equimolar_m": r_e, "W_cnt_J": w_cnt, "W_cnt_over_kT": w_cnt / kt,
                    "R_cnt_m": r_cnt,
                })
                profiles[f"{well}:{interp}:{scaling}:{x:g}"] = {
                    "r_m": (r[::20] * model.delta_s).tolist(), "phi": u[::20].tolist()}
    with (out_dir / "nucleus_table.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    (out_dir / "nucleus_profiles.json").write_text(json.dumps(profiles))
    base = calibrate(gamma_eq, d1090, "quartic", "standard")
    summary = {
        "eps2_J_per_m": base.eps2, "w_J_per_m3": base.w, "delta_s_m": base.delta_s,
        "gamma_check_J_per_m2": math.sqrt(base.eps2) / 6.0 * math.sqrt(base.w / 2.0),
        "d_10_90_check_m": 2.0 * math.log(9.0) * base.delta_s,
        "H_f_over_v_m_J_per_m3": h_f / v_m, "v_m_m3_per_mol": v_m,
    }
    (out_dir / "model_parameters.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    for row in rows:
        print(f"{row['well']:8s} {row['interpolation']:9s} {row['temperature_scaling']:12s} "
              f"x={row['reduced_undercooling']:.3f} W*/kT={row['W_star_over_kT']:8.2f} "
              f"CNT={row['W_cnt_over_kT']:8.2f} R_e={row['R_equimolar_m'] * 1e10:6.2f} A "
              f"phi0={row['phi_center']:.6f}")


# ---------------------------------------------------------------------------------------
# steady planar front (case A2): travelling wave phi(x - v t) by shooting on the speed
# ---------------------------------------------------------------------------------------
def front_speed(big_d: float, well: str = "quartic", interp: str = "standard") -> float:
    """Scaled front speed c of phi'' + c phi' = 2 (g' + D p'), phi(-inf) = 1, phi(+inf) = 0.

    Scaled variables: x / delta_s and t M w / 2, so the physical speed is
    v = c delta_s M w / 2. For p = phi^2 (3 - 2 phi) the closed form is c = -12 D.
    The orbit leaves phi = 1 along its unstable direction; a speed that is too small
    overshoots phi = 0, one that is too large turns back before reaching it.
    """
    dg = WELLS[well][1]
    dp = INTERPOLATIONS[interp][1]
    h = 1e-6
    s1 = 2.0 * ((dg(1.0) - dg(1.0 - h)) / h + big_d * (dp(1.0) - dp(1.0 - h)) / h)  # S'(1)

    def outcome(c: float) -> int:
        lam = 0.5 * (-c + math.sqrt(c * c + 4.0 * s1))
        a = 1e-8

        def rhs(x, y):
            return [y[1], -c * y[1] + 2.0 * (dg(y[0]) + big_d * dp(y[0]))]

        def cross(x, y):
            return y[0]

        cross.terminal = True
        cross.direction = -1.0

        def turn(x, y):
            return y[1]

        turn.terminal = True
        turn.direction = 1.0
        sol = solve_ivp(rhs, (0.0, 400.0), [1.0 - a, -a * lam], events=(cross, turn),
                        method="DOP853", rtol=1e-11, atol=1e-14)
        return -1 if len(sol.t_events[0]) > 0 else 1  # -1: overshoot (c too small)

    lo, hi = 0.0, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if outcome(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def front_table(deck: dict, out_dir: Path) -> None:
    mat = deck["Material"]
    t_f = float(mat["T_f"])
    v_m = float(mat["molar_mass"]) / float(mat["density"])
    study = deck["Study"]
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for x in (float(v) for v in study["reduced_undercoolings"].split()):
        for interp in study["interpolations"].split():
            model = calibrate(float(mat["gamma_eq"]), float(mat["d_10_90"]), "quartic", interp)
            ratio = 1.0 - x
            w_t = model.w * ratio
            delta_f = -float(mat["H_f"]) * x / v_m
            c = front_speed(delta_f / w_t, "quartic", interp)
            # nm per unit time for M = 1 in (nm, aJ, GPa) units
            v = c * model.delta_s * 1e9 * (w_t * 1e-9) / 2.0
            closed = -6.0 * delta_f * 1e-9 * model.delta_s * 1e9 if interp == "tanh" else float("nan")
            rows.append({"reduced_undercooling": x, "interpolation": interp, "D": delta_f / w_t,
                         "c_scaled": c, "v_nm_per_time": v, "v_closed_form": closed})
            print(f"x={x:.3f} {interp:9s} c={c:.8f} v={v:.6f} closed={closed:.6f}")
    with (out_dir / "front_table.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


def main() -> None:
    deck_path = Path(sys.argv[1]).resolve()
    deck = read_deck(deck_path)
    root = Path(__file__).resolve().parents[1]
    out_dir = root / deck["Study"]["out_dir"]
    if deck["Study"].get("kind", "nucleus") == "front":
        front_table(deck, out_dir)
    else:
        nickel_tables(deck, out_dir)


if __name__ == "__main__":
    main()
