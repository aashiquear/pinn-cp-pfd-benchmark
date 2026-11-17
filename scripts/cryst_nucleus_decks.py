#!/usr/bin/env python3
"""Write the standalone critical-nucleus decks of cases A3 and A4 from one base deck.

The base deck `decks/crystallization/single_field/a3_ni_nucleus_x0.200_fem.i` is hand written.
The other undercoolings and the other (well, interpolation) pairs differ from it only in
numbers (w, eps^2, Delta_f, the droplet of the initial guess, the expected values) and in the
free-energy expression; this script substitutes them from the reference table
(`scripts/cryst_reference.py`). Every deck it writes is a complete input: `pinnacle run <deck>`.

Usage: python3 scripts/cryst_nucleus_decks.py
"""
from __future__ import annotations

import csv
import re
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "decks/crystallization/single_field/a3_ni_nucleus_x0.200_fem.i"
TABLE = ROOT / "results/crystallization/reference/ni_nucleus_3d/nucleus_table.csv"
K_B = 1.380649e-23

EXPRESSIONS = {
    ("quartic", "standard"): "w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df",
    ("gl_fcc", "standard"): "w*phi**2*(1-phi**2)**2/6 + phi**3*(10-15*phi+6*phi**2)*df",
    ("quartic", "gl_bcc"): "w*phi**2*(1-phi)**2/4 + phi**3*(4-3*phi)*df",
    ("gl_fcc", "gl_fcc"): "w*phi**2*(1-phi**2)**2/6 + phi**4*(3-2*phi**2)*df",
}
A4_UNDERCOOLINGS = {0.15, 0.2, 0.25, 0.3}


def main() -> None:
    base = BASE.read_text()
    rows = [r for r in csv.DictReader(TABLE.open()) if r["temperature_scaling"] == "proportional"]
    profiles = json.loads((TABLE.parent / "nucleus_profiles.json").read_text())
    written = []
    for r in rows:
        pair = (r["well"], r["interpolation"])
        x = float(r["reduced_undercooling"])
        standard = pair == ("quartic", "standard")
        if not standard and x not in A4_UNDERCOOLINGS:
            continue
        w = float(r["w_J_m3"]) * 1e-9  # GPa = aJ/nm^3
        eps2 = float(r["eps2_J_m"]) * 1e9  # aJ/nm
        df = float(r["delta_f_J_m3"]) * 1e-9
        temp = float(r["T_K"])
        kt = K_B * temp * 1e18
        w_star = float(r["W_star_J"]) * 1e18
        r_e = float(r["R_equimolar_m"]) * 1e9
        phi0 = float(r["phi_center"])
        delta = float(r["delta_s_m"]) * 1e9
        case = "a3" if standard else "a4"
        tag = f"x{x:.3f}" if standard else f"{pair[0]}_{pair[1]}_x{x:.3f}"
        name = f"{case}_ni_nucleus_{tag}_fem"
        text = base
        text = text.replace("a3_ni_nucleus_x0.200_fem", name)
        # header
        text = text.replace("Part IV, case A3 (FE route)", f"Part IV, case {case.upper()} (FE route)")
        text = text.replace("reduced undercooling Delta T / T_f = 0.20 (T = 1382.64 K).",
                            f"reduced undercooling Delta T / T_f = {x:g} (T = {temp:.2f} K).")
        if not standard:
            text = text.replace(
                "#     g = phi^2 (1 - phi)^2 / 4,   p = phi^3 (10 - 15 phi + 6 phi^2)",
                f"#     w g + p Delta_f = {EXPRESSIONS[pair]}\n"
                f"#     (well '{pair[0]}', interpolation '{pair[1]}' of the review's Table 1; Fig. 5)")
        text = re.sub(r"#   derived at T = 0\.8 T_f.*?\n#     Delta_f = [^\n]*\n",
                      f"#   derived at T = {1 - x:g} T_f (eps^2 and w proportional to T; gamma and the\n"
                      f"#   10-90 % width matched for this well):\n"
                      f"#     eps^2 = {eps2:.6f} aJ/nm,  w = {w:.6f} GPa\n"
                      f"#     Delta_f = -H_f (T_f - T)/(T_f v_m) = {df:.7f} GPa\n", text, flags=re.S)
        text = re.sub(r"#   W\* = 0\.91370 aJ.*?\n#   In results",
                      f"#   W* = {w_star:.5f} aJ = {w_star / kt:.3f} kT (kT = {kt:.7f} aJ); phi(0) = {phi0:.5f};\n"
                      f"#   equimolar radius {r_e:.4f} nm.\n#   In results", text, flags=re.S)
        if not standard:
            text = re.sub(r"#   derived at T_f: eps\^2 = [^\n]*\n", "", text)
        # numbers
        # Initial guess of the Newton solve: a three-parameter droplet (amplitude, radius,
        # width) fitted to the reference profile. The nucleus is a saddle point, so Newton
        # needs a start close to it; the equimolar radius alone is not close enough for the
        # GL-fcc well, whose planar profile 1/sqrt(1 + exp(2 x/lam)) is not symmetric.
        prof = profiles[f"{pair[0]}:{pair[1]}:proportional:{x:g}"]
        rr = np.asarray(prof["r_m"]) * 1e9
        uu = np.asarray(prof["phi"])
        if pair[0] == "gl_fcc":
            shape = lambda q: q[0] / np.sqrt(1.0 + np.exp(2.0 * (rr - q[1]) / q[2]))  # noqa: E731
            a, rad, wid = least_squares(lambda q: shape(q) - uu, [phi0, r_e, (3 * eps2 / w) ** 0.5]).x
            guess = f"{a:.4f}/sqrt(1+exp(2*(x-{rad:.4f})/{wid:.5f}))"
        else:
            shape = lambda q: q[0] * 0.5 * (1.0 - np.tanh((rr - q[1]) / (2.0 * q[2])))  # noqa: E731
            a, rad, wid = least_squares(lambda q: shape(q) - uu, [phi0, r_e, delta]).x
            guess = f"{a:.4f}*0.5*(1-tanh((x-{rad:.4f})/(2*{wid:.5f})))"
        text = text.replace("expression = 0.9127*0.5*(1-tanh((x-0.9825)/(2*0.21791)))",
                            f"expression = {guess}")
        text = text.replace("expression = 'w*phi**2*(1-phi)**2/4 + phi**3*(10-15*phi+6*phi**2)*df'",
                            f"expression = '{EXPRESSIONS[pair]}'")
        text = text.replace("w  = 12.740387", f"w  = {w:.6f}")
        text = text.replace("df = -0.4923638", f"df = {df:.7f}")
        text = text.replace("0.329871", f"{eps2:.6f}")
        # domain: at least ten equimolar radii, same element size
        length = max(12.0, round(10.0 * r_e))
        text = text.replace("domain            = '0.0 12.0'      # nm; ten nucleus radii",
                            f"domain            = '0.0 {length:.1f}'      # nm; ten nucleus radii")
        text = text.replace("divisions         = 2400 ", f"divisions         = {int(length / 0.005)} ")
        out = BASE.parent / f"{name}.i"
        if out != BASE:
            out.write_text(text)
        written.append(out.name)
    print("\n".join(written))



def network_decks() -> list[str]:
    """PINN and PINN+FEM decks of case A3 for every undercooling, from the x0.200 bases.

    Their initial guess is the classical droplet (amplitude 1, R = 2 gamma(T) / |Delta_f|,
    planar width), which needs no reference solution.
    """
    rows = [r for r in csv.DictReader(TABLE.open()) if r["temperature_scaling"] == "proportional"
            and r["well"] == "quartic" and r["interpolation"] == "standard"]
    base_row = next(r for r in rows if abs(float(r["reduced_undercooling"]) - 0.2) < 1e-9)
    written = []
    for route in ("pinn", "hybrid"):
        base_path = BASE.parent / f"a3_ni_nucleus_x0.200_{route}.i"
        base = base_path.read_text()
        for r in rows:
            x = float(r["reduced_undercooling"])
            if abs(x - 0.2) < 1e-9:
                continue
            w = float(r["w_J_m3"]) * 1e-9
            eps2 = float(r["eps2_J_m"]) * 1e9
            df = float(r["delta_f_J_m3"]) * 1e-9
            temp = float(r["T_K"])
            kt = K_B * temp * 1e18
            w_star = float(r["W_star_J"]) * 1e18
            r_e = float(r["R_equimolar_m"]) * 1e9
            r_cnt = float(r["R_cnt_m"]) * 1e9
            text = base.replace(f"a3_ni_nucleus_x0.200_{route}", f"a3_ni_nucleus_x{x:.3f}_{route}")
            text = text.replace("reduced undercooling Delta T / T_f = 0.20 (T = 1382.64 K).",
                                f"reduced undercooling Delta T / T_f = {x:g} (T = {temp:.2f} K).")
            text = re.sub(r"#   derived at T = 0\.8 T_f.*?\n#     Delta_f = [^\n]*\n",
                          f"#   derived at T = {1 - x:g} T_f (eps^2 and w proportional to T):\n"
                          f"#     eps^2 = {eps2:.6f} aJ/nm,  w = {w:.6f} GPa\n"
                          f"#     Delta_f = -H_f (T_f - T)/(T_f v_m) = {df:.7f} GPa\n", text, flags=re.S)
            text = re.sub(r"#   W\* = 0\.91370 aJ.*?\n#   In results",
                          f"#   W* = {w_star:.5f} aJ = {w_star / kt:.3f} kT (kT = {kt:.7f} aJ); "
                          f"phi(0) = {float(r['phi_center']):.5f};\n"
                          f"#   equimolar radius {r_e:.4f} nm.\n#   In results", text, flags=re.S)
            text = text.replace("0.981 nm", f"{r_cnt:.3f} nm").replace("= 0.9814 nm", f"= {r_cnt:.4f} nm")
            text = text.replace("expression = 0.5*(1-tanh((x-0.9814)/(2*0.22756)))",
                                f"expression = 0.5*(1-tanh((x-{r_cnt:.4f})/(2*0.22756)))")
            text = text.replace("w  = 12.740387", f"w  = {w:.6f}").replace("df = -0.4923638", f"df = {df:.7f}")
            text = text.replace("0.329871", f"{eps2:.6f}")
            if route == "pinn":
                length = max(6.0, round(5.0 * r_e))
                text = text.replace("domain            = '0.0 6.0'       # nm; five nucleus radii",
                                    f"domain            = '0.0 {length:.1f}'       # nm; five nucleus radii")
                text = text.replace("divisions         = 1200 ", f"divisions         = {int(length / 0.005)} ")
            else:
                length = max(12.0, round(10.0 * r_e))
                text = text.replace("domain            = '0.0 12.0'      # nm; ten nucleus radii",
                                    f"domain            = '0.0 {length:.1f}'      # nm; ten nucleus radii")
                text = text.replace("divisions         = 2400 ", f"divisions         = {int(length / 0.005)} ")
            assert "0.4923638" not in text and "12.740387" not in text, (route, x)
            (BASE.parent / f"a3_ni_nucleus_x{x:.3f}_{route}.i").write_text(text)
            written.append(f"a3_ni_nucleus_x{x:.3f}_{route}.i")
    del base_row
    return written


if __name__ == "__main__":
    main()
    print("\n".join(network_decks()))
