#!/usr/bin/env python3
"""Write the parameter-series decks of the orientation-field model (Part IV, case C3).

One-at-a-time variations around the conditions of Fig. 4 of (granasy 2005 spherulites):
supersaturation S, anisotropy s0 of the interface free energy, orientational mobility chi_l.
Every deck is complete and runs alone (`pinnacle run <deck>`); this script only saves typing.
The base deck is decks/crystallization/orientation/c3_fig4c_s0.8_chi54_fem.i.

Usage: python3 scripts/cryst_orientation_decks.py        (writes decks/.../orientation/series/)
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "decks/crystallization/orientation/c3_fig4c_s0.8_chi54_fem.i"
OUT = ROOT / "decks/crystallization/orientation/series"
C_L, C_S = 0.466219, 0.399112

# (S, s0, chi_l, steps): the step count keeps the crystal inside the 500 x 500 box
STEPS = {0.7: 2600, 0.8: 2000, 0.9: 1500, 1.0: 1300}
SERIES = {
    "supersaturation": [(0.7, 0.1, 360), (0.8, 0.1, 360), (0.9, 0.1, 360), (1.0, 0.1, 360)],
    "anisotropy": [(0.8, 0.0, 360), (0.8, 0.025, 360), (0.8, 0.05, 360), (0.8, 0.1, 360)],
    "mobility_s0.8": [(0.8, 0.1, 360), (0.8, 0.1, 180), (0.8, 0.1, 108), (0.8, 0.1, 54), (0.8, 0.1, 18)],
    "mobility_s1.0": [(1.0, 0.1, 360), (1.0, 0.1, 240), (1.0, 0.1, 180), (1.0, 0.1, 108), (1.0, 0.1, 54)],
}


def name(s: float, s0: float, chi: float) -> str:
    return f"c3_s{s:.1f}_a{s0:g}_chi{chi:g}_fem"


def sub(text: str, pattern: str, repl: str) -> str:
    new, count = re.subn(pattern, repl, text, count=1, flags=re.M)
    if count != 1:
        raise SystemExit(f"pattern not found in the base deck: {pattern}")
    return new


def deck(s: float, s0: float, chi: float, member_of: list[str]) -> str:
    text = BASE.read_text()
    steps = STEPS[s]
    c_inf = C_L - s * (C_L - C_S)
    head = (
        "# Part IV, case C3 (FE route), parameter series of the orientation-field model:\n"
        f"# supersaturation S = {s:g}, anisotropy of the interface free energy s0 = {s0:g} (four-fold),\n"
        f"# orientational mobility chi_l = {chi:g} (360 = value of (granasy 2005 spherulites)).\n"
        f"# Member of the series: {', '.join(member_of)}. Conditions otherwise those of Fig. 4 of the paper.\n"
    )
    text = "\n".join(text.split("\n")[3:])
    text = head + text
    text = sub(text, r"c_inf = [0-9.]+;", f"c_inf = {c_inf:.6f};")
    text = sub(text, r"every 100 steps\.", f"every {steps // 20} steps.") if "every 100 steps." in text else text
    text = sub(text, r"^# EXPECTED RESULT\..*$",
               "# EXPECTED RESULT. See the series table of the Part IV report (orientation-field model, "
               "\"What to expect when the parameters change\").")
    text = sub(text, r"output_dir = .*", f"output_dir = results/crystallization/orientation/series/{name(s, s0, chi)}")
    text = sub(text, r"expression = 0\.399112\+[0-9.]+\*",
               f"expression = {C_S:.6f}+{c_inf - C_S:.6f}*")
    text = sub(text, r"anisotropy_strength  = .*", f"anisotropy_strength  = {s0:g}")
    text = sub(text, r"mobility_liquid    = .*",
               f"mobility_liquid    = {chi:g}       # chi_l (360 in the paper)")
    text = sub(text, r"n_steps         = \d+", f"n_steps         = {steps}")
    text = sub(text, r"snapshots_every = \d+", f"snapshots_every = {steps // 20}")
    text = text.replace("2000 steps, field series", f"{steps} steps, field series")
    return text


# Fig. 6 of the paper: fibrillar growth with fixed-angle branching. Two-fold kinetic coefficient
# with 99.5 % anisotropy, isotropic interface energy, branching energy with m = 3 (30 degrees for
# a two-fold crystal) and weight x = 0.15; the paper varies S. Added here: x and delta0 at S = 0.9.
FIBRIL_STEPS = {0.75: 2600, 0.85: 2000, 0.9: 1700, 0.95: 1500, 1.0: 1300, 1.1: 1100}
FIBRIL = {
    "supersaturation": [(s, 0.15, 0.995) for s in (0.75, 0.85, 0.9, 0.95, 1.0, 1.1)],
    "branching_weight": [(0.9, x, 0.995) for x in (0.0, 0.05, 0.15, 0.3)],
    "kinetic_anisotropy": [(0.9, 0.15, d) for d in (0.0, 0.5, 0.9, 0.995)],
}


def fibril_name(s: float, x: float, delta0: float) -> str:
    return f"c4_s{s:g}_x{x:g}_d{delta0:g}_fem"


def fibril_deck(s: float, x: float, delta0: float, member_of: list[str]) -> str:
    text = BASE.read_text()
    steps = FIBRIL_STEPS[s]
    c_inf = C_L - s * (C_L - C_S)
    head = (
        "# Part IV, case C4 (FE route): fibrillar growth with fixed-angle branching, from needle\n"
        "# crystals to category 2 spherulites, after Fig. 6 of (granasy 2005 spherulites).\n"
        f"# Supersaturation S = {s:g}; two-fold kinetic coefficient M = 1 + delta0 cos(2 vartheta - 2 pi theta)\n"
        f"# with delta0 = {delta0:g}; isotropic interface free energy; orientational energy\n"
        f"# x F0 + (1 - x) F1 with x = {x:g}, m = 3 (preferred misorientation 30 degrees), n = 1/2.\n"
        f"# Member of the series: {', '.join(member_of)}. Seed crystal at 30 degrees (theta = 1/6).\n"
        "# ASSUMED beyond the header below: chi_l = 360 (not stated for this figure).\n"
    )
    text = head + "\n".join(text.split("\n")[3:])
    text = text.replace("orientation in units of 2 pi / 4:", "orientation in units of 2 pi / 2 (two-fold crystal):")
    text = text.replace(
        "#     s = 1 + s0 cos(4 vartheta - 2 pi theta),  F(u) = |sin(pi u)| (u < 1/2),  u = xi0 |grad theta|",
        "#     s = 1 (isotropic); the phi equation is multiplied by M(vartheta, theta);\n"
        "#     F = x F0 + (1 - x) F1, F0 = |sin(2 pi m u)| (u < 3/(4m)), F1 = |sin(2 pi n u)| (u < 1/(4n)),\n"
        "#     both 1 beyond; u = xi0 |grad theta|",
    )
    text = sub(text, r"c_inf = [0-9.]+;", f"c_inf = {c_inf:.6f};")
    text = sub(text, r"^# EXPECTED RESULT\..*$",
               "# EXPECTED RESULT. See the series table of the Part IV report (fibrillar growth, "
               "\"What to expect when the parameters change\").")
    text = sub(text, r"output_dir = .*", f"output_dir = results/crystallization/orientation/fibril/{fibril_name(s, x, delta0)}")
    text = sub(text, r"expression = 0\.399112\+[0-9.]+\*", f"expression = {C_S:.6f}+{c_inf - C_S:.6f}*")
    text = sub(text, r"inside   = 0\.0", "inside   = 0.1666667     # 30 degrees in units of pi")
    text = sub(text, r"anisotropy_strength  = .*", "anisotropy_strength  = 0.0")
    text = sub(text, r"anisotropy_mode      = 4",
               "anisotropy_mode      = 2\n"
               f"    kinetic_anisotropy_strength = {delta0:g}     # delta0 (paper: 0.995)\n"
               "    kinetic_anisotropy_mode     = 2")
    text = sub(text, r"orientation          = theta .*", "orientation          = theta      # each grain has its frame")
    for _ in range(2):  # the coupling kernel and the flux kernel
        text = sub(text, r"branch_fraction    = 0\.0.*\n(\s*)branch_n           = 0\.5",
                   f"branch_fraction    = {x:g}        # x: weight of the branching energy F0\n"
                   "\\1branch_m           = 3.0\n\\1branch_n           = 0.5")
    text = sub(text, r"mobility_liquid    = .*", "mobility_liquid    = 360       # chi_l")
    text = sub(text, r"n_steps         = \d+", f"n_steps         = {steps}")
    text = sub(text, r"snapshots_every = \d+", f"snapshots_every = {steps // 20}")
    text = text.replace("2000 steps, field series", f"{steps} steps, field series")
    text = text.replace("every 100 steps.", f"every {steps // 20} steps.")
    return text


# Fig. 10 of the paper (its Table II): ten spherulitic morphologies. Columns: S, s0 (interface
# free energy), delta0 (kinetic coefficient), x, preferred misorientation in degrees (None: x = 0),
# m_phi (phase-field mobility), chi_l, grid of the paper. Both anisotropies are two-fold.
GALLERY = {
    "a": (1.0, 0.335, 0.75, 0.1, 30, 1.0, 144, 1000),
    "b": (0.75, 0.5, 0.5, 0.0, None, 0.9, 108, 2000),
    "c": (0.95, 0.0, 0.995, 0.275, 15, 1.0, 360, 1000),
    "d": (0.9, 0.0, 0.995, 0.15, 90, 1.0, 360, 2000),
    "e": (0.9, 0.75, 0.995, 0.15, 90, 1.0, 1440, 500),
    "f": (1.1, 0.0, 0.995, 0.2, 30, 1.0, 360, 1000),
    "g": (0.9, 0.5, 0.0, 0.0, None, 0.9, 360, 2000),
    "h": (1.1, 0.5, 0.0, 0.2, 15, 1.0, 360, 500),
    "i": (0.9, 0.5, 0.0, 0.0, None, 0.8, 180, 1000),
    "j": (0.875, 0.0, 0.995, 0.2, 15, 1.0, 18, 3000),
}
GALLERY_STEPS = {0.75: 2600, 0.875: 1850, 0.9: 1700, 0.95: 1500, 1.0: 1300, 1.1: 1100}


def gallery_deck(panel: str) -> str:
    s, s0, delta0, x, angle, m_phi, chi, n_paper = GALLERY[panel]
    FIBRIL_STEPS.setdefault(s, GALLERY_STEPS[s])
    text = fibril_deck(s, x, delta0, ["gallery"])
    m = 180.0 / (2.0 * angle) if angle else 3.0  # minimum of F0 at 1/(2m) of the 180 degree period
    # with s0 >= 0.5 the gradient term is up to (1 + s0)^2 times stiffer: half the step
    halve = s0 >= 0.5
    steps = FIBRIL_STEPS[s] * (2 if halve else 1)
    head = (
        f"# Part IV, case C6 (FE route): spherulitic morphology of panel ({panel}) of Fig. 10 of\n"
        "# (granasy 2005 spherulites), conditions of its Table II (reported_table):\n"
        f"# S = {s:g}; two-fold anisotropy of the interface free energy s0 = {s0:g} and of the kinetic\n"
        f"# coefficient delta0 = {delta0:g}; branching weight x = {x:g}"
        + (f" with a preferred misorientation of {angle} degrees (m = {m:g})" if angle else " (no branching minimum)")
        + f";\n# phase-field mobility m_phi = {m_phi:g}; orientational mobility chi_l = {chi:g}.\n"
        f"# The paper's grid is {n_paper} x {n_paper}; here 500 x 500 (the step count is not printed).\n"
        "# Seed crystal at 30 degrees (theta = 1/6): ASSUMED.\n"
        + ("# ASSUMED for stability of the explicit step with this strong anisotropy: half the time step\n"
           "# of the other cases (2.375e-6) and twice the steps.\n" if halve else "")
    )
    text = head + "\n".join(text.split("\n")[7:])
    text = text.replace("#     s = 1 (isotropic); the phi equation is multiplied by M(vartheta, theta);",
                        "#     s = 1 + s0 cos(2 vartheta - 2 pi theta); the phi equation is multiplied by\n"
                        "#     m_phi M(vartheta, theta), M = 1 + delta0 cos(2 vartheta - 2 pi theta);")
    text = sub(text, r"output_dir = .*", f"output_dir = results/crystallization/orientation/gallery/c6_fig10{panel}_fem")
    text = sub(text, r"anisotropy_strength  = 0\.0", f"anisotropy_strength  = {s0:g}")
    text = text.replace("    mobility             = 1.0\n    gradient_coefficient = 1.0",
                        f"    mobility             = {m_phi:g}\n    gradient_coefficient = 1.0")
    text, count = re.subn(r"(strength           = 1900\.0.*\n\s*)mobility           = 1\.0", rf"\g<1>mobility           = {m_phi:g}", text)
    assert count == 1, "coupling mobility"
    text, count = re.subn(r"(concentration = c\n\s*)mobility      = 1\.0", rf"\g<1>mobility      = {m_phi:g}", text)
    assert count == 1, "bulk mobility"
    text = text.replace("branch_m           = 3.0", f"branch_m           = {m:g}")
    text = sub(text, r"mobility_liquid    = .*", f"mobility_liquid    = {chi:g}       # chi_l (Table II)")
    text = sub(text, r"n_steps         = \d+", f"n_steps         = {steps}")
    text = sub(text, r"snapshots_every = \d+", f"snapshots_every = {steps // 20}")
    if halve:
        text = sub(text, r"dt              = 4\.75e-6", "dt              = 2.375e-6")
        text = sub(text, r"subcycles       = 'theta 20'", "subcycles       = 'theta 10'")
        # the noise amplitudes are increments per step: keep the increment per unit time
        text = sub(text, r"amplitude    = 1578\.7", "amplitude    = 2232.6")
        text = sub(text, r"amplitude             = 23679\.0", "amplitude             = 33487.2")
    return text


# Case E4: static and dynamic heterogeneities. A needle crystal (two-fold kinetic anisotropy,
# no branching minimum) against the same run with orientation pinning centres (foreign
# particles of fixed random orientation; 6.25 % of the nodes, the density of the 2000 x 2000
# run of Fig. 11 of the paper) and against a run with the orientational mobility reduced by 5.
PINNING = {"needle": (360.0, 0.0), "particles": (360.0, 0.0625), "low_mobility": (72.0, 0.0)}


def pinning_deck(kind: str, s: float = 0.85) -> str:
    chi, fraction = PINNING[kind]
    text = fibril_deck(s, 0.0, 0.995, ["pinning"])
    head = (
        "# Part IV, case E4 (FE route): duality of static and dynamic heterogeneities, after Fig. 11\n"
        "# (right block) of (granasy 2005 spherulites) and Fig. 6 of (granasy 2006 review).\n"
        f"# Variant '{kind}': supersaturation S = {s:g}; two-fold kinetic coefficient with\n"
        f"# delta0 = 0.995; no branching minimum (x = 0); chi_l = {chi:g}"
        + (f"; orientation pinning centres on {100 * fraction:g} % of the nodes\n"
           "# (random nodes that keep their initial, random orientation: foreign particles)." if fraction
           else ".")
        + "\n# reported_text: kinetic anisotropy 99.5 %, k = 2, isotropic interface free energy, 250,000\n"
        "# single-pixel pinning centres on 2000 x 2000 (6.25 %), mobility reduced by 5 in the third run.\n"
        "# assumed: S (not printed), the grid (500 x 500), the seed at 30 degrees, the step count.\n"
    )
    text = head + "\n".join(text.split("\n")[7:])
    text = sub(text, r"output_dir = .*", f"output_dir = results/crystallization/orientation/pinning/e4_{kind}_fem")
    text = sub(text, r"mobility_liquid    = .*", f"mobility_liquid    = {chi:g}       # chi_l")
    if fraction:
        text = sub(text, r"(  subcycles       = 'theta 20'.*\n)",
                   rf"\g<1>  pinned_fraction = 'theta {fraction:g}'   # orientation pinning centres\n  pinning_seed    = 7\n")
    return text


def main() -> None:
    pin_out = ROOT / "decks/crystallization/orientation/pinning"
    pin_out.mkdir(parents=True, exist_ok=True)
    FIBRIL_STEPS.setdefault(0.85, 2000)
    for kind in PINNING:
        (pin_out / f"e4_{kind}_fem.i").write_text(pinning_deck(kind))
    print(len(PINNING), "pinning decks")
    gal_out = ROOT / "decks/crystallization/orientation/gallery"
    gal_out.mkdir(parents=True, exist_ok=True)
    for panel in GALLERY:
        (gal_out / f"c6_fig10{panel}_fem.i").write_text(gallery_deck(panel))
    print(len(GALLERY), "gallery decks")
    fib_out = ROOT / "decks/crystallization/orientation/fibril"
    fib_out.mkdir(parents=True, exist_ok=True)
    fib: dict[tuple[float, float, float], list[str]] = {}
    for series, cases in FIBRIL.items():
        for case in cases:
            fib.setdefault(case, []).append(series)
    for (s, x, d), member_of in sorted(fib.items()):
        path = fib_out / f"{fibril_name(s, x, d)}.i"
        path.write_text(fibril_deck(s, x, d, member_of))
        print(path.relative_to(ROOT), member_of)
    print(len(fib), "fibril decks")
    OUT.mkdir(parents=True, exist_ok=True)
    members: dict[tuple[float, float, float], list[str]] = {}
    for series, cases in SERIES.items():
        for case in cases:
            members.setdefault(case, []).append(series)
    for (s, s0, chi), member_of in sorted(members.items()):
        path = OUT / f"{name(s, s0, chi)}.i"
        path.write_text(deck(s, s0, chi, member_of))
        print(path.relative_to(ROOT), member_of)
    print(len(members), "decks")


if __name__ == "__main__":
    main()
