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


def main() -> None:
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
