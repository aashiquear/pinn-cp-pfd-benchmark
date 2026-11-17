#!/usr/bin/env python3
"""Write the standalone decks of case B1 (heterogeneous nucleus on a wall) for every undercooling.

`decks/crystallization/walls/b1_wall_model1_dT040_fem.i` is the hand-written base (Model I,
40 K). The other undercoolings change numbers only (driving force, classical radius, domain);
the Model II decks replace the wall condition by phi = 2/3 and come in pairs: the nucleus and
the wall layer without a nucleus, whose energy is the reference level of the barrier.
Every deck written is a complete input: `pinnacle run <deck>`.

Usage: python3 scripts/cryst_wall_decks.py
"""
from __future__ import annotations

import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "decks/crystallization/walls"
BASE = DIR / "b1_wall_model1_dT040_fem.i"
PSI = math.radians(61.2)
GAMMA = 1.0 / 12.0
# Ni of (granasy 2007 heterogeneous): L / w(T_m) with w(T_m) = 12 gamma / delta
T_M = 1728.0
L_OVER_W = 2.666e9 / (12.0 * 0.364 / (2.0e-9 / (2.0 * math.log(9.0))))
MODEL1 = range(10, 91, 10)  # 100 K: the constrained solve of the classical cap does not converge
MODEL2 = range(10, 81, 10)  # the layer loses stability at 92.0 K


def numbers(d_t: int) -> dict:
    big_d = -L_OVER_W * d_t / (T_M - d_t)
    radius = GAMMA / abs(big_d)
    lx = round(2.2 * radius / 5) * 5 + 10
    ly = round(1.3 * radius * (1 - math.cos(PSI)) / 5) * 5 + 15
    h = 0.25 if radius < 30 else 0.4
    return {"D": big_d, "R": radius, "yc": radius * math.cos(PSI), "Lx": lx, "Ly": ly, "h": h,
            "nx": int(round(lx / h)), "ny": int(round(ly / h)),
            "wcap": GAMMA**2 / abs(big_d) * (PSI - math.sin(PSI) * math.cos(PSI))}


def model1(d_t: int, base: str, b: dict) -> str:
    n = numbers(d_t)
    t = base.replace("b1_wall_model1_dT040_fem", f"b1_wall_model1_dT{d_t:03d}_fem")
    t = t.replace("undercooling 40 K", f"undercooling {d_t} K")
    t = t.replace(f"#       = {b['D']:.7f} at 40 K.", f"#       = {n['D']:.7f} at {d_t} K.")
    t = t.replace(f"on [0, {b['Lx']}] x [0, {b['Ly']}]", f"on [0, {n['Lx']}] x [0, {n['Ly']}]")
    t = t.replace(f"R = gamma / |Delta_f| = {b['R']:.4f} with its centre {b['yc']:.4f} below",
                  f"R = gamma / |Delta_f| = {n['R']:.4f} with its centre {n['yc']:.4f} below")
    assert f"= {n['R']:.4f} with its centre" in t
    t = t.replace(f"W_cap = {b['wcap']:.5f}.", f"W_cap = {n['wcap']:.5f}.")
    t = t.replace(f"domain    = '0 0 {b['Lx']} {b['Ly']}'", f"domain    = '0 0 {n['Lx']} {n['Ly']}'")
    t = t.replace(f"divisions = '{b['nx']} {b['ny']}'       # h = {b['h']} delta",
                  f"divisions = '{n['nx']} {n['ny']}'       # h = {n['h']} delta")
    t = t.replace(f"(y+{b['yc']:.4f})**2)-{b['R']:.4f})", f"(y+{n['yc']:.4f})**2)-{n['R']:.4f})")
    t = t.replace(f"df = {b['D']:.7f}", f"df = {n['D']:.7f}")
    assert f"{b['D']:.7f}" not in t or d_t == 40
    return t


M2_HEADER = """# Part IV, case B1 (FE route): heterogeneous crystal nucleus on a foreign wall, Model II,
# wall value phi_0 = 2/3 (contact angle 61.2 deg), undercooling {dT} K{what}
# (review Fig. 6, central row and dashed line).
#
# PHYSICS. Single-field model (granasy 2019 nucleation, Eq. 12, 17-19; granasy 2007
# heterogeneous), phi = 0 liquid, phi = 1 crystal, two dimensions:
#
#     F = int {{ eps^2/2 |grad phi|^2 + w g(phi) + p(phi) Delta_f }} dA
#
# Model II prescribes the order parameter on the wall, phi = phi_0, which orders the liquid
# next to it: a diffuse wall layer. The contact angle follows from
# cos(Psi) = -1 + 6 phi_0^2 - 4 phi_0^3 (phi_0 = 2/3: Psi = 61.2 deg). The barrier is the
# free energy of the nucleus on the wall minus that of the wall layer alone:
#
#     W = 2 [ F(nucleus deck) - F(layer deck) ]      (2: the mirror half)
#
# The layer exists while w g(phi) + p(phi) Delta_f > 0 for 0 < phi <= phi_0; at
# Delta_f = -w g(2/3)/p(2/3) = -w/64 it becomes unstable and the barrier vanishes: the
# surface spinodal, 92.0 K in the review.
#
# UNITS AND PARAMETERS. Scaled: lengths in delta = sqrt(2 eps^2 / w) = 0.4551 nm, energies
# per unit length in w(T) delta^2; eps^2 = 0.5, w = 1, gamma = 1/12. Ni data of the original
# paper (granasy 2007 heterogeneous, footnote 16; `reported_text`): gamma_SL = 0.364 J/m^2,
# interface thickness 2 nm, L = 2.666e9 J/m^3, T_m = 1728 K, eps^2 and w multiplied by T:
#     D = Delta_f / w(T) = -0.277781 Delta T / (1728 K - Delta T) = {D:.7f} at {dT} K.
#
# HOW IT IS SOLVED. Half of the symmetric state on [0, {Lx}] x [0, {Ly}]: phi = 2/3 on the
# wall y = 0, symmetry (natural) at x = 0, liquid on the top side, natural on the right
# side (the wall layer continues there). {solver}, started from
# {guess}.
#
# EXPECTED RESULT. `free_energy` of postprocessors.csv is F of this deck; the barrier needs
# the partner deck. Classical cap: W_cap = {wcap:.5f}. Model II lies well below the cap and
# reaches zero at the surface spinodal (review Fig. 6, dashed line); values in the report.
"""


def model2(d_t: int, base: str, layer: bool) -> str:
    n = numbers(d_t)
    name = f"b1_wall_model2_dT{d_t:03d}_{'layer' if layer else 'nucleus'}_fem"
    body = model1(d_t, base, numbers(40))
    body = body[body.index("[Global]"):]
    body = body.replace(f"b1_wall_model1_dT{d_t:03d}_fem", name)
    # wall: Dirichlet 2/3; right side natural
    i = body.index("[BCs]")
    j = body.index("[Executioner]")
    bcs = """[BCs]
  # Model II wall: the order parameter is prescribed
  [./wall]
    variable   = phi
    kind       = dirichlet
    node_set   = ymin
    components = [0]
    values     = [0.6666666666667]
  []
  [./liquid_top]
    variable   = phi
    kind       = dirichlet
    node_set   = ymax
    components = [0]
    values     = [0.0]
  []
[]

"""
    body = body[:i] + bcs + body[j:]
    # the wall layer at the melting point: phi = 1/(1 + exp(y + ln 2)) = 1/(1 + 0.5 exp(y))
    layer_expr = "1/(1+0.5*exp(y))"
    cap = f"0.5*(1-tanh((sqrt(x**2+(y+{n['yc']:.4f})**2)-{n['R']:.4f})/2))"
    # union of cap and layer without max(): 1 - (1 - cap)(1 - layer)
    expr = layer_expr if layer else f"1-(1-{cap})*(1-{layer_expr})"
    body = re.sub(r"expression = 0\.5\*\(1-tanh[^\n]*", f"expression = {expr}", body, count=1)
    body = body.replace("# circular cap of the classical radius (unquoted, no blanks)",
                        "# wall layer of the melting point" + ("" if layer else
                        ", with the circular cap of the classical radius on it"))
    if layer:
        body = body.replace("type            = pf_fem_saddle", "type            = pf_fem_transient")
    else:
        # The wall layer holds crystal-like order along the whole wall. The constrained
        # amount is therefore counted in a window around the nucleus only; with the total
        # amount the search finds the saddle of a uniformly thickened film instead.
        body = body.replace("[ICs]", f"""  # window of the constrained amount: 1 around the nucleus, 0 beyond 1.5 R along the wall
  [./window]
    kind       = parsed
    expression = 0.5*(1-tanh((x-{1.5 * n['R']:.4f})/2))
  []
[]

[ICs]""".replace("  # window", "PLACEHOLDER  # window"), 1)
        body = body.replace("[]\n\nPLACEHOLDER", "")
        body = body.replace("  type            = pf_fem_saddle",
                            "  type            = pf_fem_saddle\n  constraint_weight = window")
    # postprocessors: bulk free energy only
    i = body.index("  [./wall_energy]")
    j = body.index("[Output]")
    body = body[:i] + "[]\n\n" + body[j:]
    head = M2_HEADER.format(dT=d_t, D=n["D"], Lx=n["Lx"], Ly=n["Ly"], wcap=n["wcap"],
                            what=": the wall layer WITHOUT a nucleus (reference level)" if layer
                            else ": the nucleus on the wall",
                            solver="One backward-Euler step with dt = 1e12 (Newton on the stationary weak form; the layer is a minimum)"
                            if layer else "`pf_fem_saddle` (Newton at a fixed amount of crystal, the amount moved until "
                            "the constraint multiplier vanishes, then a plain Newton polish)",
                            guess="the wall layer of the melting point" if layer else
                            "that layer with the circular cap of the classical radius on it")
    return name, head + "\n" + body


def main() -> None:
    base = BASE.read_text()
    b = numbers(40)
    for d_t in MODEL1:
        if d_t != 40:
            (DIR / f"b1_wall_model1_dT{d_t:03d}_fem.i").write_text(model1(d_t, base, b))
    for d_t in MODEL2:
        for layer in (False, True):
            name, text = model2(d_t, base, layer)
            (DIR / f"{name}.i").write_text(text)
    print(len(list(DIR.glob("*.i"))), "decks in", DIR.relative_to(ROOT))


if __name__ == "__main__":
    main()
