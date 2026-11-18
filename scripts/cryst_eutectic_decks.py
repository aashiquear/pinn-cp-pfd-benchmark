#!/usr/bin/env python3
"""Write decks of case F2: lamellar eutectic growth of Ag-Cu at a fixed temperature.

One lamella pair of spacing lambda in a box periodic across the lamellae; the liquid of
composition c0 lies above a solid slab of the two solid solutions. Every deck is complete
(`pinnacle run <deck>`); the numbers come from scripts/cryst_eutectic_thermo.py.

Usage: python3 scripts/cryst_eutectic_decks.py <T in K> <c0> <lambda in nm> [...lambda]
       optional trailing "steps=<n>" and "length=<nm>"
       python3 scripts/cryst_eutectic_decks.py grain [T=900] [c0=0.345] [size=48] [steps=150000]
              [S=0.1,...] [chi=0.03,...] [backend=torch]      (case F2b, with the orientation field)
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cryst_eutectic_thermo import eutectic, groups  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "decks/crystallization/eutectic"
DX, DT = 0.25, 0.004

TEMPLATE = """# Part IV, case F2 (FE route): lamellar eutectic growth of Ag-Cu at T = {T:g} K
# (undercooling {dT:.1f} K below the eutectic temperature {TE:.1f} K), liquid composition
# c0 = {c0:g} (mole fraction of Cu), lamellar spacing lambda = {lam:g} nm.
#
# PHYSICS. Binary phase-field model with regular-solution thermodynamics, (lewis 2004 eutectic)
# and (ratkai 2017 nanoscale); phi = 1 liquid, phi = 0 solid, c = Cu fraction; energy in R T/v:
#
#     F = int {{ kappa_phi/2 |grad phi|^2 + kappa_c/2 |grad c|^2 + w(c) g(phi)
#               + [1 - p(phi)] f_S(c) + p(phi) f_L(c) }} dV
#     f_L = omega_L c (1 - c) + c ln c + (1 - c) ln(1 - c)
#     f_S = (1 - c) f_A + c f_B + omega_S c (1 - c) + c ln c + (1 - c) ln(1 - c)
#     phi_t = -m_phi {{ w(c) g' + p' [f_L - f_S] - kappa_phi lap phi }}
#     c_t   = div{{ D(phi) c (1 - c) grad mu }},   D = r + (1 - r) p(phi)
#     mu    = (w_B - w_A) g + [1 - p] f_S' + p f_L' - kappa_c lap c
#
# omega_S > 2: the solid has a miscibility gap (c_alpha = {c_alpha:.4f}, c_beta = {c_beta:.4f} at this
# temperature); the liquid solidifies into alternating lamellae of the two solid solutions,
# exchanging Cu sideways through the liquid ahead of the front. The gradient energy of c gives
# the lamella boundaries their energy; the composition equation is fourth order.
#
# UNITS. Length xi = 1 nm, time xi^2/D_L = 7.143e-10 s, energy density R T/v = {scale:.4g} J/m^3.
#
# PARAMETERS (provenance: docs/14_phase_field_crystallization_plan.md, sections 9.3, 9.4;
# computed by scripts/cryst_eutectic_thermo.py).
#   reported_text: L_Ag = 1097 J/cm^3, L_Cu = 1728 J/cm^3, v = 8.7 cm^3/mol, T_Ag = 1234 K,
#     T_Cu = 1358 K, Omega_L v = 16,000 - 4 T, Omega_S v = 26,000 - 4 T J/mol; gamma = 172 (Ag),
#     227 (Cu) mJ/m^2; 10-90 % interface thickness 1 nm; eps_c^2 = eps_phi^2/4; D_L = 1.4e-9 m^2/s;
#     M_phi = 0.25 m^3/(J s); D_S = D_L/100.
#   derived: f_A = {solid_a:.6f}, f_B = {solid_b:.6f}, omega_S = {omega_solid:.5f}, omega_L = {omega_liquid:.5f},
#     kappa_phi = {kappa_phi:.6f}, kappa_c = {kappa_c:.6f}, w_A = {well_a:.5f}, w_B = {well_b:.5f}, m_phi = {m_phi:.6f}.
#   assumed: the initial lamella pair (solid slab of 10 nm with the lever-rule fraction
#     {frac:.3f} of the Cu-rich phase), mesh size 0.25 nm, time step 0.004.
#
# HOW IT IS SOLVED. `pf_fem_explicit` (forward Euler, lumped mass, bilinear elements); mu is an
# auxiliary variable: the lumped projection of its kernel terms each step, so no fourth
# derivative is formed. Periodic across the lamellae, no flux at the bottom and top.
#
# EXPECTED RESULT. The two-phase front advances at a steady speed that depends on lambda (a
# maximum at an intermediate spacing; too narrow a pair is overgrown by one phase).
# The front speed is the slope of the solid fraction against time times the box length:
# `phi_mean` of solver_history.csv. Report: Part IV, eutectic section.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/eutectic/{name}
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 {lam:g} {length:g}'      # nm: one lamella pair across, growth along y
  divisions = '{nx} {ny}'                   # 0.25 nm
[]

[Variables]
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
  [./c]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # solid below y = 10 nm; the Cu-rich lamella occupies the central fraction of the width
  [./front]
    kind       = parsed
    expression = 0.5*(1+tanh((y-10.0)/0.4551))
  []
  [./composition]
    kind       = parsed
    expression = {c0:g}*0.5*(1+tanh((y-10.0)/0.4551))+0.5*(1-tanh((y-10.0)/0.4551))*({c_alpha:.6f}+{dc:.6f}*0.5*(1+tanh((cos(2*pi*(x/{lam:g}-0.5))-({cosf:.6f}))/0.15)))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = front
  []
  [./c0]
    variable = c
    kind     = from_function
    function = composition
  []
  [./mu0]
    variable = mu
    kind     = constant
    value    = 0.0
  []
[]

[Kernels]
  [./dphi_dt]
    type        = time_derivative
    variable    = phi
    coefficient = 1.0
  []
  [./phi_gradient]
    type                 = ac_interface_anisotropic
    variable             = phi
    mobility             = {m_phi:.6f}
    gradient_coefficient = {kappa_phi:.6f}
    anisotropy_strength  = 0.0
  []
  [./phi_bulk]
    type          = rs_phase_bulk
    variable      = phi
    concentration = c
    mobility      = {m_phi:.6f}
    solid_a       = {solid_a:.6f}
    solid_b       = {solid_b:.6f}
    omega_solid   = {omega_solid:.5f}
    omega_liquid  = {omega_liquid:.5f}
    well_a        = {well_a:.5f}
    well_b        = {well_b:.5f}
  []
  [./dc_dt]
    type        = time_derivative
    variable    = c
    coefficient = 1.0
  []
  [./solute]
    type               = rs_solute_flux
    variable           = c
    phase              = phi
    chemical_potential = mu
    diffusivity        = 1.0
    diffusivity_ratio  = 0.01
  []
  [./potential]
    type                 = rs_chemical_potential
    variable             = mu
    phase                = phi
    concentration        = c
    gradient_coefficient = {kappa_c:.6f}
    solid_a              = {solid_a:.6f}
    solid_b              = {solid_b:.6f}
    omega_solid          = {omega_solid:.5f}
    omega_liquid         = {omega_liquid:.5f}
    well_a               = {well_a:.5f}
    well_b               = {well_b:.5f}
  []
[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = {steps}
  dt                  = {dt:g}
  periodic_axes       = 'x'
  auxiliary_variables = 'mu'
  snapshots_every     = {snap}
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c mu'
[]
"""


def main() -> None:
    args = [a for a in sys.argv[1:] if "=" not in a]
    opts = dict(a.split("=") for a in sys.argv[1:] if "=" in a)
    t, c0 = float(args[0]), float(args[1])
    steps, length = int(opts.get("steps", 100000)), float(opts.get("length", 60.0))
    te = eutectic()[3]
    gr = groups(t)
    frac = (c0 - gr["c_alpha"]) / (gr["c_beta"] - gr["c_alpha"])
    OUT.mkdir(parents=True, exist_ok=True)
    for lam in [float(a) for a in args[2:]]:
        nx = int(round(lam / DX))
        name = f"f2_lamellar_T{t:g}_c{c0:g}_lam{lam:g}_fem"
        text = TEMPLATE.format(dT=te - t, TE=te, c0=c0, lam=lam, scale=gr["energy_unit_J_m3"], frac=frac,
                               name=name, length=length, nx=nx, ny=int(round(length / DX)),
                               dc=gr["c_beta"] - gr["c_alpha"], cosf=float(np.cos(np.pi * frac)),
                               steps=steps, dt=DT, snap=max(1, steps // 20), **gr)
        (OUT / f"{name}.i").write_text(text)
        print(name, f"nx {nx} fraction {frac:.3f}")


# ---------------------------------------------------------------------------------------------
# Case F2b: a eutectic grain with the orientation field (two solids at a locked misorientation)
# ---------------------------------------------------------------------------------------------
def grain_deck(t: float, c0: float, size: float, steps: int, strength: float, chi: float,
               backend: str = "torch") -> tuple[str, str]:
    """Deck text and name: a seed of the Ag-rich solid in the liquid, with theta."""
    te = eutectic()[3]
    gr = groups(t)
    n = int(round(size / DX))
    name = f"f2b_grain_T{t:g}_c{c0:g}_S{strength:g}_chi{chi:g}_fem"
    text = TEMPLATE.format(dT=te - t, TE=te, c0=c0, lam=size, scale=gr["energy_unit_J_m3"], frac=0.0,
                           name=name, length=size, nx=n, ny=n, dc=0.0, cosf=0.0, steps=steps, dt=DT,
                           snap=max(1, steps // 20), **gr)
    head = f"""# Part IV, case F2b (FE route): growth of a eutectic grain of Ag-Cu at T = {t:g} K
# (undercooling {te - t:.1f} K) from a seed of the Ag-rich solid, liquid composition c0 = {c0:g},
# with an orientation field: the two solid solutions keep a fixed orientation relation,
# after (lewis 2004 eutectic).
#
# PHYSICS. The binary regular-solution model of case F2 (equations below) plus the orientation
# theta in [0, 1) with the composition-dependent energy
#
#     f_ori = [1 - p(phi)] S E(c, u),   u = d |grad theta|
#     E = h(c) u + [1 - h(c)] {{A + B |cos(2 pi m u)|}},   h = (1/2){{1 + cos[2 pi (c - c_alpha)/(c_beta - c_alpha)]}}
#
# Inside either solid (h = 1) E is the grain-boundary energy, linear in the misorientation;
# across a boundary between the two solids (h -> 0) it is lowest at u = 1/(4 m) = 1/2:
# a step of half a period in theta. theta_t = chi(phi) div{{dE-flux}}, chi = chi_l p(phi).
# reported_text (lewis 2004 eutectic): A = 0.9, B = 0.1, m = 1/2, psi = 0, d = dx.
# assumed (not printed anywhere): S = {strength:g}, chi_l = {chi:g}, no noise, the seed.
#
"""
    text = head + "\n".join(text.split("\n")[3:])
    text = text.replace(
        "# EXPECTED RESULT. The two-phase front advances at a steady speed that depends on lambda (a\n"
        "# maximum at an intermediate spacing; too narrow a pair is overgrown by one phase).\n"
        "# The front speed is the slope of the solid fraction against time times the box length:\n"
        "# `phi_mean` of solver_history.csv. Report: Part IV, eutectic section.",
        "# EXPECTED RESULT. The seed grows, rejects Cu, and the Cu-rich solid forms on it; the two\n"
        "# solids then grow side by side. Where they touch, theta differs by 0.5.\n"
        "# Report: Part IV, eutectic section.")
    text = text.replace(f"domain    = '0 0 {size:g} {size:g}'      # nm: one lamella pair across, growth along y",
                        f"domain    = '-{size / 2:g} -{size / 2:g} {size / 2:g} {size / 2:g}'      # nm")
    text = text.replace("backend    = numpy\n  use_gpu    = false",
                        f"backend    = {backend}\n  use_gpu    = {'true' if backend == 'torch' else 'false'}")
    a, b = text.index("[Functions]"), text.index("[Kernels]")
    seed = f"""[Functions]
  # seed of the Ag-rich solid: radius 4 nm
  [./seed_phi]
    kind       = parsed
    expression = 0.5*(1+tanh((sqrt(x**2+y**2)-4.0)/0.4551))
  []
  [./seed_c]
    kind       = parsed
    expression = {gr['c_alpha']:.6f}+{c0 - gr['c_alpha']:.6f}*0.5*(1+tanh((sqrt(x**2+y**2)-4.0)/0.4551))
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = seed_phi
  []
  [./c0]
    variable = c
    kind     = from_function
    function = seed_c
  []
  [./mu0]
    variable = mu
    kind     = constant
    value    = 0.0
  []
  # no orientational order in the liquid; the seed has theta = 0.25
  [./theta0]
    variable = theta
    kind     = uniform_random
    minimum  = 0.0
    maximum  = 1.0
    seed     = 20261118
    center   = '0 0'
    radius   = 4.0
    inside   = 0.25
  []
[]

"""
    text = text[:a] + seed + text[b:]
    text = text.replace("  [./mu]\n    family     = fem\n    components = 1\n    order      = 2\n  []\n[]",
                        "  [./mu]\n    family     = fem\n    components = 1\n    order      = 2\n  []\n"
                        "  [./theta]\n    family     = fem\n    components = 1\n    order      = 2\n  []\n[]")
    energy = (f"    strength            = {strength:g}        # S (assumed)\n"
              f"    boundary_thickness  = {DX:g}       # d = dx\n"
              f"    solid_composition_a = {gr['c_alpha']:.6f}\n"
              f"    solid_composition_b = {gr['c_beta']:.6f}\n"
              "    offset              = 0.9\n    amplitude           = 0.1\n"
              "    misorientation_index = 0.5\n    regularization      = 1.0e-3\n")
    extra = f"""  [./phi_orientation]
    type          = eo_phase_coupling
    variable      = phi
    orientation   = theta
    concentration = c
    mobility      = {gr['m_phi']:.6f}
{energy}  []
  [./mu_orientation]
    type          = eo_potential_coupling
    variable      = mu
    phase         = phi
    orientation   = theta
    concentration = c
{energy}  []
  [./dtheta_dt]
    type        = time_derivative
    variable    = theta
    coefficient = 1.0
  []
  [./theta_flux]
    type             = eo_grain_flux
    variable         = theta
    phase            = phi
    concentration    = c
    mobility_liquid  = {chi:g}       # chi_l (assumed)
    mobility_crystal = 0.0
{energy}  []
[]

[Executioner]"""
    text = text.replace("[]\n\n[Executioner]", extra, 1)
    text = text.replace("  periodic_axes       = 'x'\n", "  periodic_axes       = 'x y'\n  angular_variables   = 'theta'\n")
    text = text.replace("fields  = 'phi c mu'", "fields  = 'phi c theta'")
    return text, name


def grain_main(argv: list[str]) -> None:
    opts = dict(a.split("=") for a in argv)
    t, c0 = float(opts.get("T", 900.0)), float(opts.get("c0", 0.345))
    out = OUT / "grain"
    out.mkdir(parents=True, exist_ok=True)
    for strength in [float(v) for v in opts.get("S", "0.1").split(",")]:
        for chi in [float(v) for v in opts.get("chi", "0.03").split(",")]:
            text, name = grain_deck(t, c0, float(opts.get("size", 48.0)), int(opts.get("steps", 150000)),
                                    strength, chi, opts.get("backend", "torch"))
            seed_lam, seed_r, sub = float(opts.get("seed_lam", 0.0)), float(opts.get("seed_r", 4.0)), int(opts.get("sub", 1))
            if seed_lam > 0.0:  # two-phase seed: a disc cut out of a lamellar structure
                gr = groups(t)
                frac = (c0 - gr["c_alpha"]) / (gr["c_beta"] - gr["c_alpha"])
                old_name, name = name, name.replace("_fem", f"_lam{seed_lam:g}_fem")
                text = text.replace(old_name, name)
                hull = f"0.5*(1+tanh((sqrt(x**2+y**2)-{seed_r:g})/0.4551))"
                stripes = (f"({gr['c_alpha']:.6f}+{gr['c_beta'] - gr['c_alpha']:.6f}*0.5*(1+tanh((cos(2*pi*x/{seed_lam:g})"
                           f"-{math.cos(math.pi * frac):.6f})/0.3)))")
                lines = text.split("\n")
                for i, line in enumerate(lines):
                    if line.startswith("    expression = ") and "tanh((sqrt(x**2+y**2)" in line:
                        lines[i] = (f"    expression = {hull}" if "seed_phi" in lines[i - 2] else
                                    f"    expression = {stripes}*(1-{hull})+{c0:g}*{hull}")
                text = "\n".join(lines)
                text = text.replace("  # seed of the Ag-rich solid: radius 4 nm",
                                    f"  # two-phase seed: a disc of radius {seed_r:g} nm cut out of lamellae of spacing "
                                    f"{seed_lam:g} nm along x\n  # (Cu-rich fraction {frac:.3f}, the lever rule)")
                text = text.replace("    radius   = 4.0\n", f"    radius   = {seed_r:g}\n")
                text = text.replace("from a seed of the Ag-rich solid,", "from a two-phase (lamellar) seed,")
                text = text.replace("# EXPECTED RESULT. The seed grows, rejects Cu, and the Cu-rich solid forms on it; the two\n"
                                    "# solids then grow side by side. Where they touch, theta differs by 0.5.",
                                    "# EXPECTED RESULT. The lamellae of the seed grow outward side by side; new Cu-rich solid\n"
                                    "# takes an orientation half a period from the Ag-rich solid next to it.")
            if sub > 1:
                text = text.replace("  angular_variables   = 'theta'\n",
                                    f"  angular_variables   = 'theta'\n  subcycles           = 'theta {sub}'\n")
            text = text.replace(f"results/crystallization/eutectic/{name}", f"results/crystallization/eutectic/grain/{name}")
            (out / f"{name}.i").write_text(text)
            print("grain/" + name)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "grain":
        grain_main(sys.argv[2:])
    else:
        main()
