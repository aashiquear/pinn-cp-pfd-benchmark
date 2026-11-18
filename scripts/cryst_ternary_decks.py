#!/usr/bin/env python3
"""Write decks of case F3: ternary eutectic growth (two-dimensional sections, isothermal).

Model and numbers: (ratkai 2015 ternary), dimensionless as printed there. Every deck is complete
(`pinnacle run <deck>`).

Usage: python3 scripts/cryst_ternary_decks.py fs=<f_s,...> [lam=<lambda,...>] [c3=0.09] [width=96]
           [length=256] [steps=100000] [backend=numpy|torch]
  fs    free energy of the pure solids relative to the liquid (the same for the three components;
        the paper: -0.9640 at the cold end, -0.9140 at the dendrite tip): the driving force
  lam   spacing of the initial lamella pairs (cells)
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "decks/crystallization/ternary"

HEAD = """# Part IV, case F3 (FE route): ternary eutectic growth at a fixed driving force, two-dimensional
# section. Free energy of the pure solids f_s = {fs:g} (liquid: 0), liquid composition
# (c1, c2, c3) = ({c12:g}, {c12:g}, {c3:g}), initial lamellar spacing {lam:g} cells.
#
# PHYSICS. Ternary phase-field model of (ratkai 2015 ternary), dimensionless; phi = 1 liquid,
# phi = 0 solid (the paper uses the opposite labels), c3 = 1 - c1 - c2:
#
#     F = int {{ kappa_phi/2 |grad phi|^2 + w g(phi) + p(phi) f_l + [1 - p(phi)] [f_s + kappa_c/2 sum_i |grad c_i|^2] }} dV
#     f_l = sum_i c_i ln c_i,   f_s = sum_i c_i (f_s,i + ln c_i) + Omega_12 c1 c2 + Omega_23 c2 c3 + Omega_31 c3 c1
#     phi_t = -M_phi {{ w g' + p' [f_l - f_s - kappa_c/2 sum_i |grad c_i|^2] - kappa_phi lap phi }}
#     c_i,t = sum_j div[ p(phi) M_ij grad mu_j ],   M = 1 on the diagonal, -1/2 off it
#     mu_j  = ln c_j + [1 - p] (f_s,j + sum_k Omega_jk c_k) - kappa_c div([1 - p] grad c_j)
#
# Components 1 and 2 separate into two solid solutions (Omega_12 = 3.05 > 2); component 3 is
# soluble in neither (Omega_23 = Omega_31 = 3.0) and is rejected into the liquid, where it
# piles up ahead of the front and, in three dimensions, destabilizes it into two-phase dendrites.
#
# PARAMETERS (reported_text, the paper's reference conditions): kappa_phi = 0.75, kappa_c = 0.4,
#   w = 0.0469, M_phi = 1, Omega_12 = 3.05, Omega_23 = Omega_31 = 3.0, dx = 1, dt = 0.0025,
#   liquid composition (0.455, 0.455, 0.09). The paper pulls the sample at a fixed velocity
#   through a gradient of f_s (-0.9640 at the cold end, -0.9140 at the tip) in a 96 x 96 x 612
#   box; here the driving force is uniform and the section two-dimensional (`assumed` set-up).
#   assumed: the initial lamellae (solid slab of 40 cells, compositions 0.8225 / 0.0875 / 0.09).
#
# HOW IT IS SOLVED. `pf_fem_explicit`; the three chemical potentials are auxiliary variables
# (lumped projections), so the fourth-order composition equations need no third derivative.
# Periodic across the lamellae, no flux at the bottom and top.
#
# EXPECTED RESULT. See the Part IV report (ternary eutectic section): a lamellar front that
# rejects component 3; with rising driving force the front speeds up and traps it.
"""

BODY = """
[Global]
  backend    = {backend}
  use_gpu    = {gpu}
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/ternary/{name}
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 {width:g} {length:g}'
  divisions = '{nx} {ny}'
[]

[Variables]
{variables}[]

[Functions]
  # solid below y = 40; interface profile of width sqrt(2 kappa_phi / w) = 5.66
  [./front]
    kind       = parsed
    expression = 0.5*(1+tanh((y-40.0)/11.31))
  []
  [./comp1]
    kind       = parsed
    expression = {c12:g}+0.5*(1-tanh((y-40.0)/11.31))*{amp:g}*tanh(cos(2*pi*x/{lam:g})/0.3)
  []
  [./comp2]
    kind       = parsed
    expression = {c12:g}-0.5*(1-tanh((y-40.0)/11.31))*{amp:g}*tanh(cos(2*pi*x/{lam:g})/0.3)
  []
[]

[ICs]
  [./phi0]
    variable = phi
    kind     = from_function
    function = front
  []
  [./c10]
    variable = c1
    kind     = from_function
    function = comp1
  []
  [./c20]
    variable = c2
    kind     = from_function
    function = comp2
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
    mobility             = 1.0
    gradient_coefficient = 0.75
    anisotropy_strength  = 0.0
  []
  [./phi_bulk]
    type                             = ts_phase_bulk
    variable                         = phi
    well                             = 0.0469
    mobility                         = 1.0
    composition_gradient_coefficient = 0.4
{thermo}  []
  [./dc1_dt]
    type        = time_derivative
    variable    = c1
    coefficient = 1.0
  []
  [./c1_flux]
    type                = ts_solute_flux
    variable            = c1
    phase               = phi
    chemical_potentials = 'mu1 mu2 mu3'
    mobilities          = '1.0 -0.5 -0.5'
  []
  [./dc2_dt]
    type        = time_derivative
    variable    = c2
    coefficient = 1.0
  []
  [./c2_flux]
    type                = ts_solute_flux
    variable            = c2
    phase               = phi
    chemical_potentials = 'mu1 mu2 mu3'
    mobilities          = '-0.5 1.0 -0.5'
  []
{potentials}[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = {steps}
  dt                  = 0.0025
  periodic_axes       = 'x'
  auxiliary_variables = 'mu1 mu2 mu3'
  snapshots_every     = {snap}
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c1 c2'
[]
"""


def main() -> None:
    opts = dict(a.split("=") for a in sys.argv[1:])
    width, length = float(opts.get("width", 96)), float(opts.get("length", 256))
    steps, c3 = int(opts.get("steps", 100000)), float(opts.get("c3", 0.09))
    backend = opts.get("backend", "numpy")
    c12 = (1.0 - c3) / 2.0
    OUT.mkdir(parents=True, exist_ok=True)
    for fs in [float(v) for v in opts["fs"].split(",")]:
        for lam in [float(v) for v in opts.get("lam", "16").split(",")]:
            thermo = (f"    compositions                     = 'c1 c2'\n"
                      f"    solid_energies                   = '{fs:g} {fs:g} {fs:g}'\n"
                      f"    interactions                     = '3.05 3.0 3.0'\n")
            variables = "".join(f"  [./{v}]\n    family     = fem\n    components = 1\n    order      = 2\n  []\n"
                                for v in ("phi", "c1", "c2", "mu1", "mu2", "mu3"))
            potentials = "".join(
                f"  [./potential{j}]\n    type                 = ts_chemical_potential\n"
                f"    variable             = mu{j}\n    phase                = phi\n"
                f"    component            = {j}\n    gradient_coefficient = 0.4\n"
                f"    compositions         = 'c1 c2'\n    solid_energies       = '{fs:g} {fs:g} {fs:g}'\n"
                f"    interactions         = '3.05 3.0 3.0'\n  []\n" for j in (1, 2, 3))
            name = f"f3_ternary_fs{fs:g}_c3{c3:g}_lam{lam:g}_fem"
            text = HEAD.format(fs=fs, c12=c12, c3=c3, lam=lam) + BODY.format(
                backend=backend, gpu="true" if backend == "torch" else "false", name=name, width=width,
                length=length, nx=int(width), ny=int(length), variables=variables, c12=c12,
                amp=0.3675 * (1.0 - c3) / 0.91, lam=lam, thermo=thermo, potentials=potentials,
                steps=steps, snap=max(1, steps // 20))
            (OUT / f"{name}.i").write_text(text)
            print(name)


if __name__ == "__main__":
    main()
