#!/usr/bin/env python3
"""Write the FE decks of case F1 (hard-sphere nucleus) for a series of liquid volume fractions.

Every deck is complete (`pinnacle run <deck>`). The initial guess of the Newton solve is a tanh
droplet fitted to the shooting solution of scripts/cryst_hard_sphere_reference.py (run that
first); b and w are the calibrated values that script prints and stores.

Usage: python3 scripts/cryst_hard_sphere_decks.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "results/crystallization/hard_sphere"
OUT = ROOT / "decks/crystallization/hard_sphere"

TEMPLATE = """# Part IV, case F1 (FE route): critical nucleus of the fcc crystal in the hard-sphere liquid of
# volume fraction phi_inf = {phi_inf:.4f}, phase-field theory of (granasy 2003 hardsphere).
#
# PHYSICS. Structural order parameter m (1 liquid, 0 crystal) and volume fraction
# phi = pi sigma^3 rho / 6; units kT = 1, sigma = 1:
#
#     F = int {{ b/2 |grad m|^2 + w g(m) + [1 - p(m)] f_S(phi) + p(m) f_L(phi) }} dV
#     g = m^2 (1 - m)^2 / 4,   p = m^3 (10 - 15 m + 6 m^2)
#
# The particle number is conserved, so the nucleus is an extremum of the grand potential
# Omega = F - mu_inf N (Eq. 4-6 of the paper):
#
#     b (m'' + 2 m'/r) = w g'(m) + p'(m) [f_L(phi) - f_S(phi)]
#     [1 - p(m)] mu_S(phi) + p(m) mu_L(phi) = mu_inf       (no gradient: phi follows m locally)
#
# and its work of formation is W* = int 4 pi r^2 {{ b/2 m'^2 + Delta omega }} dr.
#
# PARAMETERS.
#   Equations of state: Hall (1972), rational form for the fluid, polynomial fit for the fcc
#   solid (`derived`: the paper cites Hall without printing the forms; this pair reproduces its
#   solid densities to four digits). Coexistence: liquid 0.494 (`reported_text`), solid at equal
#   pressure (0.5466; paper 0.546), equal chemical potential there.
#   b = {b:.5f} kT/sigma and w = {w:.5f} kT/sigma^3 (`derived`): from the interface free energy
#   {gamma:g} kT/sigma^2 and 10-90 % thickness {d1090:g} sigma of the "{interface}" interface
#   (`reported_text`) by Eq. 10-11 of the paper (scripts/cryst_hard_sphere_reference.py).
#
# HOW IT IS SOLVED. One backward-Euler step with a very large time step is Newton's method on
# the two stationary equations (kernels `hs_phase_bulk`, `hs_density_equilibrium`, `ac_interface`);
# the radial line carries the volume element 4 pi r^2. The start is a tanh droplet close to the
# nucleus: the nucleus is a saddle point and Newton needs a near start.
#
# EXPECTED RESULT (reference: shooting solution, scripts/cryst_hard_sphere_reference.py).
#   W* = {w_star:.3f} kT; m(0) = {m0:.4f}; phi(0) = {phi0:.4f}.{paper}
#   In results: `energy` of solver_history.csv is W* in kT.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/hard_sphere/{name}
[]

[Problem]
  type          = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type              = line_mesh
  domain            = '0.0 {r_max:g}'     # sigma
  divisions         = {divisions}          # h = 0.025 sigma
  coordinate_system = spherical
[]

[Variables]
  [./m]
    family     = fem
    components = 1
    order      = 2
  []
  [./phi]
    family     = fem
    components = 1
    order      = 2
  []
[]

[Functions]
  # tanh droplet fitted to the reference profile (unquoted, no blanks)
  [./droplet_m]
    kind       = parsed
    expression = 1-{amp:.5f}*0.5*(1-tanh((x-{radius:.5f})/{width:.5f}))
  []
  [./droplet_phi]
    kind       = parsed
    expression = {phi_inf:.4f}+{dphi:.5f}*0.5*(1-tanh((x-{radius:.5f})/{width:.5f}))
  []
[]

[ICs]
  [./m0]
    variable = m
    kind     = from_function
    function = droplet_m
  []
  [./phi0]
    variable = phi
    kind     = from_function
    function = droplet_phi
  []
[]

[Kernels]
  [./dm_dt]
    type        = time_derivative
    variable    = m
    coefficient = 1.0
  []
  [./m_gradient]
    type                 = ac_interface
    variable             = m
    mobility             = 1.0
    gradient_coefficient = {b:.5f}      # b
    spatial_dim          = 1
  []
  [./m_bulk]
    type                   = hs_phase_bulk
    variable               = m
    density                = phi
    well                   = {w:.5f}      # w
    liquid_volume_fraction = {phi_inf:.4f}
    fluid_eos              = hall
    solid_eos              = hall_fcc
    coexistence_liquid     = 0.494
  []
  [./dphi_dt]
    type        = time_derivative
    variable    = phi
    coefficient = 1.0
  []
  [./phi_balance]
    type                   = hs_density_equilibrium
    variable               = phi
    phase                  = m
    liquid_volume_fraction = {phi_inf:.4f}
    fluid_eos              = hall
    solid_eos              = hall_fcc
    coexistence_liquid     = 0.494
  []
[]

[BCs]
  [./far_liquid]
    variable   = m
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
  []
[]

[Executioner]
  type            = pf_fem_transient
  scheme          = bdf1
  n_steps         = 1
  dt              = 1.0e12
  max_newton_iter = 60
  newton_tol      = 1.0e-6     # the residual floor is 3e-8 (tabulated free energies, weights up to 4 pi 60^2)
  max_cutbacks    = 0
  snapshots_every = 1
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'm phi'
[]
"""


def droplet(r, amp, radius, width):
    return 1.0 - amp * 0.5 * (1.0 - np.tanh((r - radius) / width))


def main() -> None:
    profiles = np.load(REF / "reference_profiles.npz")
    OUT.mkdir(parents=True, exist_ok=True)
    with (REF / "reference.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        phi_inf, interface = float(row["phi_inf"]), row["interface"]
        b, w = float(row["b"]), float(row["w"])
        r, m, phi = profiles[f"{interface}_{phi_inf:.4f}"]
        (amp, radius, width), _ = curve_fit(droplet, r, m, p0=[1 - m[0], r[np.argmin(abs(m - 0.5 * (1 + m[0])))], 2.5])
        name = f"f1_nucleus_{interface}_phi{phi_inf:.4f}_fem"
        paper = ("\n#   Paper: W* = 40.63 kT (`reported_text`)."
                 if abs(phi_inf - 0.5207) < 1e-9 and interface == "average" else "")
        r_max = 60.0
        text = TEMPLATE.format(phi_inf=phi_inf, b=b, w=w, w_star=float(row["W_pft"]), m0=float(m[0]),
                               phi0=float(phi[0]), paper=paper, name=name, r_max=r_max,
                               divisions=int(r_max / 0.025), amp=amp, radius=radius, width=width,
                               dphi=float(phi[0] - phi_inf) * amp / (1 - m[0]),
                               gamma=float(row["gamma"]), d1090=float(row["d1090"]), interface=interface)
        (OUT / f"{name}.i").write_text(text)
    print(len(rows), "decks in", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
