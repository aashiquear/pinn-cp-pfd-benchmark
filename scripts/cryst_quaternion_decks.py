#!/usr/bin/env python3
"""Write decks of case F4: polycrystalline solidification in three dimensions (quaternions).

Binary Cu-Ni model of the two-dimensional cases (same dimensionless numbers) with the
orientation as a unit quaternion field, after (pusztai 2005 threedimensions). Every deck is
complete (`pinnacle run <deck>`).

Usage: python3 scripts/cryst_quaternion_decks.py <case> n=<cells per edge> [steps=..] [S=..]
           [chi=..] [gamma=..] [seeds=..] [backend=torch|numpy]
  case = dendrites   several seeds of random orientation, cubic kinetic anisotropy (paper Fig. 2)
         spherulite  one seed, reduced orientational mobility, large driving force (paper Fig. 4)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "decks/crystallization/quaternion"
DX, DT = 6.25e-3, 4.75e-6
C_L, C_S = 0.466219, 0.399112

TEMPLATE = """# Part IV, case F4 (FE route): {title}
# after (pusztai 2005 threedimensions). {n} x {n} x {n} cells, supersaturation S = {S:g}.
#
# PHYSICS. Binary ideal-solution model of the two-dimensional cases (phi = 1 liquid, phi = 0
# crystal, c = Cu fraction) in three dimensions; the orientation is a unit quaternion
# q = (q0, q1, q2, q3), q0 = cos(delta/2), q_i = c_i sin(delta/2) (rotation by delta about c):
#
#     f_ori = [1 - p(phi)] h F(u),   u = xi0 |grad q| k / pi,   F = |sin(pi u)| (u < 1/2), 1 beyond
#     phi_t = M(n') {{ lap phi - a_w(c) g'(phi) - p'(phi) [ b(c) - h F(u) ] }} + zeta_phi
#     c_t   = div{{ p(phi) [ grad c - c (1 - c) (q_f p' + q_w g') grad phi ] }}
#     q_i,t = chi_l p(phi) {{ div J_i - q_i sum_k q_k div J_k }},  J_i = [1 - p] (pi/2) cos(pi u) xi0 grad q_i / u
#     M(n') = 1 - 3 gamma + 4 gamma (n'_1^4 + n'_2^4 + n'_3^4),  n' = R(q)^T grad phi / |grad phi|
#
# 2 xi0 |grad q| is the misorientation angle across one correlation length; for a rotation
# about a fixed axis the model is the two-dimensional one with theta in units of 2 pi / k.
# The orientation is random in the liquid, orders in the interface and is frozen in the
# crystal. The cubic mobility is taken in the frame of the local crystal, so every grain
# grows along its own <100> directions. Crystal symmetries are not applied to the
# orientational energy (triclinic, as in the paper's runs).
#
# UNITS. Lengths in xi = 2.1e-4 cm, time in xi^2/D_l (as the two-dimensional Cu-Ni cases).
#
# PARAMETERS. Binary model: as case C2 (reported_table / derived, plan section 7.2).
#   reported_text (the paper): S = {S_paper}; M_q reduced by 3 to 3.5 for the spherulite.
#   assumed: h = 1900 and chi_l as in the two-dimensional cases (the paper prints mobilities in
#   other units and no H), k = 4, gamma = {gamma:g} (the paper prints no magnitude; 0.2 is quoted
#   for Cu-Ni by the same authors in a later proceedings paper), the seeds, the noise
#   (standard deviation 0.0025 of the nodal increment per step, liquid only), the box size
#   (paper: 300^3; here {n}^3).
#
# HOW IT IS SOLVED. `pf_fem_explicit` on trilinear elements, dx = 6.25e-3, dt = 4.75e-6,
# periodic. q is a unit-norm variable (tangent projection of the rate, normalisation, q and
# -q identified) advanced in 20 sub-steps per step. {steps} steps, field series every {snap}.
#
# EXPECTED RESULT. {expected}

[Global]
  backend    = {backend}
  use_gpu    = {gpu}
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/quaternion/{name}
[]

[Problem]
  type          = composed
  dimension     = 3
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 3
  domain    = '-{half:g} -{half:g} -{half:g} {half:g} {half:g} {half:g}'
  divisions = '{n} {n} {n}'
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
  [./q]
    family     = fem
    components = 4
    order      = 2
  []
[]

[Functions]
  # seeds of radius 8 dx: the product of the tanh profiles of the single seeds (unquoted, no blanks)
  [./seed_phi]
    kind       = parsed
    expression = {phi_expr}
  []
  [./seed_c]
    kind       = parsed
    expression = {c_s}+{dc:.6f}*({phi_expr})
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
  # random orientations in the liquid; every seed has one orientation (random, from seed + 1)
  [./q0]
    variable = q
    kind     = random_unit_quaternion
    seed     = 20261118
    centers  = '{centers}'
    radius   = 0.05
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
    gradient_coefficient = 1.0
    anisotropy_strength  = 0.0
  []
  [./phi_mobility]
    type        = cubic_kinetic_anisotropy
    variable    = phi
    strength    = {gamma:g}
    orientation = q
  []
  [./phi_orientation]
    type               = quaternion_phase_coupling
    variable           = phi
    orientation        = q
    correlation_length = 6.25e-3
    period_index       = 4
    strength           = 1900.0
    mobility           = 1.0
  []
  [./phi_bulk]
    type          = wb_phase_bulk
    variable      = phi
    concentration = c
    mobility      = 1.0
    well_a        = 10204.6
    well_b        = 10181.9
    driving_a     = 109551.0
    driving_b     = -143771.0
  []
  [./phi_noise]
    type         = langevin_noise
    variable     = phi
    amplitude    = {noise:.1f}       # 0.0025 / (dt x 0.1764): measured factor of the lumped trilinear load
    distribution = normal
    state_factor = phi**3*(10-15*phi+6*phi**2)
    seed         = 20261118
  []
  [./dc_dt]
    type        = time_derivative
    variable    = c
    coefficient = 1.0
  []
  [./solute]
    type              = wb_solute_flux
    variable          = c
    phase             = phi
    diffusivity       = 1.0
    diffusivity_ratio = 0.0
    partition_driving = 0.273838
    partition_well    = 2.46e-5
  []
  [./dq_dt]
    type        = time_derivative
    variable    = q
    coefficient = 1.0
  []
  [./q_flux]
    type               = quaternion_grain_flux
    variable           = q
    phase              = phi
    correlation_length = 6.25e-3
    period_index       = 4
    mobility_liquid    = {chi:g}
    mobility_crystal   = 0.0
    regularization     = 1.0e-3
  []
[]

[Executioner]
  type                = pf_fem_explicit
  n_steps             = {steps}
  dt                  = 4.75e-6
  periodic_axes       = 'x y z'
  unit_norm_variables = 'q'
  subcycles           = 'q 20'
  snapshots_every     = {snap}
  energy_every        = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c q'
[]
"""


def main() -> None:
    case = sys.argv[1]
    opts = dict(a.split("=") for a in sys.argv[2:])
    n = int(opts["n"])
    backend = opts.get("backend", "torch")
    half = n * DX / 2
    if case == "dendrites":
        s, chi, gamma = float(opts.get("S", 0.78)), float(opts.get("chi", 360)), float(opts.get("gamma", 0.2))
        k = int(opts.get("seeds", 4))
        rng = np.random.default_rng(7)
        # well-separated seeds: corners of a tetrahedron-like set, jittered
        base = np.array([[-1, -1, -1], [1, 1, -1], [1, -1, 1], [-1, 1, 1], [0, 0, 0], [1, 1, 1], [-1, -1, 1], [1, -1, -1]])
        centers = (base[:k] * 0.25 * 2 * half + rng.uniform(-0.03, 0.03, (k, 3)) * 2 * half)
        title = f"growth of {k} dendritic crystals of random orientation in three dimensions,"
        expected = ("Each seed grows into a dendrite with six primary arms along the <100> directions of its own "
                    "lattice; the arms of different crystals point in unrelated directions (compare Fig. 2).")
        s_paper = "0.78 for dendrites"
        steps = int(opts.get("steps", 1500))
    elif case == "spherulite":
        s, chi, gamma = float(opts.get("S", 1.0)), float(opts.get("chi", 108)), float(opts.get("gamma", 0.2))
        centers = np.zeros((1, 3))
        k = 1
        title = "a polycrystalline spherulite in three dimensions by growth front nucleation,"
        expected = ("With the orientational mobility reduced, the front cannot order the orientation: new grains "
                    "form at the front and the crystal is a compact, nearly spherical polycrystal (compare Fig. 4).")
        s_paper = "1.0 for spherulites"
        steps = int(opts.get("steps", 500))
    else:
        raise SystemExit("case must be 'dendrites' or 'spherulite'")
    terms = [f"0.5*(1+tanh((sqrt((x-({cx:.5f}))**2+(y-({cy:.5f}))**2+(z-({cz:.5f}))**2)-0.05)/0.0280149))"
             for cx, cy, cz in centers]
    name = f"f4_{case}_n{n}_s{s:g}_chi{chi:g}_g{gamma:g}_fem"
    text = TEMPLATE.format(title=title, n=n, S=s, S_paper=s_paper, gamma=gamma, steps=steps,
                           snap=max(1, steps // 10), expected=expected, backend=backend,
                           gpu="true" if backend == "torch" else "false", name=name, half=half,
                           phi_expr="*".join(terms), c_s=C_S, dc=(1 - s) * (C_L - C_S),
                           centers=" ".join(f"{v:.5f}" for v in centers.ravel()),
                           noise=0.0025 / (DT * 0.1764), chi=chi)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.i").write_text(text)
    print(name)


if __name__ == "__main__":
    main()
