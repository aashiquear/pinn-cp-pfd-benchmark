#!/usr/bin/env python3
"""Write decks of case E1: planar crystallization front of the binary Cu-Ni model in one dimension,
from partitioning to diffusionless (solute-trapping) growth with rising supersaturation.

Usage: python3 scripts/cryst_trapping_decks.py S=<S,...> [cells=3000] [steps=40000]
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "decks/crystallization/trapping"
C_L, C_S, DX = 0.466219, 0.399112, 6.25e-3

TEMPLATE = """# Part IV, case E1 (FE route): planar crystallization front in one dimension at supersaturation
# S = {S:g}: partitioning of the solute at the front against trapping, binary Cu-Ni model.
# After the transition from diffusion-controlled to diffusionless growth shown in Fig. 8(a) of
# (granasy 2021 biomineralization) for the same model (its PF1).
#
# PHYSICS. Binary ideal-solution model of case C2 (granasy 2005 spherulites, Eq. 2.3-2.4),
# phi = 1 liquid, phi = 0 crystal, c = Cu fraction, one space dimension:
#
#     phi_t = phi_xx - a_w(c) g'(phi) - p'(phi) b(c)
#     c_t   = d/dx{{ p(phi) [ c_x - c (1 - c) (q_f p'(phi) + q_w g'(phi)) phi_x ] }}
#
# The crystal rejects Cu (c_S < c_L). Below S = 1 the rejected solute has to diffuse away and
# the front slows as t^(-1/2); above S = 1 the liquid already has less Cu than the solidus and
# the front runs at a constant speed, freezing the liquid composition in (solute trapping).
#
# UNITS and PARAMETERS: as case C2 (lengths in xi = 2.1e-4 cm, time in xi^2/D_l; provenance in
# docs/14_phase_field_crystallization_plan.md, section 7.2). c_inf = c_L - S (c_L - c_S) = {c_inf:.6f}.
# assumed: the initial crystal slab of 30 cells, the run length, no noise.
#
# HOW IT IS SOLVED. `pf_fem_explicit` on {cells} linear elements, dx = 6.25e-3, dt = 4.75e-6,
# {steps} steps; no flux at both ends.
#
# EXPECTED RESULT. See the Part IV report (binary model, "What to expect"): front position
# against time, composition of the crystal that formed and the solute peak ahead of the front.

[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = 0
  output_dir = results/crystallization/trapping/{name}
[]

[Problem]
  type          = composed
  dimension     = 1
  finite_strain = false
[]

[Mesh]
  type      = line_mesh
  domain    = '0.0 {length:g}'
  divisions = {cells}
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
[]

[Functions]
  # crystal for x < 30 dx = 0.1875
  [./seed_phi]
    kind       = parsed
    expression = 0.5*(1+tanh((x-0.1875)/0.0280149))
  []
  [./seed_c]
    kind       = parsed
    expression = {c_s}+{dc:.6f}*0.5*(1+tanh((x-0.1875)/0.0280149))
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
[]

[Executioner]
  type            = pf_fem_explicit
  n_steps         = {steps}
  dt              = 4.75e-6
  snapshots_every = {snap}
  energy_every    = 0
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'phi c'
[]
"""


def main() -> None:
    opts = dict(a.split("=") for a in sys.argv[1:])
    cells, steps = int(opts.get("cells", 3000)), int(opts.get("steps", 40000))
    OUT.mkdir(parents=True, exist_ok=True)
    for s in [float(v) for v in opts["S"].split(",")]:
        c_inf = C_L - s * (C_L - C_S)
        name = f"e1_planar_s{s:g}_fem"
        (OUT / f"{name}.i").write_text(TEMPLATE.format(S=s, c_inf=c_inf, cells=cells, steps=steps, name=name,
                                                       length=cells * DX, c_s=C_S, dc=c_inf - C_S,
                                                       snap=max(1, steps // 20)))
        print(name)


if __name__ == "__main__":
    main()
