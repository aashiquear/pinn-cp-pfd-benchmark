#!/usr/bin/env python3
"""Part IV, stage D (phase-field crystal): deck generator.

Usage: python3 scripts/cryst_pfc_decks.py
Writes decks/crystallization/pfc/*.i. Every deck is standalone (`pinnacle run deck.i`).

Model (granasy 2019 nucleation, Eq. 30 and 36 without noise; elder 2004 plastic):
    F = int { psi/2 [-eps + (1 + lap)^2] psi + psi^4/4 } dA,   psi_t = lap(dF/dpsi)
Lattice of the triangular phase in the one-mode approximation: wavenumber of the three
density waves 1, lattice constant a = 4 pi / sqrt(3) = 7.2552, rows spaced a sqrt(3)/2 = 2 pi.
A box whose sides are mirror lines of the lattice holds it without strain:
    Lx = n a / 2 (rows of particles along x),  Ly = m 2 pi / ... (see `box`).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DECKS = ROOT / "decks/crystallization/pfc"
A = 4.0 * np.pi / np.sqrt(3.0)     # lattice constant
Q = np.sqrt(3.0) / 2.0             # x wavenumber of the one-mode pattern


def tri_amplitude(eps: float, psi0: float) -> float:
    """One-mode amplitude of psi0 + A [cos(Q x) cos(y/2) - cos(y)/2] (elder 2004 plastic)."""
    root = 15.0 * eps - 36.0 * psi0**2
    if root < 0:
        return 0.0
    return 4.0 / 5.0 * (psi0 - np.sqrt(root) / 3.0) if psi0 <= 0 else 4.0 / 5.0 * (psi0 + np.sqrt(root) / 3.0)


def pattern(amp: float, angle: float = 0.0, x0: float = 0.0, y0: float = 0.0) -> str:
    """Deck expression of amp [cos(Q x') cos(y'/2) - cos(y')/2], rotated by `angle` (rad)."""
    c, s = np.cos(angle), np.sin(angle)
    xr = f"({c:.12g}*(x-{x0:.10g})+{s:.12g}*(y-{y0:.10g}))"
    yr = f"({-s:.12g}*(x-{x0:.10g})+{c:.12g}*(y-{y0:.10g}))"
    return f"{amp:.10g}*(cos({Q:.15g}*{xr})*cos(0.5*{yr})-0.5*cos({yr}))"


HEAD = """# Part IV, case {case} (FE route): {title}
#
# PHYSICS. Phase-field crystal model (granasy 2019 nucleation, Eq. 30 and 36, no noise;
# elder 2004 plastic): psi is the reduced particle density,
#
#     F = int {{ psi/2 [-eps + (1 + lap)^2] psi + psi^4/4 }} dA,
#     dpsi/dt = lap(mu),   mu = dF/dpsi = -eps psi + psi^3 + (1 + lap)^2 psi.
#
# The homogeneous state is the liquid; a triangular lattice of density peaks is the crystal
# (in two dimensions), stripes the third phase. All quantities are dimensionless.
#
# HOW IT IS SOLVED. Sixth order in space: mixed form with the chemical potential mu and
# nu = lap(psi) on bilinear elements (PINNACLE ADR 0046), implicit steps with Newton's
# method. The box sides are mirror lines (natural conditions n.grad psi = n.grad nu =
# n.grad mu = 0).
#
{notes}
"""

BODY = """
[Global]
  backend    = numpy
  use_gpu    = false
  precision  = float64
  seed       = {seed}
  output_dir = results/crystallization/pfc/{name}
[]

[Problem]
  type          = composed
  dimension     = 2
  finite_strain = false
[]

[Mesh]
  type      = box_mesh
  dim       = 2
  domain    = '0 0 {lx:.12g} {ly:.12g}'
  divisions = '{nx} {ny}'
[]

[Variables]
  [./psi]
    family     = fem
    components = 1
    order      = 2
  []
  [./mu]
    family     = fem
    components = 1
    order      = 2
  []
  [./nu]
    family     = fem
    components = 1
    order      = 2
  []
[]
{ics}
[Kernels]
  [./dpsi_dt]
    type     = coupled_time_derivative
    variable = psi
  []
  [./transport]
    type     = split_ch_wres
    variable = mu
    mobility = 1.0
  []
  [./potential]
    type               = sh_potential
    variable           = psi
    chemical_potential = mu
    laplacian          = nu
    epsilon            = {eps}
  []
  [./laplacian]
    type     = weak_laplacian
    variable = nu
    of       = psi
  []
[]

[Executioner]
  type            = pf_fem_transient
  scheme          = {scheme}
  n_steps         = {n_steps}
  dt              = {dt}
  max_newton_iter = 25
  newton_tol      = 1.0e-8
  linear_solver   = gmres_lu     # GMRES preconditioned by a reused LU factorization
  snapshots_every = {every}
[]

[Output]
  formats = 'vtk csv json'
  fields  = 'psi'
[]
"""


def function_ic(expr: str) -> str:
    return f"""
[Functions]
  [./start]
    kind       = parsed
    expression = {expr}
  []
[]

[ICs]
  [./psi0]
    variable = psi
    kind     = from_function
    function = start
  []
[]
"""


def random_ic(psi0: float, amp: float) -> str:
    return f"""
[ICs]
  [./psi0]
    variable = psi
    kind     = uniform_random
    minimum  = {psi0 - amp:.10g}
    maximum  = {psi0 + amp:.10g}
  []
[]
"""


def write(name: str, case: str, title: str, notes: str, **kw) -> None:
    DECKS.mkdir(parents=True, exist_ok=True)
    kw.setdefault("seed", 0)
    kw.setdefault("scheme", "bdf1")
    text = HEAD.format(case=case, title=title, notes=notes.rstrip()) + BODY.format(name=name, **kw)
    (DECKS / f"{name}.i").write_text(text)


def box(n_a: int, n_rows: int) -> tuple[float, float]:
    """Mirror box: n_a half lattice constants along x, n_rows half row periods (2 pi) along y."""
    return n_a * A / 2.0, n_rows * 2.0 * np.pi


#: coexistence densities of the liquid and the triangular crystal (relaxed cells of
#: scripts/cryst_pfc_reference.py diagram; data/reference/cryst_pfc_phase_diagram.csv)
LIQUIDUS = {0.1: -0.20681, 0.2: -0.30109, 0.3: -0.37931, 0.4: -0.44943, 0.5: -0.51409}
SOLIDUS = {0.1: -0.19611, 0.2: -0.26764, 0.3: -0.31656, 0.4: -0.35412, 0.5: -0.38496}


def seed_amplitude(eps: float) -> float:
    """One-mode amplitude of the crystal at its coexistence (solidus) density."""
    return tri_amplitude(eps, SOLIDUS[eps])


def seed_expr(eps: float, psi0: float, radius: float, angle: float = 0.0, x0: float = 0.0,
              y0: float = 0.0, width: float = 3.0) -> str:
    """psi0 + envelope * one-mode pattern; the mean of the pattern is zero."""
    amp = seed_amplitude(eps)
    env = f"0.5*(1-tanh((sqrt((x-{x0:.10g})**2+(y-{y0:.10g})**2)-{radius:.10g})/{width:.6g}))"
    return f"{psi0:.10g}+{env}*{pattern(amp, angle, x0, y0)}"


def main() -> None:
    # D1a: the triangular crystal in a commensurate mirror box, three mesh sizes
    eps, psi0 = 0.3, -0.25
    for per_a in (8, 12, 16):
        lx, ly = box(4, 2)          # 2 a by 4 pi: two lattice constants, two row periods
        nx, ny = 2 * per_a, int(round(ly / (A / per_a)))
        write(f"d1_triangular_eps{eps}_psi{psi0}_h{per_a}_fem", "D1",
              f"the triangular crystal at eps = {eps}, mean density {psi0}, {per_a} elements per lattice constant.",
              f"""# EXPECTED RESULT. The one-mode pattern relaxes to the crystal: amplitude and free-energy
# density (`energy` of solver_history.csv divided by the box area {lx * ly:.6f}) against the
# spectral reference scripts/cryst_pfc_reference.py; the error falls as h^2.""",
              lx=lx, ly=ly, nx=nx, ny=ny, eps=eps, n_steps=60, dt=2.0, every=10,
              ics=function_ic(f"{psi0}+{pattern(tri_amplitude(eps, psi0))}"))
    # D1b: which phase forms from a random state, series in the mean density at eps = 0.3
    for psi0 in (0.0, -0.06, -0.12, -0.2, -0.3, -0.35, -0.4):
        lx, ly = box(16, 7)
        write(f"d1_phase_eps0.3_psi{psi0}_fem", "D1",
              f"which phase forms from a random state at eps = 0.3, mean density {psi0}.",
              """# EXPECTED RESULT (one-mode phase diagram, elder 2004 plastic; values in the report):
# stripes around psi0 = 0, the triangular phase at intermediate |psi0|, the liquid (uniform
# field) beyond the liquidus.""",
              lx=lx, ly=ly, nx=96, ny=int(round(ly / (lx / 96))), eps=0.3, n_steps=150, dt=2.0, every=10,
              ics=random_ic(psi0, 0.05))
    # D2: growth of a crystal from a seed in the metastable liquid; a quarter of the crystal
    # (the seed sits in the corner, x = 0 and y = 0 are mirror lines of the lattice)
    for eps, psi0, per_a, dt, t_end, tag, radius in GROWTH:
        grow_deck(eps, psi0, per_a, dt, t_end, tag, radius=radius)
    polycrystal()


#: (eps, psi0, elements per lattice constant, dt, end time, name tag, seed radius)
#: Supersaturation S = (psi0 - liquidus)/(solidus - liquidus) with the coexistence densities of
#: data/reference/cryst_pfc_phase_diagram.csv (relaxed cells).
GROWTH = [
    (0.3, -0.34, 10, 2.0, 1500.0, "", 12.0),            # first case, S = 0.63
    (0.3, -0.34, 10, 1.0, 1500.0, "_dt1", 12.0),        # time-step check
    (0.3, -0.34, 14, 2.0, 1500.0, "_h14", 12.0),        # mesh check
    # series in the mean density at eps = 0.3 (S = 0.31, 0.47, 0.79, 0.95)
    (0.3, -0.36, 10, 2.0, 2000.0, "", 12.0),
    (0.3, -0.35, 10, 2.0, 2000.0, "", 12.0),
    (0.3, -0.36, 10, 2.0, 3000.0, "_r24", 24.0),        # a larger seed at low supersaturation
    (0.3, -0.35, 10, 2.0, 3000.0, "_r24", 24.0),
    (0.3, -0.33, 10, 2.0, 1500.0, "", 12.0),
    (0.3, -0.32, 10, 2.0, 1500.0, "", 12.0),
    # series in eps at S = 0.63, seed radius 16
    (0.1, -0.2001, 10, 4.0, 6000.0, "_r16", 16.0),
    (0.1, -0.2001, 10, 4.0, 8000.0, "_r32", 32.0),
    (0.2, -0.2800, 10, 2.0, 3000.0, "_r16", 16.0),
    (0.3, -0.3400, 10, 2.0, 1500.0, "_r16", 16.0),
    (0.4, -0.3894, 10, 2.0, 1500.0, "_r16", 16.0),
    (0.5, -0.4328, 10, 2.0, 1500.0, "_r16", 16.0),
]


def grow_deck(eps: float, psi0: float, per_a: int, dt: float, t_end: float, tag: str = "",
              n_a: int = 26, n_rows: int = 15, radius: float = 12.0) -> None:
    lx, ly = box(n_a, n_rows)
    nx = n_a * per_a // 2
    ny = int(round(ly / (lx / nx)))
    steps = int(round(t_end / dt))
    write(f"d2_growth_eps{eps}_psi{psi0}{tag}_fem", "D2",
          f"growth of a triangular crystal from a seed, eps = {eps}, mean density {psi0}.",
          f"""# SET-UP. A quarter of the crystal: box {lx:.3f} x {ly:.3f} ({n_a / 2:g} lattice constants by
# {n_rows} row periods), {per_a} elements per lattice constant, seed of radius {radius:g} (one-mode
# pattern of amplitude {seed_amplitude(eps):.4f}, that of the crystal at its coexistence density, under a tanh
# envelope of width 3) in the corner. The density is conserved: the
# crystal is denser than the liquid, so a depletion zone forms ahead of the front and the
# growth stops at the lever-rule fraction.
#
# EXPECTED RESULT. The crystalline fraction against time follows the spectral reference
# (scripts/cryst_pfc_reference.py grow {eps} {psi0} {n_a} {n_rows} {radius:g} {t_end:g} 0.125 {seed_amplitude(eps):.6f}); values in the report.""",
          lx=lx, ly=ly, nx=nx, ny=ny, eps=eps, n_steps=steps, dt=dt, every=max(steps // 50, 1),
          ics=function_ic(seed_expr(eps, psi0, radius)))


#: seeds of the polycrystal: (x, y, rotation in degrees); the lattice has six-fold symmetry,
#: so 30 degrees is the largest misorientation
SEEDS = [(30.0, 28.0, 0.0), (98.0, 34.0, 12.0), (62.0, 70.0, 27.0), (24.0, 104.0, 18.0), (104.0, 100.0, -9.0)]


def polycrystal() -> None:
    eps, psi0, per_a = 0.3, -0.32, 10
    lx, ly = box(36, 21)
    nx = 36 * per_a // 2
    ny = int(round(ly / (lx / nx)))
    amp = seed_amplitude(eps)
    terms = []
    for x0, y0, deg in SEEDS:
        env = f"0.5*(1-tanh((sqrt((x-{x0:g})**2+(y-{y0:g})**2)-12)/3))"
        terms.append(f"{env}*{pattern(amp, np.deg2rad(deg), x0, y0)}")
    write("d3_polycrystal_eps0.3_psi-0.32_fem", "D3",
          "five seeds of different lattice orientation grow and impinge, eps = 0.3, mean density -0.32.",
          f"""# SET-UP. Box {lx:.3f} x {ly:.3f} (18 lattice constants), {per_a} elements per lattice constant;
# seeds of radius 12 at {", ".join(f"({x:g}, {y:g}) rotated {d:g} deg" for x, y, d in SEEDS)}.
# The box sides are mirror lines: a grain that reaches a side meets its mirror image there.
#
# EXPECTED RESULT. The grains grow until they meet; because their lattices do not match,
# the boundaries between them are rows of dislocations (particles with five or seven
# neighbours) whose spacing falls with growing misorientation. Counts in the report.""",
          lx=lx, ly=ly, nx=nx, ny=ny, eps=eps, n_steps=1000, dt=2.0, every=20,
          ics=function_ic(f"{psi0}+" + "+".join(terms)))


if __name__ == "__main__":
    main()
