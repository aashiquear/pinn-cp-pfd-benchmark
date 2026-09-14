#!/usr/bin/env python3
"""Relax a Voronoi polycrystal into a natural morphology by grain growth.

A Voronoi tessellation gives straight boundaries, flat-faced polygons and
triple junctions at arbitrary angles -- none of which a real polycrystal has.
The standard remedy, which is what this script does, is to use the Voronoi
structure as the INITIAL CONDITION for a curvature-driven phase-field grain
growth simulation and take the morphology after the boundary network has
relaxed: boundaries curve, triple junctions approach 120 degrees, and small
many-sided grains shrink.

SCOPE. This is a MORPHOLOGY GENERATOR, not a simulation of Al2219 grain
growth. Curvature-driven growth is self-similar, so the material constants set
the timescale and the boundary width, not the shapes; any FCC metal gives
comparable morphology. Nothing here should be reported as Al2219 kinetics.

RESOLUTION. Growth runs on its OWN fine grid, not the mechanics mesh. The
diffuse boundary needs h <= w_GB/3 (gb_evolution's own requirement) while
w_GB must stay well below the grain size; at 32 grains that needs h ~ 0.5 um
against the mechanics mesh's 2.08 um, a factor of four. Morphology and
mechanics resolutions are independent, so the relaxed grain IDs are mapped
onto the mechanics element centroids afterwards. The plate is one element
thick, so in-plane 2-D growth is the exact analogue of growth on the
pseudo-3-D domain rather than an approximation of it.

UNITS. Normalized: sigma = 1, M_GB = 1, lengths in micrometers. Only the
morphology is used, so the absolute timescale is irrelevant to the product;
the mapping to the reference deck's physical values is recorded in the output
for traceability. Note the reference deck's ``Q = 0.23`` is in eV, not J/mol
-- reading it as J/mol makes the mobility 377x too large.

GRAIN IDS. With graph colouring several NON-ADJACENT grains share an order
parameter, so an argmax over order parameters gives the COLOUR, not the grain.
Grain identity therefore comes from connected-component labelling within each
colour, which is what MOOSE's FeatureFloodCount does internally.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

STUDY = Path(__file__).resolve().parent.parent
PLATE = (100.0, 50.0)          # um, the plate's in-plane extent


def build_deck(*, n_grains: int, n_ops: int, seeds: np.ndarray, div: tuple[int, int],
               w_gb: float, n_steps: int, dt: float, out_dir: Path) -> Path:
    """A 2-D grain-growth deck seeded from an existing Voronoi structure."""
    flat = " ".join(f"{v:.6f}" for v in seeds.reshape(-1))
    deck = f"""# Grain-growth morphology relaxation, seeded from the study's Voronoi
# structure so the relaxed polycrystal is demonstrably the SAME polycrystal
# evolved, not a different random one.
#
# Normalized units: sigma = 1, M_GB = 1, lengths in micrometers.
# w_GB = {w_gb:.4f} um, chosen as grain_width/8 with the grid at h <= w_GB/3.

[Global]
  backend          = numpy
  use_gpu          = false
  precision        = float64
  seed             = 0
  output_dir       = {out_dir / "solver"}
[]

[Mesh]
  type             = box_mesh
  dim              = 2
  domain           = '0 0 {PLATE[0]:.6g} {PLATE[1]:.6g}'
  divisions        = '{div[0]} {div[1]}'
[]

[Variables]
  [./eta]
    family         = fem
    components     = {n_ops}
    # NOTE: ``order`` here is the Gauss QUADRATURE rule, not the polynomial
    # order of the element. Only order = 2 is implemented; order = 1 raises
    # from deep inside quadrature construction.
    order          = 2
  []
[]

[Materials]
  [./gb]
    kind           = gb_evolution
    M_GB           = 1.0
    w_GB           = {w_gb:.6g}
    sigma          = 1.0
    gamma          = 1.5
  []
[]

[Kernels]
  [./deta_dt]
    type           = time_derivative
    variable       = eta
    coefficient    = 1.0
  []
  [./bulk]
    type           = ac_bulk
    variable       = eta
    material       = gb
    spatial_dim    = 2
  []
  [./interface]
    type           = ac_interface
    variable       = eta
    material       = gb
    spatial_dim    = 2
  []
[]

[ICs]
  [./grains]
    variable       = eta
    kind           = polycrystal_coloring
    n_grains       = {n_grains}
    n_order_params = {n_ops}
    # The study's own Voronoi seeds, so the relaxed structure is the same
    # polycrystal the single-crystal orientation sets were built on.
    seeds          = '{flat}'
    periodic       = false
    int_width      = {w_gb:.6g}
    inside         = 1.0
    outside        = 0.0
  []
[]

[Utilities]
  [./count]
    kind           = feature_flood_count
    threshold      = 0.5
  []
  # NOTE: ``grain_tracker`` is deliberately NOT used here. Its ``n_order_params``
  # is the size of an OP POOL that bounds how many grains can coexist, so with
  # grains coloured onto fewer order parameters than there are grains it fails
  # with "no free order parameter left for a new grain". The tracker's persistent
  # identity and remap bookkeeping are for long coarsening runs; this is a short
  # relaxation at constant grain count, and ``feature_flood_count`` already
  # reports the connected-component count that is the stopping criterion.
[]

[Executioner]
  type             = pf_fem_transient
  scheme           = bdf1
  n_steps          = {n_steps}
  dt               = {dt:.6g}
  linear_solver    = direct
[]

[Postprocessors]
  [./gb_area]
    type            = grain_boundary_area
    fields          = 'eta'
    grains_per_side = 2
    delta_eta       = 1.0
  []
[]

[Output]
  formats          = 'csv json'
  fields           = 'eta'
[]
"""
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / "deck.i"
    p.write_text(deck)
    return p


def label_grains(eta: np.ndarray, div: tuple[int, int], threshold: float = 0.5):
    """Grain ids by connected-component labelling within each colour.

    ``eta`` is ``(n_nodes, n_ops)`` on the structured grid. The colour a node
    belongs to is its argmax; grains are the connected components of each
    colour's own indicator, because colouring deliberately puts several
    non-adjacent grains on one order parameter.
    """
    nx, ny = div[0] + 1, div[1] + 1
    col = np.argmax(eta, axis=1).reshape(nx, ny, order="F")
    val = np.max(eta, axis=1).reshape(nx, ny, order="F")
    gid = np.full((nx, ny), -1, dtype=int)
    nxt = 0
    for k in range(eta.shape[1]):
        mask = (col == k) & (val >= threshold)
        if not mask.any():
            continue
        lab, n = ndimage.label(mask)
        for i in range(1, n + 1):
            gid[lab == i] = nxt
            nxt += 1
    return gid, nxt


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, required=True)
    ap.add_argument("--density", type=int, default=3)
    ap.add_argument("--n-ops", type=int, default=8)
    ap.add_argument("--steps", type=int, default=40)
    ap.add_argument("--wgb-over-grain", type=float, default=1 / 8.0)
    ap.add_argument("--h-over-wgb", type=float, default=1 / 3.0)
    ap.add_argument("--out", default="results/graingrowth")
    args = ap.parse_args()

    sets = {r["n_grains"]: r for r in
            json.loads((STUDY / f"meshes/plate/polycrystal_sets_d{args.density}.json").read_text())}
    rec = sets[args.grains]
    gw = rec["grain_width_um_median"]
    w_gb = gw * args.wgb_over_grain
    h = w_gb * args.h_over_wgb
    div = (int(np.ceil(PLATE[0] / h)), int(np.ceil(PLATE[1] / h)))

    z = np.load(STUDY / f"meshes/plate/poly_d{args.density}_n{args.grains}.npz")
    seeds = np.asarray(z["seeds"], dtype=float)

    # Characteristic curvature-driven time: grain_width^2 / (M_GB sigma) with
    # both normalized to 1. Relax a modest fraction of it -- enough to curve
    # the boundaries and settle the junctions, not enough to coarsen the grain
    # count away, which is the stopping criterion.
    tau = gw * gw
    dt = 0.02 * tau / args.steps * 10.0

    out = STUDY / args.out / f"n{args.grains}"
    deck = build_deck(n_grains=args.grains, n_ops=args.n_ops, seeds=seeds, div=div,
                      w_gb=w_gb, n_steps=args.steps, dt=dt, out_dir=out)
    print(json.dumps({"n_grains": args.grains, "grain_width_um": gw, "w_GB_um": w_gb,
                      "h_um": h, "divisions": list(div), "n_nodes": (div[0]+1)*(div[1]+1),
                      "tau_um2": tau, "dt": dt, "n_steps": args.steps,
                      "deck": str(deck)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
