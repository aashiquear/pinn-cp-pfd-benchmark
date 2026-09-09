#!/usr/bin/env python3
"""Generate the graded all-hex void meshes for the Al2219 benchmark (gap G5).

The source paper never reports its element count or any mesh-convergence
evidence (data-sheet gaps C.2 and C.7), so this study has to establish both.
This script builds the refinement family the convergence study runs on.

WHY NOT ``box_with_void``: PINNACLE's bundled generator deletes background
elements by centroid test, which staircases the void surface and puts no extra
resolution near it. The void surface is where damage nucleates, so a staircased
boundary is a defect exactly where the answer lives. Hexaform's cubed-sphere
O-grid instead places nodes ON the analytic sphere and grades outward into the
Cartesian background, which is why the study meshes come from there and the
route is hexaform -> Exodus -> ``[Mesh] type = file``.

HOW THE VOID IS MADE: hexaform meshes the sphere as an *inclusion* carrying its
own block id. Dropping that block leaves a cavity whose boundary is the
conforming O-grid surface -- the void is the absence of the inclusion, not a
separate construction, so the surface keeps the geometric convergence the
inclusion mesh was verified to have.

GRADING TARGET, AND WHY IT NEEDS TWO CONTROLS, NOT ONE. The refined O-grid
macroblock must span the void radius plus three to five phase-field lengths, so
that a damage band of width ~l0 forming at the void surface sits inside the
refined zone rather than straddling the transition. With r = 5 um and l0 = 1 um
that is a refined radius of 8-10 um -- and setting only that is not enough.

The macroblock's radial layer count (``MacroBlockOptions.n_outer``) defaults to
2, so the shell from the void surface out to the macro-box is spanned by two
layers WHATEVER the background density. Raising the density then refines the
background and leaves the first element layer at the void surface almost
unchanged. Measured on an earlier version of this script that set only the
refined radius: the background element size fell 3.125 -> 0.781 um across the
family while the first layer at the void surface went only 3.15 -> 1.93 um,
giving 0.3-0.5 elements across l0 at every level. A convergence study run on
that family would have converged the background and never resolved the damage
band -- reporting mesh convergence for the one region where it had not been
achieved.

So the radial resolution is set explicitly here: ``n_outer`` layers graded by
``outer_grading`` toward the interface (its own docstring notes that >1
"clusters at the interface (where plastic strain and phase-field gradients
concentrate)"). The achieved first-layer thickness is measured and reported per
level rather than assumed, and the script refuses a family that fails the
requested points-per-l0 target.

``n_inner`` is left at its default: those layers are inside the sphere and are
carved away.

UNITS: micrometre, matching the source paper. The written mesh is in um and the
deck reads it with ``scale = 1.0``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from hexaform.assembly.builder import assemble
from hexaform.background.cartesian import CartesianGrid
from hexaform.core.mesh import HexMesh
from hexaform.core.quality import assess
from hexaform.core.topology import signed_volumes
from hexaform.io import exodus
from hexaform.packing.rsa import PackedInclusion
from hexaform.shapes.library import Sphere
from hexaform.templates.macroblock import MacroBlockOptions

EDGE_UM = 50.0
VOID_RADIUS_UM = 5.0
L0_UM = 1.0
REFINED_RADII_IN_L0 = (3.0, 5.0)


def carve_void(mesh: HexMesh, inclusion_block_min: int = 100) -> tuple[HexMesh, int]:
    """Return the mesh with the inclusion block removed, plus the count dropped.

    Nodes are NOT renumbered: an unused node is harmless on import and
    renumbering would invalidate the node sets hexaform wrote. What matters is
    that no element remains inside the void, so the cavity is genuinely empty
    rather than filled with a soft material.
    """
    keep = mesh.block_id < inclusion_block_min
    dropped = int((~keep).sum())
    carved = HexMesh(
        points=mesh.points,
        cells=mesh.cells[keep],
        block_id=mesh.block_id[keep],
        block_names={k: v for k, v in (mesh.block_names or {}).items() if k < inclusion_block_min},
    )
    return carved, dropped


def first_layer_thickness(mesh: HexMesh, centre: np.ndarray, radius: float) -> float:
    """Radial thickness of the element layer touching the void surface.

    This, not the background element size, is what the damage band has to be
    resolved by: the band forms at the void surface, so the relevant measure is
    how finely the first layer off that surface is divided.
    """
    r_node = np.linalg.norm(mesh.points - centre, axis=1)
    on_surface = np.abs(r_node - radius) < 1.0e-6
    touching = on_surface[mesh.cells].any(axis=1)
    if not touching.any():
        return float("nan")
    r_touch = np.linalg.norm(mesh.points[mesh.cells[touching]].reshape(-1, 3) - centre, axis=1)
    outward = r_touch[r_touch > radius + 1.0e-9]
    return float(np.median(outward) - radius) if outward.size else float("nan")


def build_level(
    density: int,
    refined_radius_um: float,
    n_outer: int,
    outer_grading: float,
) -> dict:
    """Mesh one refinement level and return its verification record."""
    grid = CartesianGrid.cube(EDGE_UM, density)
    h = float(grid.h if np.isscalar(grid.h) else np.max(grid.h))

    # O-grid macroblock half-width in cells, sized so the refined zone reaches
    # ``refined_radius_um``. Even ``k`` keeps the block centred on the sphere.
    k = 2 * int(np.ceil(refined_radius_um / h))
    centre_target = np.full(3, EDGE_UM / 2.0)
    lo = grid.point_to_index(centre_target) - k // 2
    lo = np.clip(lo, 1, np.array(grid.divisions) - k - 1)
    centre = grid.origin + (lo + k / 2.0) * h

    inclusion = PackedInclusion(Sphere(centre, VOID_RADIUS_UM), lo, k, 1)
    template = MacroBlockOptions(n_outer=n_outer, outer_grading=outer_grading)
    full = assemble(grid, [inclusion], template_options=template)

    # Geometric fidelity of the sphere BEFORE carving, which is the only stage
    # at which the inclusion volume can be compared with the analytic value.
    incl = full.block_id >= 100
    vol = float(signed_volumes(full.points, full.cells)[incl].sum())
    exact = 4.0 / 3.0 * np.pi * VOID_RADIUS_UM**3
    on_surface = np.abs(np.linalg.norm(full.points - centre, axis=1) - VOID_RADIUS_UM) < 1.0e-6
    node_err = float(
        np.abs(np.linalg.norm(full.points[on_surface] - centre, axis=1) - VOID_RADIUS_UM).max()
    )

    carved, dropped = carve_void(full)
    q = assess(carved)
    vols = signed_volumes(carved.points, carved.cells)

    return {
        "density": density,
        "h_background_um": h,
        "k_cells": k,
        "refined_radius_um": refined_radius_um,
        "refined_radius_in_l0": (refined_radius_um - VOID_RADIUS_UM) / L0_UM,
        "void_centre_um": [float(c) for c in centre],
        "centre_offset_um": float(np.linalg.norm(centre - centre_target)),
        "n_nodes": int(carved.points.shape[0]),
        "n_elements_total_before_carve": int(full.cells.shape[0]),
        "n_elements": int(carved.cells.shape[0]),
        "n_elements_carved": dropped,
        "min_scaled_jacobian": float(q.min_scaled_jacobian),
        "mean_scaled_jacobian": float(q.mean_scaled_jacobian),
        "p01_scaled_jacobian": float(q.p01_scaled_jacobian),
        "n_inverted": int(q.n_inverted),
        "n_below_threshold": int(q.n_below_threshold),
        "quality_threshold": float(q.threshold),
        "max_edge_ratio": float(q.max_edge_ratio),
        "min_element_volume_um3": float(vols.min()),
        "negative_jacobian_elements": int((vols <= 0.0).sum()),
        "inclusion_volume_um3": vol,
        "inclusion_volume_error": abs(vol - exact) / exact,
        "surface_node_radial_error_um": node_err,
        "n_outer_layers": int(n_outer),
        "outer_grading": float(outer_grading),
        "first_layer_um": first_layer_thickness(carved, centre, VOID_RADIUS_UM),
        "_mesh": carved,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--densities", type=int, nargs="+", default=[16, 24, 32, 48, 64])
    ap.add_argument("--refined-radius-um", type=float, default=8.0)
    ap.add_argument("--n-outer", type=int, nargs="+", default=[8],
                    help="radial layers between the void surface and the macro-box; "
                         "give one value per density to CO-REFINE the shell with the "
                         "background, or a single value to hold it fixed")
    ap.add_argument("--outer-grading", type=float, default=4.0,
                    help="largest/smallest radial layer; >1 clusters at the void surface")
    ap.add_argument("--min-points-per-l0", type=float, default=4.0,
                    help="required elements across l0 in the first layer at the void surface")
    ap.add_argument("--out-dir", default="meshes/al2219")
    args = ap.parse_args()

    lo_l0, hi_l0 = REFINED_RADII_IN_L0
    target_lo, target_hi = VOID_RADIUS_UM + lo_l0 * L0_UM, VOID_RADIUS_UM + hi_l0 * L0_UM
    if not target_lo <= args.refined_radius_um <= target_hi:
        raise SystemExit(
            f"--refined-radius-um {args.refined_radius_um} is outside the "
            f"r + 3..5*l0 window [{target_lo}, {target_hi}] this study requires; "
            "pass a value inside it or change the documented target deliberately."
        )

    # A refinement family has to refine EVERYWHERE. Holding n_outer fixed while
    # the background density rises leaves the first layer at the void surface
    # essentially unchanged (measured: 0.164-0.210 um across the family at
    # n_outer = 12), so the local fields where damage nucleates would be
    # computed on the same radial discretisation at every level and the study
    # would again be converging only the background.
    n_outer_list = list(args.n_outer)
    if len(n_outer_list) == 1:
        n_outer_list = n_outer_list * len(args.densities)
    elif len(n_outer_list) != len(args.densities):
        raise SystemExit(
            f"--n-outer takes one value or one per density; got "
            f"{len(n_outer_list)} for {len(args.densities)} densities."
        )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for density, n_outer in zip(args.densities, n_outer_list, strict=True):
        rec = build_level(density, args.refined_radius_um, n_outer, args.outer_grading)
        rec["points_across_l0"] = L0_UM / rec["first_layer_um"]
        mesh = rec.pop("_mesh")
        path = out_dir / f"void_r5_L50_d{density}.e"
        exodus.write(mesh, path)
        rec["path"] = str(path)
        rec["round_trip"] = exodus.round_trip_check(mesh, path)
        records.append(rec)
        print(
            f"d={density:3d}  h={rec['h_background_um']:.3f} um  k={rec['k_cells']:3d}  "
            f"elem={rec['n_elements']:8d}  minSJ={rec['min_scaled_jacobian']:.4f}  "
            f"negJ={rec['negative_jacobian_elements']}  "
            f"volerr={rec['inclusion_volume_error']:.3e}  "
            f"dr1={rec['first_layer_um']:.4f}um  l0/dr1={rec['points_across_l0']:.2f}  -> {path.name}"
        )

    drs = [r["first_layer_um"] for r in records]
    if len(drs) > 1 and not all(a > b for a, b in zip(drs, drs[1:])):
        print(
            "\nWARNING: the first-layer thickness at the void surface is not "
            f"strictly decreasing across the family ({[round(v, 4) for v in drs]}). "
            "The near-void region is then not being refined in step with the "
            "background, so a convergence study on this family would measure "
            "background convergence only. Scale --n-outer with the density."
        )

    # The resolution gate applies to the FINEST level, not to every level. A
    # convergence study needs levels that are demonstrably too coarse -- that is
    # how convergence becomes observable -- so an under-resolved coarse level is
    # a feature. What would invalidate the family is the finest level failing to
    # resolve the band, because then nothing in it is converged.
    finest = max(records, key=lambda r: r["density"])
    if finest["points_across_l0"] < args.min_points_per_l0:
        raise SystemExit(
            f"REFUSED: the finest level (d = {finest['density']}) gives only "
            f"{finest['points_across_l0']:.2f} elements across l0 = {L0_UM:g} um at "
            f"the void surface, below the required {args.min_points_per_l0:g}. "
            "Nothing in this family would be resolved; raise --n-outer for it."
        )
    under = [r for r in records if r["points_across_l0"] < args.min_points_per_l0]
    for r in under:
        r["under_resolved_at_l0"] = True
    if under:
        levels = ", ".join(f"d={r['density']} ({r['points_across_l0']:.2f})" for r in under)
        print(
            f"\nDeliberately under-resolved at l0 = {L0_UM:g} um "
            f"(fewer than {args.min_points_per_l0:g} elements across the band): {levels}. "
            "These anchor the coarse end of the convergence study and must be "
            "reported as under-resolved wherever their local fields are quoted."
        )

    summary = out_dir / "mesh_family.json"
    summary.write_text(json.dumps(records, indent=1))
    print(f"\nwrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
