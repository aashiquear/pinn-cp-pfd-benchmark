#!/usr/bin/env python3
"""The standard per-run output contract. Every harness in this study calls it.

Two products per run, both required:

* ``results.csv`` -- one row per load step: step index and load factor,
  damage volume for c > 0.85, volume-averaged stress, applied strain at the
  loaded boundary, and volume-averaged strain.
* ``fields/`` -- a VTK series (one ``.vtu`` per step plus a ``.pvd``) carrying
  displacement, damage, stress and strain on the mesh, so a quantity nobody
  anticipated can be extracted later without re-running. Re-running is the
  expensive branch here: a coupled solve is hours, reading a ``.vtu`` is
  milliseconds.

WHY THE FULL TENSORS ARE STORED. "Average stress" is ambiguous -- axial, von
Mises, hydrostatic and the full tensor are all reasonable readings, and which
one matters depends on the question being asked months from now. The CSV
carries all six independent components of the volume-averaged stress and
strain, plus von Mises and the axial component as named columns. Storing the
tensor costs six floats a row; guessing wrong costs a re-run.

DAMAGE VOLUME IS ELEMENT-BASED, NOT NODE-BASED. Damage is a nodal field, so a
node count is mesh-dependent and has no volume attached. An element counts
toward the damaged volume when the mean of its eight nodal damage values
exceeds the threshold, and it contributes its own volume -- consistent with
the element-nodal-average criterion agreed for fracture initiation, and
mesh-independent in the limit.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

DAMAGE_THRESHOLD = 0.85


def hex_volumes(nodes: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Volume of every HEX8 element, by 2x2x2 Gauss quadrature on the Jacobian.

    Not the cross-product shortcut: that is exact only for a parallelepiped,
    and the O-grid elements around the hole are not. Integrating |J| over the
    reference cube is exact for a trilinear map.
    """
    g = 1.0 / np.sqrt(3.0)
    sgn = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)
    vol = np.zeros(cells.shape[0], dtype=float)
    xe = nodes[cells]                                     # (ne, 8, 3)
    for a in (-g, g):
        for b in (-g, g):
            for c in (-g, g):
                xi = np.array([a, b, c])
                # dN/dxi for the trilinear hex, shape (8, 3)
                dN = np.empty((8, 3))
                for i in range(8):
                    s = sgn[i]
                    dN[i, 0] = 0.125 * s[0] * (1 + s[1] * xi[1]) * (1 + s[2] * xi[2])
                    dN[i, 1] = 0.125 * (1 + s[0] * xi[0]) * s[1] * (1 + s[2] * xi[2])
                    dN[i, 2] = 0.125 * (1 + s[0] * xi[0]) * (1 + s[1] * xi[1]) * s[2]
                J = np.einsum("eai,aj->eij", xe, dN)      # (ne, 3, 3)
                vol += np.abs(np.linalg.det(J))           # weight 1 per point
    return vol


def _elem_tensor(field: np.ndarray, n_elem: int) -> np.ndarray | None:
    """Reduce a per-quadrature-point tensor field to one tensor per element."""
    if field is None:
        return None
    a = np.asarray(field, dtype=float)
    if a.size % 9 != 0:
        return None
    t = a.reshape(-1, 3, 3)
    if t.shape[0] % n_elem != 0:
        return None
    return t.reshape(n_elem, -1, 3, 3).mean(axis=1)        # (ne, 3, 3)


def _von_mises(t: np.ndarray) -> float:
    dev = t - np.eye(3) * np.trace(t) / 3.0
    return float(np.sqrt(1.5 * np.sum(dev * dev)))


def collect_step_rows(
    snapshots: list[dict[str, Any]],
    nodes: np.ndarray,
    cells: np.ndarray,
    *,
    load_axis: int = 0,
    gauge_length: float,
    boundary_nodes: np.ndarray,
    threshold: float = DAMAGE_THRESHOLD,
) -> list[dict[str, float]]:
    """One CSV row per snapshot."""
    vol = hex_volumes(nodes, cells)
    total_vol = float(vol.sum())
    rows: list[dict[str, float]] = []
    for sn in snapshots:
        row: dict[str, float] = {
            "step": int(sn.get("step", -1)),
            "load_factor": float(sn.get("step_factor", np.nan)),
            "total_volume": total_vol,
        }

        d = sn.get("damage")
        if d is not None:
            dn = np.asarray(d, dtype=float).ravel()
            d_elem = dn[cells].mean(axis=1)
            dmg = vol[d_elem > threshold].sum()
            row["damage_volume"] = float(dmg)
            row["damage_volume_fraction"] = float(dmg / total_vol)
            row["peak_damage_nodal"] = float(dn.max())
            row["peak_damage_element"] = float(d_elem.max())
            # Volume-weighted mean damage: the diffuse counterpart of the
            # thresholded volume, and the pair distinguishes a localized band
            # from damage smeared over the whole domain.
            row["mean_damage"] = float(np.sum(d_elem * vol) / total_vol)

        u = sn.get("displacement")
        if u is not None:
            un = np.asarray(u, dtype=float).reshape(len(nodes), -1)
            # Applied strain at the loaded boundary, from the DISPLACEMENT the
            # solver produced rather than from the prescribed value, so a
            # boundary condition that failed to apply shows up here.
            row["boundary_strain"] = float(un[boundary_nodes, load_axis].mean() / gauge_length)
            row["boundary_displacement"] = float(un[boundary_nodes, load_axis].mean())

        for name, key in (("stress", "stress"), ("strain", "strain")):
            te = _elem_tensor(sn.get(key), len(cells))
            if te is None:
                continue
            avg = np.einsum("eij,e->ij", te, vol) / total_vol
            for i, j, lab in ((0, 0, "xx"), (1, 1, "yy"), (2, 2, "zz"),
                              (0, 1, "xy"), (0, 2, "xz"), (1, 2, "yz")):
                row[f"avg_{name}_{lab}"] = float(avg[i, j])
            row[f"avg_{name}_vm"] = _von_mises(avg)
            row[f"avg_{name}_axial"] = float(avg[load_axis, load_axis])
        rows.append(row)
    return rows


def write_vtk_series(
    out_dir: Path,
    snapshots: list[dict[str, Any]],
    nodes: np.ndarray,
    cells: np.ndarray,
) -> list[str]:
    """Write one ``.vtu`` per step plus a ``.pvd`` index."""
    import meshio

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for sn in snapshots:
        step = int(sn.get("step", len(written)))
        point_data: dict[str, np.ndarray] = {}
        cell_data: dict[str, list[np.ndarray]] = {}
        d = sn.get("damage")
        if d is not None:
            point_data["damage"] = np.asarray(d, dtype=float).ravel()
        u = sn.get("displacement")
        if u is not None:
            point_data["displacement"] = np.asarray(u, dtype=float).reshape(len(nodes), -1)
        for key in ("stress", "strain"):
            te = _elem_tensor(sn.get(key), len(cells))
            if te is None:
                continue
            # Element averages, written as cell data: the raw field is per
            # quadrature point and VTK has nowhere to put that.
            cell_data[key] = [te.reshape(len(cells), 9)]
            cell_data[f"{key}_vm"] = [
                np.array([_von_mises(t) for t in te], dtype=float)
            ]
        m = meshio.Mesh(points=nodes, cells=[("hexahedron", cells)],
                        point_data=point_data, cell_data=cell_data)
        path = out_dir / f"step_{step:04d}.vtu"
        meshio.write(path, m)
        written.append(path.name)

    pvd = ['<?xml version="1.0"?>',
           '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
           "  <Collection>"]
    for sn, name in zip(snapshots, written):
        pvd.append(f'    <DataSet timestep="{float(sn.get("step_factor", 0.0)):.6f}" '
                   f'file="{name}"/>')
    pvd += ["  </Collection>", "</VTKFile>"]
    (out_dir / "fields.pvd").write_text("\n".join(pvd))
    # Opening an individual .vtu shows ONE frame with no way to advance time,
    # which reads as "only the last step was written". The .pvd is the time
    # series; say so next to the files rather than in a document nobody has
    # open at the time.
    (out_dir / "README.txt").write_text(
        "Open fields.pvd (NOT an individual step_XXXX.vtu) to animate.\n"
        "Each step_XXXX.vtu is one load step; fields.pvd indexes them in time.\n"
        f"{len(written)} steps written.\n"
        "Point data: damage, displacement.  Cell data: stress, strain (3x3, row-major),\n"
        "stress_vm, strain_vm.\n"
    )
    return written


def write_run_results(
    out_dir: Path | str,
    snapshots: list[dict[str, Any]],
    nodes: np.ndarray,
    cells: np.ndarray,
    *,
    boundary_nodes: np.ndarray,
    gauge_length: float,
    load_axis: int = 0,
    threshold: float = DAMAGE_THRESHOLD,
    write_vtk: bool = True,
) -> dict[str, Any]:
    """Write ``results.csv`` and the VTK series. Returns a small summary."""
    import csv

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = collect_step_rows(snapshots, nodes, cells, load_axis=load_axis,
                             gauge_length=gauge_length,
                             boundary_nodes=boundary_nodes, threshold=threshold)
    if not rows:
        raise RuntimeError(
            "no snapshots to write: the executioner needs snapshots_every > 0. "
            "Failing rather than writing an empty results file, which would "
            "look like a run that produced no damage."
        )
    cols = list(rows[0].keys())
    with (out / "results.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    files = write_vtk_series(out / "fields", snapshots, nodes, cells) if write_vtk else []
    return {"csv": str(out / "results.csv"), "n_steps": len(rows),
            "vtk_files": len(files), "columns": cols,
            "damage_threshold": threshold}
