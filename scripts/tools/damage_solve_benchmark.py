#!/usr/bin/env python3
"""Dense against sparse damage solve (PINNACLE ADR 0031) on the plate meshes.

For each Exodus mesh given, assembles the AT2 damage system with a random
history field once densely and once sparsely, applies a Dirichlet set,
solves both, and reports wall time, operator memory and the difference
between the two solutions. The dense path is skipped above ``--dense-max``
nodes (its memory grows as n^2). Usage::

    damage_solve_benchmark.py results/pfd_study/mesh/h40_hole5/mesh.e ... [--out csv]

A mesh path of the form ``hexaform:h`` generates a 100 x 50 x 1 plate with
a 5 um hole at background size h in memory (needs hexaform on the path).
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
import tracemalloc

import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
sys.path.insert(0, "/mnt/extended/projects/hexaform/src")


def load(path: str):
    from pinnacle.mesh.mesh import Mesh

    if path.startswith("hexaform:"):
        from hexaform import HoleSpec, PlateOptions, plate_mesh

        h = float(path.split(":", 1)[1])
        m = plate_mesh(PlateOptions(h=h, hole=HoleSpec("circle", 2.5)))
        mesh = Mesh(nodes=m.points, connectivity=m.cells, dim=3, node_sets={k: v for k, v in m.node_sets.items()})
        return mesh, f"hexaform h={h}"
    from pinnacle.io.mesh_readers import read_mesh

    return read_mesh("exodus", path), path


def run(mesh, label: str, dense_max: int) -> dict:
    from pinnacle.fem.elements import Hex8
    from pinnacle.phasefield.damage_field import (
        apply_damage_dirichlet,
        assemble_damage_system,
        nodal_to_quadrature,
        solve_damage_system,
    )

    element = Hex8(order=2)
    rng = np.random.default_rng(0)
    H = nodal_to_quadrature(mesh, element, rng.uniform(0.0, 2.0, mesh.n_nodes))
    mask = np.zeros(mesh.n_nodes, dtype=bool)
    xs = mesh.nodes[:, 0]
    mask[np.isclose(xs, xs.min())] = True
    mask[np.isclose(xs, xs.max())] = True
    rec = {"mesh": label, "n_nodes": mesh.n_nodes, "n_elements": mesh.n_elements}
    for sparse in (False, True):
        if not sparse and mesh.n_nodes > dense_max:
            rec.update({"dense_assemble_s": np.nan, "dense_solve_s": np.nan, "dense_mem_MB": np.nan})
            continue
        tracemalloc.start()
        t0 = time.perf_counter()
        K, r = assemble_damage_system(mesh, element, np.eye(3), Gc=1.0, l0=2.5, model="at2", H_qp=H, sparse=sparse)
        t1 = time.perf_counter()
        K, r = apply_damage_dirichlet(K, r, mask, 0.0)
        d = solve_damage_system(K, r, "direct")
        t2 = time.perf_counter()
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        tag = "sparse" if sparse else "dense"
        rec[f"{tag}_assemble_s"] = t1 - t0
        rec[f"{tag}_solve_s"] = t2 - t1
        rec[f"{tag}_mem_MB"] = peak / 1e6
        rec[f"{tag}_nnz" if sparse else "dense_entries"] = int(K.nnz) if sparse else int(K.size)
        rec[f"{tag}_d"] = d
    if "dense_d" in rec:
        rec["max_abs_diff"] = float(np.max(np.abs(rec["dense_d"] - rec["sparse_d"])))
    rec.pop("dense_d", None); rec.pop("sparse_d", None)
    return rec


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("meshes", nargs="+")
    ap.add_argument("--dense-max", type=int, default=12000)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv[1:])
    rows = []
    for p in a.meshes:
        mesh, label = load(p)
        rec = run(mesh, label, a.dense_max)
        rows.append(rec)
        print({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in rec.items()})
    if a.out and rows:
        keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in rows[0], k))
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
        print("wrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
