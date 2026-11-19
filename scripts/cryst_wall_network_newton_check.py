#!/usr/bin/env python3
"""Case B1: plain Newton started from saved fields (which stationary state does it reach?).

Usage: python3 scripts/cryst_wall_network_newton_check.py
For every pair (plain-Newton deck, run whose saved field starts the iteration) below, the
deck is run with `pf_fem_transient.set_initial_state` and the result is compared with the
field of the FE saddle search. Also prints the discrete L2 norm of the assembled residual
(`sum R_i^2/m_i / sum m_i`, free unknowns) of the starting field.
Writes results/crystallization/walls/newton_check.csv.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vtk_legacy import read_vtk  # noqa: E402

from pinnacle.backend import get_backend  # noqa: E402
from pinnacle.io import load_input  # noqa: E402
from pinnacle.plugins import discover_components  # noqa: E402
from pinnacle.system.problem import build_composed_problem  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/crystallization/walls"
DECKS = ROOT / "decks/crystallization/walls"

M2 = "b1_wall_model2_dT{:03d}_nucleus"
PAIRS = (
    [(M2.format(40), f"{M2.format(40)}_pinn_seed{s}") for s in (0, 1, 2)]
    + [(M2.format(d), f"{M2.format(d)}_fem") for d in (20, 40, 60, 80)]
    + [(M2.format(80), f"{M2.format(80)}_hybrid_seed0"), (M2.format(40), f"{M2.format(40)}_hybrid_seed2"),
       (M2.format(20), f"{M2.format(20)}_hybrid150_seed0")]
)


def field(run: str) -> np.ndarray:
    return np.asarray(read_vtk(RES / run / "fields.vtk")["point"]["phi"], dtype=float)


def main() -> None:
    discover_components(force=True)
    rows = []
    for stem, start_run in PAIRS:
        if not (RES / start_run / "fields.vtk").exists():
            continue
        u_fe, u0 = field(f"{stem}_fem"), field(start_run)
        config = load_input(str(DECKS / f"{stem}_newton.i"))
        config.setdefault("outputs", {})
        backend = get_backend("numpy", precision="float64")
        problem = build_composed_problem(backend, config)
        ev = problem.executioner.residual_evaluator(problem)
        r, _ = ev.evaluate(u0)
        free = ~ev.fixed
        norm = float(np.sqrt(np.sum(r[free] ** 2 / ev.weights[free]) / np.sum(ev.weights[free])))
        row = {"deck": f"{stem}_newton", "start": start_run, "start_residual_l2": norm,
               "start_minus_saddle_max": float(np.abs(u0 - u_fe).max())}
        problem = build_composed_problem(backend, config)
        problem.executioner.set_initial_state({"phi": u0})
        try:
            out = problem.run()
            u = np.asarray(out["field_phi"], dtype=float).reshape(-1)
            hist = out["history"]["pf_fem_transient"]
            row.update(converged=1, newton=int(hist[0]["newton_iters"]),
                       end_minus_saddle_max=float(np.abs(u - u_fe).max()),
                       end_minus_start_max=float(np.abs(u - u0).max()))
        except Exception as err:  # noqa: BLE001  (the deck's refusal is the result)
            row.update(converged=0, newton="", end_minus_saddle_max="", end_minus_start_max="")
            print("  ", str(err)[:120])
        rows.append(row)
        print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()), flush=True)
    with (RES / "newton_check.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
