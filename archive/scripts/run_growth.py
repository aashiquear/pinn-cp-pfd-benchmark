#!/usr/bin/env python3
"""Run one grain-growth relaxation and capture the field at solve resolution.

WHY NOT Runner.from_config. The output layer RESAMPLES fields onto a small
grid -- ``fields['coordinates']`` came back as 14x14 for a 159x80 solve -- so
the morphology cannot be measured from a Runner result. Building the problem
directly and calling ``prob.run()`` returns ``result['field']`` at (n_nodes,
n_ops), which is the actual nodal solution.

HORIZON. Set as a fraction of the GRAIN time tau = grain_width^2/(M_GB sigma),
not the interface time. The interfaces are initialized at their equilibrium
width (int_width = w_GB), so no profile relaxation remains; everything wanted
here -- curved boundaries, 120-degree junctions -- requires boundary
MIGRATION, which is curvature-driven and lives on tau. A first attempt used
4.2 interface times, which is 6.6% of tau, and moved the free energy by 0.17%:
the structure did not relax at all.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, required=True)
    ap.add_argument("--deck", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import warnings
    warnings.filterwarnings("ignore")
    from pinnacle.plugins import discover_components
    discover_components(force=True)
    from pinnacle.backend import get_backend
    from pinnacle.io.loader import load_input
    from pinnacle.system.problem import build_composed_problem

    cfg = load_input(args.deck)
    prob = build_composed_problem(get_backend("numpy", precision="float64"), cfg)
    t0 = time.time()
    res = prob.run()
    wall = time.time() - t0

    field = np.asarray(res["field"], dtype=float)
    hist = res.get("history") or {}
    steps = hist.get("pf_fem_transient") if isinstance(hist, dict) else None
    traj = []
    if isinstance(steps, list):
        for st in steps:
            cnt = ((st.get("services") or {}).get("count") or {})
            vols = cnt.get("volumes") or []
            traj.append({
                "step": st.get("step"), "t": st.get("t"),
                "energy": st.get("energy"),
                "n_features": cnt.get("n_features"),
                "vol_min": (min(vols) if vols else None),
                "vol_max": (max(vols) if vols else None),
                "vol_cv": (float(np.std(vols) / np.mean(vols)) if vols else None),
                "newton_iters": st.get("newton_iters"),
                "converged": st.get("converged"),
            })

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "final_field.npz", field=field,
                        nodes=np.asarray(prob.mesh.nodes, dtype=float))
    e = [r["energy"] for r in traj if r["energy"] is not None]
    rec = {
        "n_grains": args.grains, "wall_s": wall, "n_steps": len(traj),
        "field_shape": list(field.shape), "n_nodes": int(prob.mesh.n_nodes),
        "energy_first": (e[0] if e else None), "energy_last": (e[-1] if e else None),
        "energy_drop_pct": (100.0 * (e[0] - e[-1]) / e[0] if e else None),
        "energy_monotone": (bool(all(b <= a * (1 + 1e-12) for a, b in zip(e, e[1:])))
                            if e else None),
        "n_features_first": (traj[0]["n_features"] if traj else None),
        "n_features_last": (traj[-1]["n_features"] if traj else None),
        "trajectory": traj,
    }
    (out / "growth_history.json").write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: v for k, v in rec.items() if k != "trajectory"}, indent=1,
                     default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
