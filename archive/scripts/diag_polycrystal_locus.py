#!/usr/bin/env python3
"""Where does the polycrystal exhaust the CP sub-step budget?

Refining the load program from 25 to 120 steps barely moved the cutback
count, and the failure happens on the CONVERGED displacement, so it is the
solution that demands the resolution rather than the increment. Two
candidate loci, which imply different remedies:

* grain boundaries -- a misorientation discontinuity across a single
  element, forced by the one-element-through-thickness pseudo-3D geometry;
* free surfaces -- with 12 grains in a 100 x 50 um plate every grain
  touches a free surface, which a bulk specimen's grains do not.

This reports the element that fails, its distance to the nearest grain
boundary and to the nearest free surface, and the misorientation across
that boundary. Speculating between the two costs nothing; measuring
settles it.
"""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
STUDY = Path(__file__).resolve().parents[1]

from pinnacle.plugins import discover_components; discover_components(force=True)
from pinnacle.io.loader import load_input
from pinnacle.runner.runner import Runner
from pinnacle.exceptions import CPSubstepBudgetError
import pinnacle.materials.cp_field_callback as CB
import run_plate as RP


def main() -> int:
    z = np.load(STUDY / "meshes/plate/plate_plain_d1.npz")
    nodes, cells = z["nodes"], z["cells"]
    grain, seeds = RP._voronoi_grains(nodes, cells, 12, seed=RP.GRAIN_SEED)
    eulers = RP._random_orientations(12, seed=RP.ORIENTATION_SEED)
    centroid = nodes[cells].mean(axis=1)

    # per-QP euler, 8 QP per hex, same ordering the harness uses
    euler_qp = np.repeat(eulers[grain], 8, axis=0)

    # Intercept the budget error to learn WHICH points failed.
    failed: list[int] = []
    orig = CB.CPFieldCallback._raise_substep_budget
    def spy(self, *a, **kw):
        # the signature carries the offending point indices where available
        failed.append(kw.get("index", a[0] if a else -1))
        return orig(self, *a, **kw)
    if hasattr(CB.CPFieldCallback, "_raise_substep_budget"):
        CB.CPFieldCallback._raise_substep_budget = spy

    out = STUDY / "results/diag/poly_locus"; out.mkdir(parents=True, exist_ok=True)
    deck = RP.build_deck(mesh=STUDY / "meshes/plate/plate_plain_d1.e",
                         orient="[1.0, 0.0, 0.0]", gc=1e9, l0=1.5, nsteps=25,
                         mechsub=9, maxalt=1, strain=0.17, out_dir=out,
                         cal=json.load(open(STUDY / "results/calibration/fe/cp_taylor.json"))
                             and RP.CAL_AL2219,
                         maxcutbacks=3, t_end=170.0)
    cfg = load_input(str(deck)); cfg["orientation_field"] = {"euler_angles": euler_qp}
    try:
        Runner.from_config(cfg)
        print("ran to completion -- no budget exhaustion to locate")
        return 0
    except CPSubstepBudgetError as exc:
        print("budget exhausted:", str(exc)[:200])
    except Exception as exc:
        print(f"{type(exc).__name__}: {str(exc)[:200]}")

    # Geometry of the grain structure, independent of which point failed.
    L, W = 100.0, 50.0
    on_surface = (
        (np.abs(centroid[:, 0]) < 5) | (np.abs(centroid[:, 0] - L) < 5)
        | (np.abs(centroid[:, 1]) < 5) | (np.abs(centroid[:, 1] - W) < 5)
    )
    # an element is "at a boundary" if any face-neighbour has a different grain
    from collections import defaultdict
    face_owner = defaultdict(list)
    for e, c in enumerate(cells):
        for f in ((0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)):
            face_owner[tuple(sorted(c[list(f)]))].append(e)
    at_gb = np.zeros(len(cells), bool)
    for owners in face_owner.values():
        if len(owners) == 2 and grain[owners[0]] != grain[owners[1]]:
            at_gb[owners[0]] = at_gb[owners[1]] = True
    rec = {
        "n_elements": int(len(cells)), "n_grains": 12,
        "elements_at_grain_boundary": int(at_gb.sum()),
        "fraction_at_grain_boundary": float(at_gb.mean()),
        "elements_touching_free_surface": int(on_surface.sum()),
        "fraction_touching_free_surface": float(on_surface.mean()),
        "elements_both": int((at_gb & on_surface).sum()),
        "elements_neither": int((~at_gb & ~on_surface).sum()),
        "grain_sizes": np.bincount(grain, minlength=12).tolist(),
    }
    (out / "locus.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
