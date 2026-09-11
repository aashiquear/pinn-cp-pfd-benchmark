#!/usr/bin/env python3
"""Wall-time assessment for the pseudo-3D plate, one density level per call."""
import json, re, sys, time, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
from pinnacle.plugins import discover_components; discover_components(force=True)
from pinnacle.io.loader import load_input
from pinnacle.runner.runner import Runner

STUDY = Path(__file__).resolve().parent.parent
CAL = {"TAUC0": "60.31771394231891", "RHO0": "2.393387495846754",
       "K1": "0.08114503403820844", "YC2": "0.0020581167904802733"}

def build(d, *, gc, l0, nsteps, mechsub, maxalt, strain, tag, orient="[1.0, 0.0, 0.0]"):
    src = (STUDY / "decks/al2219/plate_template.i").read_text()
    sub = dict(CAL, MESH=str(STUDY / f"meshes/plate/plate_hole_d{d}.e"),
               ORIENT=orient, GC=f"{gc:g}", L0=f"{l0:g}", NSTEPS=str(nsteps),
               MECHSUB=str(mechsub), MAXALT=str(maxalt), STRAIN=f"{strain:g}",
               OUTDIR=str(STUDY / f"runs/plate_{tag}"))
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    assert "@" not in re.sub(r"#.*", "", src), "unsubstituted placeholder"
    p = STUDY / f"runs/plate_{tag}.i"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(src)
    return p

if __name__ == "__main__":
    d = int(sys.argv[1])
    fam = {r["density"]: r for r in json.load(open(STUDY / "meshes/plate/plate_family.json"))}
    ne = fam[d]["n_elements"]
    deck = build(d, gc=10.0, l0=1.0, nsteps=1, mechsub=1, maxalt=1, strain=0.01, tag=f"t{d}")
    t0 = time.perf_counter(); Runner.from_config(load_input(str(deck))); dt = time.perf_counter() - t0
    print(json.dumps({"density": d, "n_elements": ne, "n_dofs": fam[d]["n_dofs_mech"],
                      "seconds_one_step": dt, "ms_per_element": 1e3 * dt / ne}))

def timed(d, *, nsteps, mechsub, maxalt, strain, gc=10.0, l0=1.0):
    """Run and report wall time WITH the Newton iteration count.

    A single-load-step timing is not a size measurement: iteration count
    varies with how hard the step is, and at density 4 it made a smaller mesh
    three times slower than a larger one. Cost per element is only meaningful
    divided by the work actually done.
    """
    import pinnacle.solvers.linear as L
    from pinnacle.fem import assembly as A
    n = {"tangent": 0}
    orig = A.Assembler.tangent
    def spy(self, *a, **k):
        n["tangent"] += 1
        return orig(self, *a, **k)
    A.Assembler.tangent = spy
    deck = build(d, gc=gc, l0=l0, nsteps=nsteps, mechsub=mechsub, maxalt=maxalt,
                 strain=strain, tag=f"m{d}")
    t0 = time.perf_counter(); Runner.from_config(load_input(str(deck))); dt = time.perf_counter() - t0
    A.Assembler.tangent = orig
    return dt, n["tangent"]
