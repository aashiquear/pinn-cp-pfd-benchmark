#!/usr/bin/env python3
"""Why does the load step at crack propagation refuse to solve?

Three remedies have already failed to move it: more alternations (+1 step),
a 3x finer load-increment floor (no change at all), and the degradation fix
(got further, still fails). This records what the solver actually sees at
each alternation so the cause is measured rather than inferred.

The hypothesis under test: the strain increment the CP integrator must
swallow is driven by the DAMAGE jump between alternations, not by the load
increment. If so, shrinking the load step cannot help -- which is exactly
what the cutback experiment showed -- and the fix belongs in how the
staggered scheme advances damage.
"""
import json, sys, warnings
from pathlib import Path
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pinnacle.plugins import discover_components  # noqa: E402

discover_components(force=True)
from pinnacle.io.loader import load_input  # noqa: E402
from pinnacle.runner.runner import Runner  # noqa: E402
from pinnacle.materials.cp_field_callback import CPFieldCallback  # noqa: E402
import run_plate as RP  # noqa: E402

STRAIN, NSTEPS, ALTS, GC, L0 = 0.12, 60, 2, 30.0, 1.5
out = RP.STUDY / "results/plate/diag_snapback"
out.mkdir(parents=True, exist_ok=True)
z = np.load(RP.STUDY / "meshes/plate/plate_hole_d1.npz")
mechsub = max(1, int(np.ceil((STRAIN / NSTEPS) / 8.0e-4)))
deck = RP.build_deck(mesh=RP.STUDY / "meshes/plate/plate_hole_d1.e",
                     orient="[1.0, 0.0, 0.0]", gc=GC, l0=L0, nsteps=NSTEPS,
                     mechsub=mechsub, maxalt=ALTS, strain=STRAIN,
                     out_dir=out, cal=RP.CAL_AL2219, maxcutbacks=6)

log = []
prev = {"d": None, "eps": None}

_orig_set = CPFieldCallback.set_damage_field
def set_damage_field(self, d):
    d = np.asarray(d, dtype=float)
    pv = prev["d"]
    log.append({"event": "damage_push", "max_d": float(d.max()),
                "max_jump": None if pv is None else float(np.abs(d - pv).max()),
                "n_above_085": int((d > 0.85).sum())})
    prev["d"] = d.copy()
    return _orig_set(self, d)

_orig_update = CPFieldCallback.update_state
def update_state(self, strain, dt):
    e = np.asarray(self.backend.to_numpy(strain), dtype=float)
    pe = prev["eps"]
    rec = {"event": "update_state", "dt": float(dt),
           "max_abs_eps": float(np.abs(e).max()),
           "max_d_eps": None if pe is None or pe.shape != e.shape
                        else float(np.abs(e - pe).max())}
    prev["eps"] = e.copy()
    try:
        r = _orig_update(self, strain, dt)
        rec["ok"] = True
    except Exception as exc:                      # noqa: BLE001
        rec["ok"] = False
        rec["error"] = f"{type(exc).__name__}: {exc}"
        log.append(rec)
        raise
    log.append(rec)
    return r

CPFieldCallback.set_damage_field = set_damage_field
CPFieldCallback.update_state = update_state

try:
    Runner.from_config(load_input(str(deck)))
    status = "completed"
except Exception as exc:                          # noqa: BLE001
    status = f"{type(exc).__name__}: {exc}"

(out / "trace.json").write_text(json.dumps({"status": status, "log": log}, indent=1))

pushes = [x for x in log if x["event"] == "damage_push"]
fails = [x for x in log if x["event"] == "update_state" and not x.get("ok", True)]
print(f"status: {status[:180]}")
print(f"damage pushes: {len(pushes)} | update_state calls: {len(log)-len(pushes)} | failures: {len(fails)}")
print(f"{'push#':>6} {'max_d':>8} {'max_jump':>10} {'n>0.85':>7}")
for i, p in enumerate(pushes):
    if p["max_jump"] is None or p["max_jump"] > 0.05 or i >= len(pushes) - 4:
        print(f"{i:6d} {p['max_d']:8.4f} "
              f"{-1 if p['max_jump'] is None else p['max_jump']:10.4f} {p['n_above_085']:7d}")
if fails:
    print("FIRST FAILURE:", json.dumps(fails[0], indent=1)[:700])
