#!/usr/bin/env python3
"""The Fig. 3 polycrystal calibration specimen: 30 grains, 50 x 50 x 100 um.

TWO JOBS, selected by --gc.

  --gc 0   no damage. Runs the calibrated CP set through a finite-element
           polycrystal and compares against the Fig. 3 experiment. This
           measures the Taylor-to-FE offset, which the CP calibration
           carries as an explicit caveat: Taylor is an isostrain upper
           bound (compatible, not equilibrated), so the finite-element
           aggregate must come out softer, and by how much is a
           measurement rather than a guess.

  --gc >0  damage on. Sweeps Gc against the experimental fracture onset:
           the test was stopped at 17% strain by catastrophic fracture,
           so 17% is the initiation target.

INITIATION CRITERION. Primary is c >= 0.85 on the ELEMENT NODAL AVERAGE,
the common phase-field convention. The paper's own criterion -- 0.65 at any
node -- is recorded alongside it, because the published initiation strains
this study compares against were produced with that one, and the two are not
interchangeable: a nodal peak crosses any threshold earlier than an element
average of the same field.

WHY THE MESH IS COARSER THAN THE VOID STUDY'S, AND WHY THAT IS NOT A
SHORTCUT. Resolving l0 = 1 um at h <= l0/2 across this domain needs 2e6
elements, about 795 h per run at the measured cost -- not affordable, and no
amount of scheduling changes that. The specimen is therefore run at a larger
l0 with h/l0 held at the same 0.5, so the resolution rule is satisfied and
only the regularization length differs. That matters for transferring Gc: in
the standard phase-field model the critical stress scales as
sqrt(Gc*E/l0), so Gc calibrated at one l0 does NOT transfer to another at
fixed Gc -- it transfers at fixed Gc/l0. That prediction is testable and is
tested, rather than assumed, against the single-crystal calibration.

GRAIN GEOMETRY. Voronoi tessellation over a structured hex mesh, so grain
boundaries are stepped rather than conforming. Recorded as a limitation. At
30 grains in this domain the equivalent grain edge is about 20 um, so at
h = 2.5 um a grain spans roughly 8 elements.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent
PINNACLE = Path("/mnt/extended/projects/pinnacle")

DOMAIN = (50.0, 50.0, 100.0)          # um, Sec. 2.4
N_GRAINS = 30                          # Sec. 2.4
ORIENTATION_SEED = 20260910            # same realization as the CP calibration
GRAIN_SEED = 20260911

CALIBRATED_CP = {
    "tau_c0": 60.317714, "rho0": 2.393387,
    "k1": 0.081145, "annihilation_diameter": 0.002058,
}
FIXED_CP = {
    "C11": 108.0e3, "C12": 61.3e3, "C44": 28.5e3, "mu": 25.7e3,
    "gamma0": 1.0e-5, "rate_exponent_n": 1.0 / 0.12,
    "q_lat": 1.4, "alpha0": 0.4, "burgers": 0.000286,
}
EXPERIMENTAL_FRACTURE_STRAIN = 0.17    # test stopped here by catastrophic fracture


def random_orientations(n: int, seed: int) -> np.ndarray:
    """Bunge Euler triples (degrees), uniform on SO(3).

    cos(Phi) uniform, not Phi uniform: the naive draw clusters near the poles
    and biases the aggregate stiffness.
    """
    rng = np.random.default_rng(seed)
    phi1 = rng.uniform(0.0, 2.0 * np.pi, n)
    Phi = np.arccos(rng.uniform(-1.0, 1.0, n))
    phi2 = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.degrees(np.column_stack([phi1, Phi, phi2]))


def voronoi_grains(centroids: np.ndarray, n_grains: int, seed: int) -> np.ndarray:
    """Nearest-seed grain id per element."""
    rng = np.random.default_rng(seed)
    seeds = np.column_stack([rng.uniform(0.0, L, n_grains) for L in DOMAIN])
    d2 = ((centroids[:, None, :] - seeds[None, :, :]) ** 2).sum(axis=2)
    return np.argmin(d2, axis=1), seeds


def element_average(field_nodal: np.ndarray, conn: np.ndarray) -> np.ndarray:
    return field_nodal[conn].mean(axis=1)


def main() -> int:
    global STUDY, PINNACLE
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=2.5, help="element size, um")
    ap.add_argument("--l0", type=float, default=None,
                    help="phase-field length; defaults to 2*h (the h <= l0/2 rule)")
    ap.add_argument("--gc", type=float, default=0.0, help="0 disables damage")
    ap.add_argument("--strain", type=float, default=0.20)
    ap.add_argument("--load-steps", type=int, default=20)
    ap.add_argument("--max-alternations", type=int, default=3)
    ap.add_argument("--out", default=None)
    ap.add_argument("--study-root", default=None)
    ap.add_argument("--pinnacle-root", default=None)
    ap.add_argument("--pinnacle-rev", default=None)
    args = ap.parse_args()

    if args.study_root:
        STUDY = Path(args.study_root).resolve()
    if args.pinnacle_root:
        PINNACLE = Path(args.pinnacle_root).resolve()
    out_dir = Path(args.out or (STUDY / "results" / "polycrystal"))
    out_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(PINNACLE / "src"))
    sys.path.insert(0, str(Path(__file__).resolve().parent))

    from pinnacle import registry
    from pinnacle.plugins import discover_components
    from pinnacle.runner.runner import Runner

    discover_components(force=True)

    l0 = args.l0 if args.l0 is not None else 2.0 * args.h
    damage_on = args.gc > 0.0
    div = [max(1, int(round(L / args.h))) for L in DOMAIN]
    # NOTE the two spellings: the DECK key is ``domain`` (flat min/max corner
    # pair), while the factory keyword is ``bounds`` (per-axis pairs). Passing
    # the deck spelling here raises rather than being silently ignored.
    mesh = registry.create("mesh", "box_mesh",
                           bounds=[[0.0, L] for L in DOMAIN], divisions=div)
    nodes = np.asarray(mesh.nodes, dtype=float).reshape(-1, 3)
    conn = np.asarray(mesh.connectivity, dtype=int)
    centroids = nodes[conn].mean(axis=1)
    grain_id, seeds = voronoi_grains(centroids, N_GRAINS, GRAIN_SEED)
    eulers = random_orientations(N_GRAINS, ORIENTATION_SEED)

    # One orientation per element, repeated over the element's quadrature
    # points. Assigning per ELEMENT rather than per QP means the mapping does
    # not depend on the local quadrature ordering, which is not part of the
    # public interface.
    n_qp_per_elem = 8
    euler_qp = np.repeat(eulers[grain_id], n_qp_per_elem, axis=0)

    tag = f"h{args.h:g}_l0{l0:g}_Gc{args.gc:g}"
    deck_src = (STUDY / "decks" / "al2219" / "polycrystal_template.i").read_text()
    axial_um = args.strain * DOMAIN[2]
    subs = {
        "@OUTDIR@": str(out_dir / f"run_{tag}"),
        "@DIVISIONS@": "[" + ", ".join(str(d) for d in div) + "]",
        "@AXIAL_EXPR@": f"{axial_um:.9g}*t/10.0",
        "@TAU_C0@": f"{CALIBRATED_CP['tau_c0']:.9g}",
        "@RHO0@": f"{CALIBRATED_CP['rho0']:.9g}",
        "@K1@": f"{CALIBRATED_CP['k1']:.9g}",
        "@YC2@": f"{CALIBRATED_CP['annihilation_diameter']:.9g}",
        # --gc 0 means "no damage", and the way to get that is a Gc so large
        # the driving force never reaches it -- NOT a Gc near zero, which is
        # the opposite: zero fracture toughness fails every point immediately
        # and degrades the stress the run is supposed to measure. A first
        # version used 1e-12 and reported peak damage 1.0 at the first step.
        "@GC@": f"{(1.0e9 if not damage_on else args.gc):.9g}",
        "@L0@": f"{l0:.9g}",
        "@LOAD_STEPS@": str(args.load_steps),
        "@MAX_ALT@": str(args.max_alternations if damage_on else 1),
    }
    for key, val in subs.items():
        if key not in deck_src:
            raise SystemExit(f"placeholder {key} missing from the template")
        deck_src = deck_src.replace(key, val)
    deck = out_dir / f"polycrystal_{tag}.i"
    deck.write_text(deck_src)

    record = {
        "h_um": args.h, "l0_um": l0, "h_over_l0": args.h / l0,
        "Gc_MPa_um": args.gc, "damage": damage_on,
        "divisions": div, "n_elements": int(conn.shape[0]),
        "n_nodes": int(nodes.shape[0]), "n_grains": N_GRAINS,
        "grain_seed": GRAIN_SEED, "orientation_seed": ORIENTATION_SEED,
        "grain_sizes": np.bincount(grain_id, minlength=N_GRAINS).tolist(),
        "cp_parameters": {**FIXED_CP, **CALIBRATED_CP},
        "experimental_fracture_strain": EXPERIMENTAL_FRACTURE_STRAIN,
        "applied_strain": args.strain, "load_steps": args.load_steps,
        "grain_boundaries": "stepped (Voronoi over a structured hex mesh), not conforming",
    }
    # ---- run ----
    from pinnacle.io.loader import load_input

    cfg = load_input(str(deck))
    # The orientation field is injected here rather than written into the
    # deck: it is one Euler triple per quadrature point (~128k numbers at
    # h = 2.5 um). ``problem.py`` reads it from this exact key.
    cfg["orientation_field"] = {"euler_angles": euler_qp}

    t0 = time.time()
    result = Runner.from_config(cfg)
    record["wall_seconds"] = time.time() - t0

    def _find(obj, key, depth=0):
        if depth > 7:
            return None
        if isinstance(obj, dict):
            if key in obj:
                return obj[key]
            for v in obj.values():
                got = _find(v, key, depth + 1)
                if got is not None:
                    return got
        return None

    snapshots = _find(result, "snapshots") or []
    history = _find(result, "staggered") or []
    curve, init_085, init_065 = [], None, None
    for snap in snapshots:
        eps = float(snap["step_factor"]) * args.strain
        row = {"step": int(snap["step"]), "strain": eps}
        sig = snap.get("stress")
        if sig is not None:
            row["stress_entries"] = int(np.asarray(sig).size)
            sig = np.asarray(sig, dtype=float)
            # Volume-averaged axial Cauchy stress. Equal element volumes on a
            # structured mesh, so a plain mean is the volume average.
            if sig.size % 9 == 0:
                t = sig.reshape(-1, 3, 3)
                row["sigma_zz_MPa"] = float(t[:, 2, 2].mean())
                # The transverse averages are the check that this is a
                # uniaxial-STRESS test. The lateral faces are free, so both
                # should be small against sigma_zz; if they are not, the
                # roller faces are constraining lateral contraction and the
                # run is closer to uniaxial strain, which would over-predict
                # the axial stress for reasons that have nothing to do with
                # the material.
                row["sigma_xx_MPa"] = float(t[:, 0, 0].mean())
                row["sigma_yy_MPa"] = float(t[:, 1, 1].mean())
                row["transverse_over_axial"] = float(
                    max(abs(t[:, 0, 0].mean()), abs(t[:, 1, 1].mean()))
                    / max(1.0e-30, abs(t[:, 2, 2].mean())))
            else:
                row["sigma_zz_MPa"] = float(np.mean(sig))
        d = snap.get("damage")
        if d is not None:
            d = np.asarray(d, dtype=float).ravel()
            d_elem = element_average(d, conn)
            row["peak_damage_nodal"] = float(d.max())
            row["peak_damage_element_avg"] = float(d_elem.max())
            if init_085 is None and d_elem.max() >= 0.85:
                init_085 = eps
            if init_065 is None and d.max() >= 0.65:
                init_065 = eps
        curve.append(row)

    record["damage_suppressed_by_large_Gc"] = (not damage_on)
    record["curve"] = curve
    record["initiation_strain_element_avg_0p85"] = init_085
    record["initiation_strain_nodal_0p65"] = init_065
    record["initiation_resolution"] = args.strain / max(1, args.load_steps)
    record["all_steps_converged"] = (
        all(h.get("converged") for h in history) if history else None)
    record["alternations_per_step"] = [int(h.get("iterations", 0)) for h in history]

    if damage_on and snapshots:
        from damage_evidence import write_evidence

        try:
            record["evidence"] = write_evidence(
                snapshots, nodes, out_dir / f"evidence_{tag}",
                label=f"Al2219 polycrystal h={args.h:g} l0={l0:g} Gc={args.gc:g}",
                slice_axis=0)
        except Exception as exc:
            record["evidence"] = {"error": f"{type(exc).__name__}: {exc}"}

    rev = args.pinnacle_rev
    if rev is None:
        try:
            rev = subprocess.run(["git", "-C", str(PINNACLE), "rev-parse", "--short", "HEAD"],
                                 capture_output=True, text=True, timeout=20).stdout.strip()
        except Exception:
            rev = "unknown"
    record["provenance"] = {"pinnacle_rev": rev, "deck": deck.name,
                            "backend": "numpy", "precision": "float64"}

    (out_dir / f"polycrystal_{tag}.json").write_text(
        json.dumps(record, indent=1, default=str))
    print(json.dumps({k: record[k] for k in
                      ("h_um", "l0_um", "h_over_l0", "n_elements", "n_nodes",
                       "divisions", "Gc_MPa_um", "damage")}, indent=1, default=str))
    print(f"grain element counts: min {min(record['grain_sizes'])}, "
          f"max {max(record['grain_sizes'])}, mean "
          f"{np.mean(record['grain_sizes']):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
