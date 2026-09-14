#!/usr/bin/env python3
"""Run one pseudo-3D plate case and write the standard outputs.

Every plate run in the study goes through here, so that no result can exist
without its ``results.csv`` and full-field VTK series. The G20 constitutive
strain-increment limit is derived from the loading rather than passed in, so
a caller cannot under-resolve a run by forgetting a flag.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import warnings
from pathlib import Path

import sys

import numpy as np

warnings.filterwarnings("ignore")
# The harness is run by path from the study root, so its own directory is not
# necessarily on sys.path; the shared results writer lives beside it.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pinnacle.plugins import discover_components  # noqa: E402

discover_components(force=True)
from pinnacle.io.loader import load_input  # noqa: E402
from pinnacle.runner.runner import Runner  # noqa: E402

from run_results import write_run_results  # noqa: E402

STUDY = Path(__file__).resolve().parent.parent
MAX_STRAIN_INCREMENT = 8.0e-4      # per constitutive update; study gap G20
PLATE_LENGTH = 100.0               # um, the loaded (x) dimension

# Calibrated Al2219 CP set (plain-plate refit supersedes this in Phase 1).
CAL_AL2219 = {
    "NEXP": "8.333333333333334",
    "TAUC0": "60.31771394231891", "RHO0": "2.393387495846754",
    "K1": "0.08114503403820844", "YC2": "0.0020581167904802733",
}


def _random_orientations(n: int, seed: int) -> np.ndarray:
    """Uniformly random orientations as Bunge Euler triples (degrees).

    Uniform in SO(3), not uniform in the Euler angles: phi is drawn from
    arccos of a uniform variate, because sampling all three angles uniformly
    clusters orientations at the poles.
    """
    rng = np.random.default_rng(seed)
    phi1 = rng.uniform(0.0, 360.0, n)
    Phi = np.degrees(np.arccos(rng.uniform(-1.0, 1.0, n)))
    phi2 = rng.uniform(0.0, 360.0, n)
    return np.column_stack([phi1, Phi, phi2])


def _voronoi_grains(pts2: np.ndarray, n_grains: int, seed: int):
    """Assign each point to the nearest of n_grains random seed points."""
    rng = np.random.default_rng(seed)
    lo, hi = pts2.min(axis=0), pts2.max(axis=0)
    seeds = rng.uniform(lo, hi, size=(n_grains, pts2.shape[1]))
    d2 = ((pts2[:, None, :] - seeds[None, :, :]) ** 2).sum(axis=2)
    return np.argmin(d2, axis=1), seeds


def _equiaxed_grains(pts2: np.ndarray, target_size: float, extent, seed: int):
    """Voronoi grains of near-equal size, seeded on a jittered triangular lattice.

    Uniformly random seeds (Poisson) give a broad grain-size distribution:
    some grains come out a few elements across while others span a quarter of
    the specimen, which is what made the first 12-grain specimen pathological
    (smallest grain 3 elements, 69% of elements on a misorientation
    discontinuity). A jittered lattice fixes the mean spacing, so every grain
    has area ~ target_size**2 and the element count per grain is roughly
    uniform.

    ``target_size`` should be set to the plate THICKNESS for a genuinely
    equiaxed structure: the pseudo-3D mesh is one element through thickness,
    so a grain is a prism of height 10 um, and only an in-plane size of about
    10 um makes it equiaxed rather than columnar.
    """
    L, W = float(extent[0]), float(extent[1])
    d = float(target_size)
    row_h = d * np.sqrt(3.0) / 2.0          # triangular lattice row spacing
    rng = np.random.default_rng(seed)
    seeds = []
    n_rows = max(1, int(np.ceil(W / row_h)))
    for j in range(n_rows):
        y = (j + 0.5) * row_h
        offset = 0.5 * d if (j % 2) else 0.0
        n_cols = max(1, int(np.ceil(L / d)))
        for i in range(n_cols):
            x = offset + (i + 0.5) * d
            if x > L + 0.5 * d:
                continue
            # jitter by a quarter spacing: keeps the structure irregular
            # (grain boundaries are not all aligned with the mesh) without
            # reintroducing the near-degenerate grains a Poisson process gives
            seeds.append((x + rng.uniform(-0.25, 0.25) * d,
                          y + rng.uniform(-0.25, 0.25) * row_h))
    seeds = np.asarray(seeds, dtype=float)
    d2 = ((pts2[:, None, :] - seeds[None, :, :]) ** 2).sum(axis=2)
    gid = np.argmin(d2, axis=1)
    # drop empty seeds and relabel consecutively
    used, gid = np.unique(gid, return_inverse=True)
    return gid, seeds[used]


def _nvl(value):
    """Absent status is reported as such, never as a falsy-looking None."""
    return "not_reported" if value is None else value


BC_MINIMAL = """[BCs]
  # Statically minimal: exactly enough to remove the six rigid-body modes.
  # Both end faces stay free to contract laterally and through-thickness, so
  # the only stress concentrator is the hole and no grip state is imposed.
  # Measured against the clamped set on an identical mesh, this holds
  # |sigma_zz|/sigma_xx_max to <= 6% everywhere against 46% for clamped ends.
  [./xmin_fixed]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./pin_y]
    variable   = u
    kind       = dirichlet
    node_set   = pin_y
    components = [1]
    values     = [0.0]
  []
  [./zmin_roller]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
  [./xmax_tension]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
    function   = axial
  []
[]
"""

BC_CLAMPED = """[BCs]
  # Clamped ends: xmin fully fixed, xmax prescribed in x with u_y = u_z = 0.
  # Physically a gripped specimen, but it suppresses Poisson contraction at
  # both ends and so adds a triaxial boundary layer; measured decay is 22% of
  # peak sigma_zz at the end faces falling to 0.76% at the hole band.
  [./xmin_fixed]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0, 1, 2]
    values     = [0.0, 0.0, 0.0]
  []
  [./xmax_tension]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0, 1, 2]
    values     = [1.0, 0.0, 0.0]
    function   = axial
  []
[]
"""

BC_FACEROLLER = """[BCs]
  # Face rollers with one fixed end: xmin holds u_x = 0, ymin holds u_y = 0,
  # zmin holds u_z = 0, xmax carries the prescribed u_x.
  #
  # Differs from BC_MINIMAL only in the transverse restraint -- a whole FACE
  # of y-rollers instead of a single pin line. It removes the same rigid-body
  # rotation but also forbids the ymin face from narrowing, so the specimen
  # cannot contract freely on that side. Whether that is benign is
  # measurable, not assumed: a thin plate in uniaxial tension must carry
  # sigma_yy ~ 0, so transverse stress near ymin is the test.
  #
  # Differs from BC_CLAMPED in NOT holding u_y or u_z on the LOADED face, so
  # no triaxial grip state is created where the load is applied.
  [./xmin_fixed]
    variable   = u
    kind       = dirichlet
    node_set   = xmin
    components = [0]
    values     = [0.0]
  []
  [./ymin_roller]
    variable   = u
    kind       = dirichlet
    node_set   = ymin
    components = [1]
    values     = [0.0]
  []
  [./zmin_roller]
    variable   = u
    kind       = dirichlet
    node_set   = zmin
    components = [2]
    values     = [0.0]
  []
  [./xmax_tension]
    variable   = u
    kind       = dirichlet
    node_set   = xmax
    components = [0]
    values     = [1.0]
    function   = axial
  []
[]
"""

BC_SETS = {"minimal": BC_MINIMAL, "clamped": BC_CLAMPED,
           "faceroller": BC_FACEROLLER}

def build_deck(*, mesh: Path, orient: str, gc: float, l0: float, nsteps: int,
               mechsub: int, maxalt: int, strain: float, out_dir: Path,
               cal: dict[str, str], maxcutbacks: int = 15,
               t_end: float = 10.0, bc: str = "minimal") -> Path:
    src = (STUDY / "decks/al2219/plate_template.i").read_text()
    sub = dict(cal, MESH=str(mesh), ORIENT=orient, GC=f"{gc:.9g}", L0=f"{l0:.9g}",
               NSTEPS=str(nsteps), MECHSUB=str(mechsub), MAXALT=str(maxalt),
               MAXCUTBACKS=str(maxcutbacks), TEND=f'{t_end:.10g}',
               STRAIN=f"{strain:.9g}", OUTDIR=str(out_dir / "solver"),
               BCBLOCK=BC_SETS[bc])
    for k, v in sub.items():
        src = src.replace(f"@{k}@", v)
    leftover = set(re.findall(r"@([A-Z0-9_]+)@", re.sub(r"#.*", "", src)))
    if leftover:
        raise SystemExit(f"unsubstituted placeholders: {sorted(leftover)}")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "deck.i"
    path.write_text(src)
    return path


def _find_list(obj, key, depth=0):
    """Recursive search for ``key`` whose value is a LIST.

    Distinguishes the wanted leaf from same-named container dicts on the way
    down, which a plain key search cannot do.
    """
    if depth > 7:
        return None
    if isinstance(obj, dict):
        got = obj.get(key)
        if isinstance(got, list):
            return got
        for value in obj.values():
            found = _find_list(value, key, depth + 1)
            if found is not None:
                return found
    return None


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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--density", type=int, required=True)
    ap.add_argument("--geometry", choices=("hole", "plain"), default="hole")
    ap.add_argument("--orientation", default="[1.0, 0.0, 0.0]")
    ap.add_argument("--nexp", type=float, default=0.0,
                    help="Rate exponent n = 1/m_s. 0 keeps the set's value.")
    ap.add_argument("--equiaxed", type=float, default=0.0,
                    help="Target in-plane grain size (um). Overrides --grains "
                         "and builds equiaxed grains on a jittered triangular "
                         "lattice. Set this to the plate thickness (10 um) for "
                         "a genuinely equiaxed structure.")
    ap.add_argument("--grains", type=int, default=0,
                    help="0 = single crystal. N > 0 builds an N-grain columnar "
                         "polycrystal by Voronoi tessellation of the element "
                         "centroids, with one random orientation per grain.")
    ap.add_argument("--grain-seed", type=int, default=20260911)
    ap.add_argument("--orientation-seed", type=int, default=20260910)
    ap.add_argument("--cp-set", default=None,
                    help="JSON file with TAUC0/RHO0/K1/YC2 overriding the built-in "
                         "calibrated set (used by the CP calibration loop).")
    ap.add_argument("--gc", type=float, default=10.0)
    ap.add_argument("--l0", type=float, default=1.0)
    ap.add_argument("--strain", type=float, default=0.12)
    ap.add_argument("--strain-rate", type=float, default=1.0e-3,
                    help="Imposed engineering strain rate (1/s). The constitutive "
                         "law is rate dependent (power-law slip, reference rate "
                         "gamma0 = 1e-5 /s), so this MUST match the rate the "
                         "parameters were calibrated at or the set does not "
                         "transfer. The deck's t_end is derived as strain/rate.")
    # Denser by default than the smoke runs: at 20 steps to 12% strain the
    # first recorded point is already at 0.6%, past yield, so the elastic
    # branch is one point wide and fracture initiation is bracketed to a
    # whole load step. 60 steps puts the first point at 0.2%. Load steps
    # multiply damage solves but NOT constitutive work, since the sub-step
    # count per step falls as the increment shrinks.
    ap.add_argument("--load-steps", type=int, default=60)
    ap.add_argument("--alternations", type=int, default=2)
    ap.add_argument("--max-cutbacks", type=int, default=15,
                    help="Load-increment halvings before a step is abandoned. "
                         "The framework default of 6 bottoms out at 1/64 of a "
                         "load step, which is too coarse at the crack tip.")
    ap.add_argument("--bc", choices=tuple(BC_SETS), default="minimal",
                    help="Boundary-condition set; see the BC_* blocks.")
    ap.add_argument("--no-vtk", action="store_true")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", default="results/plate")
    args = ap.parse_args()

    # G20: the constitutive increment per update must stay at or below 8e-4.
    # Derived, not passed, so a run cannot be silently under-resolved.
    d_eps = args.strain / max(1, args.load_steps)
    mechsub = max(1, int(np.ceil(d_eps / MAX_STRAIN_INCREMENT)))

    tag = args.tag or f"{args.geometry}_d{args.density}_Gc{args.gc:g}_l0{args.l0:g}"
    out_dir = STUDY / args.out / tag
    mesh_npz = STUDY / f"meshes/plate/plate_{args.geometry}_d{args.density}.npz"
    mesh_e = STUDY / f"meshes/plate/plate_{args.geometry}_d{args.density}.e"
    z = np.load(mesh_npz)
    nodes, cells = z["nodes"], z["cells"]

    cal = dict(CAL_AL2219)
    if args.cp_set:
        cal.update({k: f"{float(v):.17g}" for k, v in json.load(open(args.cp_set)).items()})
    if args.nexp > 0.0:
        cal["NEXP"] = f"{args.nexp:.17g}"

    # t_end follows from the strain and the imposed rate, instead of being a
    # fixed number in the template. With t_end = 10 s a 17% ramp runs at
    # 1.7e-2 /s, 17x the 1e-3 /s the parameters were fitted at, which at
    # n = 8.33 is a 1.4x stress difference -- and drives the overstress that
    # exhausts the CP sub-step budget.
    t_end = args.strain / args.strain_rate

    deck = build_deck(mesh=mesh_e, orient=args.orientation, gc=args.gc, l0=args.l0,
                      nsteps=args.load_steps, mechsub=mechsub,
                      maxalt=args.alternations, strain=args.strain,
                      out_dir=out_dir, cal=cal,
                      maxcutbacks=args.max_cutbacks, t_end=t_end, bc=args.bc)

    cfg = load_input(str(deck))
    grain_info = None
    if args.grains > 0 or args.equiaxed > 0.0:
        # Columnar polycrystal: the tessellation is done on the in-plane
        # element centroids, so every grain spans the full thickness. With one
        # element through z that is the only representation available, and it
        # is NOT the experiment's equiaxed 3-D grain structure -- a difference
        # that belongs in the calibration's caveats, not buried here.
        cent = nodes[cells].mean(axis=1)
        if args.equiaxed > 0.0:
            extent = (float(nodes[:, 0].max()), float(nodes[:, 1].max()))
            gid, seeds = _equiaxed_grains(cent[:, :2], args.equiaxed, extent,
                                          args.grain_seed)
            n_grains = int(gid.max()) + 1
        else:
            gid, seeds = _voronoi_grains(cent[:, :2], args.grains, args.grain_seed)
            n_grains = args.grains
        eul = _random_orientations(n_grains, args.orientation_seed)
        n_gauss = 8
        euler_qp = np.repeat(eul[gid], n_gauss, axis=0)
        cfg["orientation_field"] = {"euler_angles": euler_qp}
        counts = np.bincount(gid, minlength=args.grains)
        grain_info = {"n_grains": int(args.grains),
                      "grain_seed": args.grain_seed,
                      "orientation_seed": args.orientation_seed,
                      "elements_per_grain_min": int(counts.min()),
                      "elements_per_grain_max": int(counts.max()),
                      "elements_per_grain_mean": float(counts.mean()),
                      "structure": "columnar (one element through thickness)"}

    t0 = time.perf_counter()
    result = Runner.from_config(cfg)
    wall = time.perf_counter() - t0

    # Diagnostics kept for BOTH outcomes. A failed run is when the
    # alternation history is most wanted, and re-running to obtain it costs
    # what the run cost. ``converged``/``failed_steps`` are the top-level
    # flags the mechanics executioner now reports; ``staggered`` is the
    # per-load-step alternation record.
    diag = {
        "staggered": _find(result, "staggered"),
        # NOT ``_find(result, "history").get("mechanics")``: _find returns the
        # FIRST "history" it meets, which is the runner's outer wrapper, and
        # the mechanics list lives one level further in
        # (history -> history -> mechanics). Asking the wrapper for
        # "mechanics" yields None, which then reads as "the solver reported
        # nothing" when in fact the record was right there.
        "target_strain": args.strain,
        "strain_rate_per_s": args.strain_rate,
        "mech_history": _find_list(result, "mechanics"),
        "mechanics_cutbacks_max": _find(result, "mechanics_cutbacks_max"),
        "mechanics_converged": _find(result, "mechanics_converged"),
        # The coupled executioner does not propagate the mechanics
        # sub-result's top-level keys, so these are frequently absent. They
        # are recorded as "not_reported" rather than None, because a null
        # here would read as "nothing failed" -- the exact confusion G25 is
        # about. The authoritative check remains the physical invariant
        # below, which does not depend on any status flag.
        "mech_converged": _nvl(_find(result, "converged")),
        "mech_failed_steps": _nvl(_find(result, "failed_steps")),
    }
    snaps = _find(result, "snapshots") or []
    summary = write_run_results(
        out_dir, snaps, nodes, cells,
        boundary_nodes=np.asarray(z["set_xmax"], dtype=int),
        gauge_length=PLATE_LENGTH, load_axis=0, write_vtk=not args.no_vtk,
    )
    _early = dict(diag)
    _early.update({"tag": tag, "n_elements": int(cells.shape[0]),
                   "n_nodes": int(nodes.shape[0]), "wall_seconds": wall,
                   "steps_written": summary.get("n_steps"),
                   "status": "pending_validation"})
    (out_dir / "run.json").write_text(json.dumps(_early, indent=1, default=str))

    # VERIFY THE BOUNDARY CONDITION WAS ACTUALLY APPLIED (study gap G25).
    #
    # Checked against the physical invariant rather than a status flag. A
    # mechanics step that hits the cutback floor is recorded
    # ``converged: False`` and then ABANDONED -- step_idx advances, the load
    # increment is never solved, and the run returns success at the top level
    # with a displacement field that does not satisfy its own Dirichlet
    # condition. Guarding on the reported flag does not work: with a single
    # alternation the staggered record is ``converged: False`` by
    # construction, so it cannot distinguish a lagged coupling from an
    # abandoned load step. The prescribed displacement is unambiguous.
    import csv as _csv

    with (out_dir / "results.csv").open() as fh:
        _rows = list(_csv.DictReader(fh))
    bad = []
    for _r in _rows:
        want = args.strain * float(_r["load_factor"])
        got = float(_r["boundary_strain"])
        if abs(got - want) > max(1e-6, 1e-3 * abs(want)):
            bad.append(int(_r["step"]))

    # Distinguish a run that FAILED from one that FINISHED FRACTURING.
    #
    # Once the specimen is severed it carries no load and is a mechanism:
    # there is no equilibrium state at a larger prescribed displacement, so
    # the solve legitimately has nowhere to go. That is a terminal condition
    # of the physics, not a defect, and the load steps before it are a valid
    # fracture curve that must not be thrown away. The discriminator is the
    # load-carrying capacity at the last good step, measured against the peak
    # the run reached: a genuine failure aborts while the specimen is still
    # carrying load.
    status = "complete"
    valid = _rows
    init = summary.get("initiation")
    if bad:
        first_bad = min(bad)
        valid = [r for r in _rows if int(r["step"]) < first_bad]
        peak = max((abs(float(r["avg_stress_axial"])) for r in _rows), default=0.0)
        last_ok = abs(float(valid[-1]["avg_stress_axial"])) if valid else float("inf")
        severed = bool(valid) and peak > 0 and last_ok <= 0.05 * peak
        # Did damage significantly initiate BEFORE the solve gave up?
        #
        # On a coarse mesh the load increment the solver can carry collapses
        # once a damage band forms. A run that stops shortly after initiation
        # is therefore the expected outcome rather than a failure, and the
        # practice is to report initiation instead of insisting on
        # propagation the discretization cannot carry. Such a run is a VALID
        # initiation measurement and its pre-failure steps are kept. A run
        # that quits BEFORE initiation measured nothing and is rejected.
        initiated = init is not None and int(init["step"]) < first_bad
        prop = summary.get("propagation")
        propagated = bool(prop and prop.get("propagated"))
        # Why did it stop? "Damage exists" is not a reason to stop -- a
        # resolved mesh should carry propagation. The reason is the load
        # increment collapsing under repeated cutback, which the mechanics
        # history now reports per step as increment_fraction.
        mech_hist = diag.get("mech_history") or []
        depths = [int(h.get("cutbacks", 0)) for h in mech_hist if h.get("converged")]
        if not mech_hist:
            # No mechanics history means the question was not answered, not
            # that the increment held up. Saying "no deep cutback" here would
            # report an absence of evidence as evidence of absence.
            stop_reason = "load-step behavior not reported by the solver"
            worst = -1
        else:
            worst = max(depths) if depths else 0
        stop_reason = stop_reason if not mech_hist else (
            "load increment collapsed under cutback "
            f"(deepest {worst} halvings = 1/{2 ** worst} of the nominal step)"
            if worst >= 3 else
            "solve failed without deep cutback -- not an increment-size limit"
        )
        if severed or initiated:
            status = ("fractured" if severed
                      else ("propagated" if propagated else "initiated"))
            with (out_dir / "results.csv").open("w", newline="") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(valid[0].keys()))
                w.writeheader()
                w.writerows(valid)
            if severed:
                (out_dir / "FRACTURED.txt").write_text(
                    f"Terminated at complete fracture after step {first_bad - 1}. "
                    f"Axial stress fell to {last_ok:.4g} MPa from a peak of "
                    f"{peak:.4g} MPa ({100 * last_ok / peak:.2f}%).\n"
                )
            else:
                (out_dir / f"{status.upper()}.txt").write_text(
                    (f"Terminated after damage PROPAGATION: the damaged fraction grew "
                     f"from {prop['damaged_fraction_at_initiation']:.4g} to "
                     f"{prop['damaged_fraction_at_last_step']:.4g} over "
                     f"{prop['steps_after_initiation']} steps and "
                     f"{prop['strain_span_after_initiation']:.5f} strain after "
                     f"initiation.\nStop reason: {stop_reason}.\n\n"
                     if propagated else
                     "Terminated after damage initiation, the expected outcome on a "
                    "mesh this coarse: the solvable load increment collapses once a "
                    "damage band forms.\n\n"
                    f"Initiation: step {init['step']}, strain "
                    f"{init['boundary_strain']:.5f}, axial stress "
                    f"{init['avg_stress_axial_MPa']:.4g} MPa, damaged fraction "
                    f"{init['damage_volume_fraction']:.4g}.\n"
                    f"Criterion: {init['criterion']}.\n"
                    f"Solve stopped at step {first_bad}, "
                    f"{first_bad - int(init['step'])} step(s) after initiation, "
                    f"still carrying {last_ok:.4g} MPa.\n\n"
                    f"Steps 0-{first_bad - 1} retained. INITIATION metrics are "
                    "usable; PROPAGATION metrics are not.\n")
                    + (f"\nStop reason: {stop_reason}.\n" if not propagated else "")
                )
        else:
            (out_dir / "INVALID.txt").write_text(
                f"Failed at step {first_bad} carrying {last_ok:.4g} MPa "
                f"({100 * last_ok / max(peak, 1e-30):.1f}% of peak {peak:.4g} MPa) "
                "and BEFORE damage initiation, so it measured nothing.\n"
            )
            raise SystemExit(
                f"{tag}: the solve failed at step {first_bad} BEFORE damage "
                "initiated (element nodal-average damage never reached 0.85), "
                f"while still carrying {last_ok:.4g} MPa "
                f"({100 * last_ok / max(peak, 1e-30):.1f}% of peak). Nothing was "
                "measured. See G25."
            )

    rec = {
        "tag": tag, "geometry": args.geometry, "density": args.density,
        "n_elements": int(cells.shape[0]), "n_nodes": int(nodes.shape[0]),
        "orientation": args.orientation, "Gc_MPa_um": args.gc, "l0_um": args.l0,
        "target_strain": args.strain,
        "strain_rate_per_s": args.strain_rate,
        "t_end_s": t_end, "load_steps": args.load_steps,
        "mech_substeps": mechsub, "alternations": args.alternations,
        "strain_increment_per_update": d_eps / mechsub,
        "wall_seconds": wall,
        "diagnostics": diag,
        "grains": grain_info,
        "cp_set": cal,
        "initiation": init,
        "propagation": summary.get("propagation"),
        "stop_reason": stop_reason if bad else "ran to the target strain",
        "status": status,
        "steps_retained": len(valid),
        "vtk_series": str(out_dir / "fields" / "fields.pvd"),
        **summary,
    }
    (out_dir / "run.json").write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: rec[k] for k in
                      ("tag", "n_elements", "mech_substeps", "n_steps",
                       "vtk_files", "wall_seconds")}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
