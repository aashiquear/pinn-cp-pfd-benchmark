#!/usr/bin/env python3
"""Round a Voronoi polycrystal by threshold dynamics (MBO), no phase field.

WHY THIS IS NOT AD HOC SMOOTHING. Merriman-Bence-Osher threshold dynamics is a
provable approximation to multiphase mean-curvature flow: diffuse each grain's
indicator function for a short time tau, then reassign every point to whichever
grain has the largest diffused indicator. Iterating that converges to
curvature flow with the correct Herring condition at triple junctions, i.e.
120-degree angles for isotropic boundary energy. So it produces the same
morphology the phase-field relaxation is supposed to produce, by a different
and much cheaper route -- and with no mobility, interface width or time-step
parameter to get wrong.

WHAT IT IS NOT. It carries no kinetics: the diffusion length sqrt(2 tau) sets
how far boundaries move per sweep, but nothing maps a sweep to physical time.
It is a MORPHOLOGY GENERATOR only, which is exactly the stated purpose here,
and it must not be reported as grain growth.

The one parameter, the diffusion length per sweep, is bounded by the physics
it approximates: it must be small against the grain size (or curvature flow is
not resolved) and large against the grid spacing (or nothing moves). Grains
smaller than the diffusion length are consumed, which is the correct behaviour
of curvature flow, so the grain count is recorded per sweep rather than
assumed constant.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grain_morphology import summarize  # noqa: E402
from map_grains import voronoi_map  # noqa: E402

STUDY = Path(__file__).resolve().parent.parent
PLATE = (100.0, 50.0)


def mbo_sweep(gid: np.ndarray, sigma_px: tuple[float, float], n_grains: int):
    """One MBO step: diffuse each indicator, then take the argmax.

    Gaussian filtering IS the heat kernel, so a Gaussian of standard deviation
    s applies diffusion for tau = s**2 / 2. Reflecting boundaries keep grains
    that touch the domain edge from being eroded by an implicit outside phase.
    """
    best = np.full(gid.shape, -np.inf)
    out = np.zeros_like(gid)
    for g in range(n_grains):
        ind = (gid == g).astype(float)
        if not ind.any():
            continue
        d = ndimage.gaussian_filter(ind, sigma=sigma_px, mode="nearest")
        m = d > best
        best[m] = d[m]
        out[m] = g
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, required=True)
    ap.add_argument("--grid", type=int, nargs=2, default=[200, 100])
    ap.add_argument("--sweeps", type=int, default=12)
    # Diffusion length per sweep as a fraction of the mean grain width. Small
    # enough that curvature flow is resolved, large enough to move boundaries.
    ap.add_argument("--length-frac", type=float, default=0.06)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    div = (args.grid[0], args.grid[1])
    hx, hy = PLATE[0] / div[0], PLATE[1] / div[1]
    seeds = np.asarray(np.load(STUDY / f"meshes/plate/poly_d3_n{args.grains}.npz")["seeds"],
                       dtype=float)
    gid = voronoi_map(seeds, div)
    n = int(gid.max()) + 1

    grain_w = float(np.sqrt(PLATE[0] * PLATE[1] / n))
    L = args.length_frac * grain_w           # diffusion length per sweep, um
    sig = (L / hx, L / hy)                   # in pixels, per axis
    if min(sig) < 1.0:
        raise SystemExit(
            f"diffusion length {L:.3f} um is under one grid cell "
            f"({min(hx, hy):.3f} um): nothing would move. Refine the grid or "
            "raise --length-frac.")

    traj = [dict(sweep=0, **{k: v for k, v in summarize(gid, hx, hy, "sweep0").items()
                             if k in ("n_grains", "n_junctions",
                                      "junction_rms_dev_from_120",
                                      "junction_frac_within_15deg",
                                      "total_boundary_length_um",
                                      "circularity_mean")})]
    hist = [gid.copy()]
    for s in range(1, args.sweeps + 1):
        gid = mbo_sweep(gid, sig, n)
        st = summarize(gid, hx, hy, f"sweep{s}")
        traj.append(dict(sweep=s, **{k: v for k, v in st.items()
                                     if k in ("n_grains", "n_junctions",
                                              "junction_rms_dev_from_120",
                                              "junction_frac_within_15deg",
                                              "total_boundary_length_um",
                                              "circularity_mean")}))
        hist.append(gid.copy())

    out = Path(args.out or STUDY / f"results/rounded/n{args.grains}")
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "sweeps.npz", **{f"s{i}": g for i, g in enumerate(hist)},
                        seeds=seeds, div=np.array(div))
    rec = dict(n_grains_requested=args.grains, grid=list(div),
               diffusion_length_um=L, diffusion_length_px=list(sig),
               mean_grain_width_um=grain_w, sweeps=args.sweeps,
               length_frac=args.length_frac, trajectory=traj)
    (out / "rounding.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: rec[k] for k in
                      ("grid", "diffusion_length_um", "mean_grain_width_um")},
                     indent=1, default=float))
    print(f"{'sweep':>6} {'grains':>7} {'RMS dev 120':>12} {'within15':>9} {'GB len':>9} {'circ':>8}")
    for t in traj:
        print(f"{t['sweep']:>6} {t['n_grains']:>7} {t['junction_rms_dev_from_120']:11.2f}  "
              f"{100*t['junction_frac_within_15deg']:8.1f}% "
              f"{t['total_boundary_length_um']:9.1f} "
              f"{t.get('circularity_mean', float('nan')):8.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
