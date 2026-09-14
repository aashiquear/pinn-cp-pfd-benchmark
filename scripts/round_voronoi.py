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
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grain_morphology import summarize  # noqa: E402
from map_grains import voronoi_map  # noqa: E402

STUDY = Path(__file__).resolve().parent.parent
PLATE = (100.0, 50.0)


def padded_seeds(n_target: int, overseed: float, pad_um: float,
                 seed: int = 0) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Seeds over a domain PADDED beyond the plate, at the plate's density.

    Two reasons for the pad, both physical rather than numerical.

    * The plate is a SECTION cut from a larger polycrystal. Grains at its edge
      have neighbours outside it, and boundaries meet the cut at whatever angle
      the outside network dictates. Reflecting the diffusion at the edge
      instead forces those boundaries to meet it at 90 degrees, which is what
      makes edge grains look wrong -- and it dominates at high grain counts,
      where most grains touch an edge.
    * Grains are consumed by curvature flow, so the count inside the plate
      FALLS. Seeding above the target absorbs that, which is the usual practice
      when a target count matters.

    Seeds sit on a jittered lattice rather than at random: a Poisson process
    gives a broad size distribution, which is not the equiaxed structure
    wanted here.
    """
    rng = np.random.default_rng(seed)
    x0, y0 = -pad_um, -pad_um
    x1, y1 = PLATE[0] + pad_um, PLATE[1] + pad_um
    area_plate = PLATE[0] * PLATE[1]
    n_plate = int(round(n_target * overseed))
    density = n_plate / area_plate
    n_pad = max(1, int(round(density * (x1 - x0) * (y1 - y0))))
    # jittered lattice over the padded box, aspect-matched
    ar = (x1 - x0) / (y1 - y0)
    ny = max(1, int(round(np.sqrt(n_pad / ar))))
    nx = max(1, int(round(n_pad / ny)))
    gx = (np.arange(nx) + 0.5) / nx * (x1 - x0) + x0
    gy = (np.arange(ny) + 0.5) / ny * (y1 - y0) + y0
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    jx, jy = 0.38 * (x1 - x0) / nx, 0.38 * (y1 - y0) / ny
    pts = np.column_stack([
        (X + rng.uniform(-jx, jx, X.shape)).ravel(),
        (Y + rng.uniform(-jy, jy, Y.shape)).ravel()])
    return pts, (x0, y0, x1, y1)


def voronoi_on_box(seeds: np.ndarray, div: tuple[int, int],
                   box: tuple[float, float, float, float]) -> np.ndarray:
    x0, y0, x1, y1 = box
    xs = np.linspace(x0, x1, div[0] + 1)
    ys = np.linspace(y0, y1, div[1] + 1)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return cKDTree(seeds).query(np.column_stack([X.ravel(), Y.ravel()]))[1].reshape(
        div[0] + 1, div[1] + 1)


def crop_to_plate(gid: np.ndarray, box: tuple[float, float, float, float]):
    """Cut the plate out of the padded field and relabel 0..k-1."""
    x0, y0, x1, y1 = box
    nx, ny = gid.shape
    xs = np.linspace(x0, x1, nx); ys = np.linspace(y0, y1, ny)
    ix = np.where((xs >= 0.0) & (xs <= PLATE[0]))[0]
    iy = np.where((ys >= 0.0) & (ys <= PLATE[1]))[0]
    sub = gid[np.ix_(ix, iy)]
    ids = np.unique(sub)
    remap = {int(v): i for i, v in enumerate(ids)}
    out = np.vectorize(remap.get)(sub).astype(int)
    return out, len(ids)


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
        # "nearest" is correct HERE because the pad, not the plate edge, is
        # where it applies -- the plate interior never sees it.
        d = ndimage.gaussian_filter(ind, sigma=sigma_px, mode="nearest")
        m = d > best
        best[m] = d[m]
        out[m] = g
    return out


def edge_band_junctions(gid, hx, hy, band_um=8.0):
    """Junction-angle spread for junctions NEAR the plate edge, and interior.

    The edge grains are what a reflecting diffusion condition distorts, so they
    need their own measure: a domain-wide average hides a defect confined to
    one band.
    """
    from grain_morphology import triple_junction_angles
    j = triple_junction_angles(gid, hx, hy)
    a = np.asarray(j.get("angles") or [], dtype=float)
    xy = np.asarray(j.get("angle_xy") or [], dtype=float)
    if a.size == 0 or xy.size == 0:
        return {}
    near = ((xy[:, 0] < band_um) | (xy[:, 0] > PLATE[0] - band_um) |
            (xy[:, 1] < band_um) | (xy[:, 1] > PLATE[1] - band_um))
    out = {}
    for lab, m in (("edge", near), ("interior", ~near)):
        out[f"{lab}_n"] = int(m.sum())
        if m.sum():
            out[f"{lab}_rms_deg"] = float(np.sqrt(np.mean((a[m] - 120.0) ** 2)))
            out[f"{lab}_within15_frac"] = float(np.mean(np.abs(a[m] - 120.0) <= 15.0))
    return out


def substantive_count(gid, hx, hy, n_target):
    """Grains with real area in the plate, and the sliver count separately.

    A section cuts pad grains, leaving fragments with a few cells of area.
    Those are genuine grains of the sectioned polycrystal and each still needs
    an orientation, but counting them as equiaxed grains would misreport the
    structure. The threshold is a quarter of the target mean grain area.
    """
    cell = hx * hy
    mean_area = PLATE[0] * PLATE[1] / max(1, n_target)
    ids, cnt = np.unique(gid, return_counts=True)
    area = cnt * cell
    keep = area >= 0.25 * mean_area
    return dict(n_total=int(len(ids)), n_substantive=int(keep.sum()),
                n_slivers=int((~keep).sum()),
                sliver_area_frac=float(area[~keep].sum() /
                                       (PLATE[0] * PLATE[1])) if (~keep).any() else 0.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, required=True,
                    help="TARGET grain count inside the plate after rounding.")
    ap.add_argument("--overseed", type=float, default=1.0,
                    help="Seed this multiple of the target inside the plate, to "
                         "absorb grains consumed by curvature flow. 1.0 = none.")
    ap.add_argument("--pad-grains", type=float, default=1.0,
                    help="Buffer outside the plate, in mean grain widths. The "
                         "plate is a SECTION of a larger polycrystal, so its "
                         "edge grains must have real neighbours rather than a "
                         "reflecting wall.")
    ap.add_argument("--grid", type=int, nargs=2, default=[320, 160])
    ap.add_argument("--sweeps", type=int, default=24)
    ap.add_argument("--length-frac", type=float, default=0.06)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    grain_w = float(np.sqrt(PLATE[0] * PLATE[1] / args.grains))
    pad = args.pad_grains * grain_w
    seeds, box = padded_seeds(args.grains, args.overseed, pad, seed=args.seed)
    # grid over the PADDED box at the plate's requested resolution
    px = args.grid[0] / PLATE[0]
    div = (int(round((box[2] - box[0]) * px)), int(round((box[3] - box[1]) * px)))
    hx = (box[2] - box[0]) / div[0]
    hy = (box[3] - box[1]) / div[1]
    L = args.length_frac * grain_w
    sig = (L / hx, L / hy)
    if min(sig) < 1.0:
        raise SystemExit(
            f"diffusion length {L:.3f} um is under one grid cell "
            f"({min(hx, hy):.3f} um): nothing would move.")

    gid = voronoi_on_box(seeds, div, box)
    n = int(gid.max()) + 1

    def record(sweep, g_full):
        g, k = crop_to_plate(g_full, box)
        st = summarize(g, hx, hy, f"sweep{sweep}")
        eb = edge_band_junctions(g, hx, hy)
        sc = substantive_count(g, hx, hy, args.grains)
        return dict(sweep=sweep, grains_in_plate=k, **sc,
                    **{kk: st[kk] for kk in
                       ("n_junctions", "junction_rms_dev_from_120",
                        "junction_frac_within_15deg",
                        "total_boundary_length_um", "circularity_mean")
                       if kk in st},
                    **eb), g

    traj = []; crops = []
    t, g = record(0, gid); traj.append(t); crops.append(g)
    for s_i in range(1, args.sweeps + 1):
        gid = mbo_sweep(gid, sig, n)
        t, g = record(s_i, gid); traj.append(t); crops.append(g)

    out = Path(args.out or STUDY / f"results/rounded/n{args.grains}_pad")
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "sweeps.npz",
                        **{f"s{i}": g for i, g in enumerate(crops)},
                        seeds=seeds, box=np.array(box),
                        div=np.array(div))
    rec = dict(n_target=args.grains, overseed=args.overseed,
               n_seeded_padded=int(len(seeds)), pad_um=pad,
               grid_padded=list(div), grid_plate=list(args.grid),
               diffusion_length_um=L, mean_grain_width_um=grain_w,
               length_frac=args.length_frac, sweeps=args.sweeps,
               trajectory=traj)
    (out / "rounding.json").write_text(json.dumps(rec, indent=1, default=float))
    print(f"target {args.grains} grains in plate | seeded {len(seeds)} over a "
          f"padded box (pad {pad:.1f} um) | L = {L:.3f} um")
    print(f"{'sweep':>6} {'total':>6} {'real':>5} {'sliv':>5} {'RMS':>7} {'edgeRMS':>8} "
          f"{'intRMS':>7} {'within15':>9} {'circ':>7} {'GB len':>8}")
    for t in traj:
        print(f"{t['sweep']:>6} {t['n_total']:>6} {t['n_substantive']:>5} {t['n_slivers']:>5} "
              f"{t['junction_rms_dev_from_120']:7.2f} "
              f"{t.get('edge_rms_deg', float('nan')):8.2f} "
              f"{t.get('interior_rms_deg', float('nan')):7.2f} "
              f"{100*t['junction_frac_within_15deg']:8.1f}% "
              f"{t.get('circularity_mean', float('nan')):7.4f} "
              f"{t['total_boundary_length_um']:8.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
