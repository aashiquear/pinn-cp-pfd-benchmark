#!/usr/bin/env python3
"""Turn a coupled run's per-step damage snapshots into evidence.

Every phase-field damage run in this study must ship visual evidence of
initiation AND propagation, because a converged scalar cannot tell localized
void-initiated fracture from diffuse damage over the whole domain -- a
distinction this study has already been caught by once, on a run whose damage
maximum sat 32.9 um from the void surface with three quarters of the nodes past
the initiation threshold.

Produces, from a snapshot list:

* ``damage_evolution.png`` -- a panel per recorded load step on the analysis
  slice, on a SHARED color scale so panels are comparable (a per-panel scale
  would make an early, faint field look identical to a fully damaged one);
* ``damage_evolution.gif`` -- the animation over the load path;
* ``damage_history.csv`` -- long-form (step, load_factor, node, x, y, z,
  damage) on the slice, so the animation can be rebuilt without re-running;
* ``damage_history.npz`` -- the full nodal fields at every step, since the CSV
  is restricted to the slice.

The CSV is the fallback the project guidelines require when an animation cannot
be produced; here it is written unconditionally, because it also carries the
data the figures were made from.
"""

from __future__ import annotations

import numpy as np


def _slice_indices(nodes: np.ndarray, axis: int, value: float, tol: float) -> np.ndarray:
    return np.nonzero(np.abs(nodes[:, axis] - value) <= tol)[0]


def write_evidence(
    snapshots: list[dict],
    nodes: np.ndarray,
    out_dir,
    *,
    label: str = "",
    slice_axis: int = 2,
    slice_value: float | None = None,
    void_centre: np.ndarray | None = None,
    void_radius: float | None = None,
) -> dict:
    """Write image, animation and history files. Returns a manifest dict."""
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    nodes = np.asarray(nodes, dtype=float).reshape(-1, 3)

    fields = [np.asarray(s["damage"], dtype=float).ravel() for s in snapshots]
    steps = [int(s["step"]) for s in snapshots]
    lfs = [float(s["step_factor"]) for s in snapshots]
    if not fields or fields[0].size != nodes.shape[0]:
        raise SystemExit(
            f"damage field size {fields[0].size if fields else 0} does not match "
            f"{nodes.shape[0]} mesh nodes; refusing to plot a mismatched field"
        )

    # Slice through the void centre by default: that is where initiation is
    # expected, and a slice elsewhere could miss the band entirely.
    if slice_value is None:
        slice_value = float(void_centre[slice_axis]) if void_centre is not None \
            else float(np.median(nodes[:, slice_axis]))
    coords = np.unique(nodes[:, slice_axis])
    nearest = float(coords[np.argmin(np.abs(coords - slice_value))])
    tol = 1.0e-6 * max(1.0, np.ptp(nodes[:, slice_axis]))
    idx = _slice_indices(nodes, slice_axis, nearest, tol)
    if idx.size < 8:
        raise SystemExit(f"slice at {slice_axis}={nearest} holds only {idx.size} nodes")

    ax_ids = [a for a in (0, 1, 2) if a != slice_axis]
    px, py = nodes[idx, ax_ids[0]], nodes[idx, ax_ids[1]]
    names = "xyz"

    vmax = max(1.0e-12, max(float(f.max()) for f in fields))
    n = len(fields)

    # --- panel figure ---
    ncol = min(n, 5)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.0 * nrow),
                             squeeze=False, constrained_layout=True)
    sc = None
    for k in range(nrow * ncol):
        ax = axes[k // ncol][k % ncol]
        if k >= n:
            ax.axis("off")
            continue
        sc = ax.tripcolor(px, py, fields[k][idx], vmin=0.0, vmax=vmax,
                          cmap="inferno", shading="gouraud")
        if void_centre is not None and void_radius is not None:
            r2 = void_radius**2 - (nearest - void_centre[slice_axis]) ** 2
            if r2 > 0:
                ax.add_patch(plt.Circle(
                    (void_centre[ax_ids[0]], void_centre[ax_ids[1]]), np.sqrt(r2),
                    fill=False, ec="cyan", lw=1.0))
        ax.set_title(f"step {steps[k]}  load {lfs[k]:.2f}", fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{names[ax_ids[0]]} (um)", fontsize=8)
        ax.set_ylabel(f"{names[ax_ids[1]]} (um)", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.colorbar(sc, ax=axes, shrink=0.8, label="damage d")
    fig.suptitle(f"Damage initiation and propagation{': ' + label if label else ''}"
                 f"  (slice {names[slice_axis]} = {nearest:g} um)", fontsize=11)
    png = out_dir / "damage_evolution.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    # --- animation ---
    gif = out_dir / "damage_evolution.gif"
    animated = True
    try:
        figa, axa = plt.subplots(figsize=(4.4, 4.0), constrained_layout=True)
        tp = axa.tripcolor(px, py, fields[0][idx], vmin=0.0, vmax=vmax,
                           cmap="inferno", shading="gouraud")
        figa.colorbar(tp, ax=axa, label="damage d")
        axa.set_aspect("equal")
        axa.set_xlabel(f"{names[ax_ids[0]]} (um)")
        axa.set_ylabel(f"{names[ax_ids[1]]} (um)")
        ttl = axa.set_title("")

        def _frame(k):
            axa.clear()
            axa.tripcolor(px, py, fields[k][idx], vmin=0.0, vmax=vmax,
                          cmap="inferno", shading="gouraud")
            if void_centre is not None and void_radius is not None:
                r2 = void_radius**2 - (nearest - void_centre[slice_axis]) ** 2
                if r2 > 0:
                    axa.add_patch(plt.Circle(
                        (void_centre[ax_ids[0]], void_centre[ax_ids[1]]), np.sqrt(r2),
                        fill=False, ec="cyan", lw=1.0))
            axa.set_aspect("equal")
            axa.set_xlabel(f"{names[ax_ids[0]]} (um)")
            axa.set_ylabel(f"{names[ax_ids[1]]} (um)")
            axa.set_title(f"step {steps[k]}   load {lfs[k]:.2f}   peak d {fields[k].max():.4f}")
            return ()

        FuncAnimation(figa, _frame, frames=n, blit=False).save(
            gif, writer=PillowWriter(fps=max(1, min(6, n // 2 or 1))))
        plt.close(figa)
    except Exception as exc:  # pragma: no cover - environment dependent
        animated = False
        (out_dir / "animation_not_produced.txt").write_text(
            f"Animation could not be produced: {type(exc).__name__}: {exc}\n"
            "The damage history CSV and NPZ beside this file carry the per-step "
            "fields needed to build it.\n"
        )

    # --- history files ---
    csv = out_dir / "damage_history.csv"
    with open(csv, "w") as fh:
        fh.write("# Per-load-step damage on the analysis slice "
                 f"({names[slice_axis]} = {nearest:g} um, {idx.size} nodes).\n")
        fh.write("# The FULL nodal field at every step is in damage_history.npz "
                 "(arrays: damage [n_steps, n_nodes], nodes [n_nodes, 3], "
                 "step, load_factor).\n")
        fh.write("step,load_factor,node_id,x_um,y_um,z_um,damage\n")
        for k in range(n):
            dk = fields[k]
            for j in idx:
                fh.write(f"{steps[k]},{lfs[k]:.6f},{j},{nodes[j,0]:.6f},"
                         f"{nodes[j,1]:.6f},{nodes[j,2]:.6f},{dk[j]:.8e}\n")
    npz = out_dir / "damage_history.npz"
    np.savez_compressed(npz, damage=np.vstack(fields), nodes=nodes,
                        step=np.asarray(steps), load_factor=np.asarray(lfs))

    return {
        "image": str(png), "animation": str(gif) if animated else None,
        "animation_produced": animated, "history_csv": str(csv),
        "history_npz": str(npz), "n_steps": n,
        "slice_axis": names[slice_axis], "slice_value_um": nearest,
        "slice_nodes": int(idx.size), "damage_max_over_path": float(vmax),
    }
