#!/usr/bin/env python3
"""Case C6 (gallery of spherulitic morphologies, Table II of granasy 2005 spherulites): figures.

Usage: python3 scripts/cryst_gallery_figures.py
Reads results/crystallization/orientation/gallery/c6_fig10<panel>_fem/frames.npz (from
cryst_binary_collect.py); writes figures/crystallization/c6_gallery.png (orientation map and
composition of every finished panel at the last frame before the crystal reaches the box
edge), c6_gallery_<panel>.gif and results/crystallization/orientation/c6_gallery_measures.csv
(the measures of cryst_orientation_series.py, two-fold symmetry).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cryst_orientation_decks import GALLERY  # noqa: E402
from cryst_orientation_figures import FIG, ROOT, orientation_image  # noqa: E402
from cryst_orientation_series import last_free_frame, measures  # noqa: E402

RES = ROOT / "results/crystallization/orientation"


def main() -> None:
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 200, "axes.linewidth": 0.8, "mathtext.fontset": "cm"})
    panels = [p for p in GALLERY if (RES / f"gallery/c6_fig10{p}_fem/frames.npz").exists()]
    cols = (len(panels) + 1) // 2
    fig, axes = plt.subplots(4, cols, figsize=(2.5 * cols, 10.4), squeeze=False)
    for a in axes.ravel():
        a.set_xticks([]); a.set_yticks([])
        a.set_frame_on(False)
    rows = []
    for j, p in enumerate(panels):
        s, s0, delta0, x, angle, m_phi, chi, _ = GALLERY[p]
        d = np.load(RES / f"gallery/c6_fig10{p}_fem/frames.npz")
        i = last_free_frame(d["phi"])
        m, _ = measures(d, 2)
        rows.append({"panel": p, "S": s, "s0": s0, "delta0": delta0, "x": x, "angle_deg": angle or 0,
                     "m_phi": m_phi, "chi_l": chi, **m})
        r, c = 2 * (j // cols), j % cols
        axes[r, c].imshow(orientation_image(d["phi"][i], d["theta"][i]), origin="lower")
        axes[r, c].set_title(f"({p}) $S$ = {s:g}, $s_0$ = {s0:g}, $\\delta_0$ = {delta0:g}\n"
                             f"$x$ = {x:g}" + (f" ({angle}°)" if angle else "") + f", $\\chi_l$ = {chi:g}", fontsize=7)
        axes[r + 1, c].imshow(d["c"][i], origin="lower", cmap="viridis", vmin=float(d["c_s"]) - 0.005,
                              vmax=float(d["c_l"]) + 0.005)
        frames = [orientation_image(d["phi"][k], d["theta"][k]) for k in range(i + 1)]
        gfig, gax = plt.subplots(figsize=(3.2, 3.2))
        gax.set_xticks([]); gax.set_yticks([])
        im = gax.imshow(frames[0], origin="lower")
        FuncAnimation(gfig, lambda k, im=im, frames=frames: (im.set_data(frames[k]),), frames=len(frames)).save(
            FIG / f"c6_gallery_{p}.gif", writer=PillowWriter(fps=4), dpi=80)
        plt.close(gfig)
    fig.tight_layout()
    fig.savefig(FIG / "c6_gallery.png", bbox_inches="tight")
    plt.close(fig)
    with (RES / "c6_gallery_measures.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(" ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
    print("wrote figures/crystallization/c6_gallery.png")


if __name__ == "__main__":
    main()
