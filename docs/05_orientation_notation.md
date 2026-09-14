# Orientation notation and the Fig. 5 IPF marker

## Notation

The fourth Al2219 single-crystal orientation is written **`[-1 4 9]`** in
(Rezwan 2025 intragranular) and in this study's decks and scripts. The author
confirms it denotes **[1&#772; 4 9]** — bar over the 1, i.e. h = -1 — and that
the published text had no way to typeset the overbar.

**This is a typography note, not a numerical one.** 1&#772; = -1, so the Miller
index is (-1, 4, 9) either way, which is the value every calculation in this
study has used. `loading_direction = [-1.0, 4.0, 9.0]` in the plate decks is
correct as written and needs no change.

American English and the ASCII spelling `[-1 4 9]` are retained throughout the
study files, since the overbar does not survive plain-text decks, CSV headers,
or filenames.

## The bar placement cannot explain the Fig. 5 marker discrepancy

A reasonable hypothesis was that lost overbar typography might mean the bar
belonged on a different component. It does not help: an inverse pole figure
plots the *symmetry representative* of a direction, and for cubic symmetry that
representative is |components| sorted descending. Every one of the eight sign
variants of (1, 4, 9) therefore lands at exactly the same point.

| index | triangle representative | (X, Y) | angle from <001> |
|---|---|---|---|
| all 8 of [±1 ±4 ±9] | (0.404061, 0.101015, 0.909137) | (0.211646, 0.052911) | 24.6136 deg |

Distinct IPF positions among all eight variants: **1**.

The digitized marker sits at (X, Y) = (0.111330, 0.050834), i.e. 13.955 deg
from <001> — **28.60 px** from where any member of the (1, 4, 9) family belongs,
against a ~1 px centroid uncertainty on a 9 px marker.

So the discrepancy is in the *magnitudes*, not the signs: the drawn marker
corresponds to magnitudes near (2.19, 1, 9.69), not (4, 1, 9). No choice of
overbar placement, and no sign convention, can move the marker onto [1&#772; 4 9].

## What remains open

Two readings are consistent with the evidence, and the published figure alone
cannot separate them:

1. The simulation used the Euler triple (168.5740, -14.032, 24.6187), which
   under PINNACLE's active Bunge ZXZ crystal-to-sample convention denotes a
   direction 14.032 deg from <001> — matching the drawn marker's polar angle to
   0.08 deg. Then the fourth case is not [1&#772; 4 9] and the label is wrong.
2. The simulation used [1&#772; 4 9] and only the IPF marker was plotted from a
   mis-converted triple. Then the label is right and the figure's dot is
   misplaced.

The two codes' triples for this orientation differ by a one-slot shift, which is
the likely origin of the mis-conversion in either reading:

| | phi1 | Phi | phi2 |
|---|---|---|---|
| PINNACLE, [1&#772; 4 9] about z | **-14.0362** | **24.6136** | 14.0362 |
| MOOSE triple as supplied | 168.5740 | **-14.032** | **24.6187** |

The same magnitudes 14.03 and 24.61 appear in both, shifted one position.

**Consequence for the benchmark.** Which reading holds decides whether a
PINNACLE run configured with `loading_direction = [-1.0, 4.0, 9.0]` reproduces
the paper's fourth case or a different crystal. The three corner orientations
are unaffected: [001], [101] and [111] are invariant under any axis relabeling
or sign convention. Only the author's original input deck can settle it.
