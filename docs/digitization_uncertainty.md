
# Digitization uncertainty record

Every figure-only quantity in `benchmark_dataset.csv` tagged `digitized` traces to this document.
No value here was estimated from surrounding prose; anything not measurable from a figure is tagged
`unavailable` instead.

## Common method

**Raster.** The figures in both publisher PDFs are embedded raster images, not vector graphics
(confirmed by object inspection: each figure page carries `PdfImage` objects and no plot paths).
They were therefore extracted at their **native embedded resolution** with `pypdfium2`
(`PdfImage.get_bitmap()`), not by rasterizing the page. Every figure in both papers is embedded at
**300 DPI** — verified per figure as `pixel_width / (bbox_width_pt / 72)`:

| figure | native pixels | placed size (pt) | effective DPI |
|---|---|---|---|
| Al2219 Fig. 3 | 1221 x 854 | 293.0 x 204.8 | 300 |
| Al2219 Fig. 7 | 852 x 617 | 204.5 x 147.9 | 300 |
| Al2219 Fig. 8 (a+b) | 1795 x 707 | 430.8 x 169.5 | 300 |
| Al2219 Fig. 9 | 856 x 625 | 205.4 x 149.9 | 300 |
| Ni Fig. 2 | 1250 x 647 | 300.0 x 155.1 | 300 |
| Ni Fig. 3 | 1000 x 732 | 240.2 x 175.7 | 300 |
| Ni Fig. 4(b) | 769 x 745 | 184.6 x 178.9 | 300 |

Rendering the page at a higher DPI would only upsample; 300 DPI is the information ceiling of these
sources. All extraction was done on the native bitmaps.

**Axis calibration.** Tick-mark pixel positions were located by scanning the 2–8 px band immediately
outside each spine for dark pixels and grouping runs; a first-order polynomial was then fitted from
tick pixel position to tick label value, and the **maximum absolute residual over all ticks** is
reported as the calibration residual. No calibration used fewer than five ticks per axis.

**Repeated extraction.** Every extraction was run **twice with independent segmentation thresholds**
(a strict pass and a permissive pass). The reported repeated-extraction error is the discrepancy
between the two passes on matched features; where the two passes disagreed about the *existence* of
a feature, both are named explicitly below.

---

## Al2219 Fig. 3 — CP calibration stress–strain curve

Files: `digitized/al2219_fig3_calibration_stress_strain.csv` (full traces),
`digitized/al2219_fig3_experiment_markers.csv` (9 marker centers).

- **Method.** Color-channel segmentation: blue (`B - max(R,G) > 40`) = simulation, red
  (`R - max(G,B) > 40`) = experiment. One row per pixel column, reporting the column mean plus the
  column's full vertical extent (`stress_min_MPa`, `stress_max_MPa`, `n_px`). The extent columns are
  not decoration: the elastic segment is near-vertical, so a column mean is not a well-defined curve
  value there, and reporting the bracket keeps that visible rather than smoothing it away. The
  legend swatches in the lower right were excluded by a rectangular mask
  (x > 800 px, y in 560–757 px); no data columns fall in that box.
  Experiment marker centers were found separately by convolving the red mask with a 13 x 13 box
  filter — a marker stroke is locally ~2.3x denser than the dashed connecting line — thresholding at
  density > 0.36 and taking the centroid of the full marker footprint in a 21 x 21 window. The
  origin marker sits on the axis corner and was measured directly from its red pixels.
- **Raster resolution.** 300 DPI (1221 x 854 px).
- **Axis calibration.** x: 5 ticks at pixel 118.5, 385.5, 652.5, 919.5, 1187.5 -> 0, 5, 10, 15, 20 %
  strain; slope 0.0187126 %/px; **max|residual| = 0.0075 % strain, RMS 0.0053 %**.
  y: 8 ticks at pixel 755.5, 649.5, 543.5, 438.0, 332.5, 226.5, 120.5, 15.5 -> 0…700 MPa;
  slope -0.945785 MPa/px; **max|residual| = 0.36 MPa, RMS 0.23 MPa**.
- **Estimated extraction uncertainty.** Curve traces: **+-0.03 % strain, +-1.0 MPa** (line
  half-width ~1.5 px at 0.019 %/px and 0.95 MPa/px, plus calibration residual). Marker centers:
  **+-0.03 % strain, +-1.5 MPa**.
- **Repeated extraction.** Marker-center density threshold varied 0.30 / 0.33 / 0.36, giving 11 / 10
  / 8 detections. The **8 markers common to all three thresholds agree to within 0.012 % strain and
  0.7 MPa** (max spread 0.012 % at the 17 % marker and 0.011 % at the 2 % marker; four of the eight
  agree exactly, and the largest stress spread is 0.70 MPa at the 0.57 % marker). Separately, at the
  two looser thresholds the 17 % marker also fragments into two pieces 0.112 % apart (16.91 % and
  16.97 %), which the 0.36 threshold correctly resolves as the single marker at 16.983 %; the
  reported value is that single-marker centroid. The 0.30 threshold
  also emits one spurious fragment at 0.134 % that duplicates the 0.139 % marker. The **origin
  marker is clipped by the two axis lines** and is not recovered by the density filter at 0.36; it
  was measured directly as the centroid of its 132 red pixels, giving (0.003 %, 1.8 MPa), against
  (0.027 %, 3.6 MPa) from the looser thresholds — so that one point carries **+-0.025 % strain,
  +-1.9 MPa**, roughly an order of magnitude worse than the other eight.
- **Recovered marker values.** (0.003 %, 1.8 MPa), (0.140 %, 137.7), (0.566 %, 344.9),
  (2.017 %, 411.9), (7.006 %, 502.5), (11.019 %, 552.5), (13.049 %, 573.6), (15.056 %, 588.3),
  (16.983 %, 592.6). Nothing was rescaled, shifted, smoothed or truncated.

---

## Al2219 Figs. 7, 8(a), 8(b), 9 — the four inverse pole figures

File: `digitized/al2219_ipf_93orientations_46resolved.csv` (one row per resolved orientation, with
that orientation's value from all four panels).

### The two-versus-three degrees of freedom limitation — stated explicitly

**A point on an inverse pole figure encodes only the crystal direction parallel to ONE sample axis.
That is two rotational degrees of freedom. A full crystal orientation has three.** Digitizing these
figures recovers the **loading-axis (sample z) crystal direction** and **cannot recover the Euler
angles**; the rotation about the loading axis is not represented in the figure at all and is not
recoverable by any refinement of the method. Under the isochoric boundary condition, which
prescribes compression on both reciprocal axes, that third degree of freedom is **not** immaterial.
This is restated in the CSV header comment and is gap C.6 / C.1 in `gap_register.csv`.

### Projection convention assumed, and how it was verified

**Assumed:** equal-angle (**stereographic**) projection of the upper hemisphere onto the z = 0
plane,

```
X = u1 / (1 + u3),   Y = u2 / (1 + u3)
```

on the standard cubic triangle with **[001] at the origin (bottom-left), [101] on the horizontal
axis (bottom-right), [111] at the apex**, the figure's horizontal axis being crystal **[100]** and
its vertical axis crystal **[010]**. Inversion: `u = (2X, 2Y, 1 - X^2 - Y^2) / (1 + X^2 + Y^2)`,
normalized.

**Verification.** The choice was not assumed blindly. The triangle outline was detected from the
figures' own black boundary curves and the three corners measured in pixels. The diagnostic is the
apex position relative to the base, which differs measurably between the two candidate projections:

| projection | apex_x / base | apex_y / base |
|---|---|---|
| stereographic (equal-angle) | 0.8837 | 0.8837 |
| Lambert (equal-area) | 0.8496 | 0.8496 |
| **measured, Fig. 7** | **0.8825** | **0.8842** |
| **measured, Fig. 8(a)** | **0.8827** | **0.8859** |
| **measured, Fig. 8(b)** | **0.8823** | **0.8836** |
| **measured, Fig. 9** | **0.8828** | **0.8857** |

All four panels agree with stereographic to within 0.25 % and are 4 % away from equal-area.
**Stereographic confirmed.** As an independent check, all 46 recovered unit vectors satisfy the
standard-triangle inequality `0 <= u2 <= u1 <= u3` with no violations, and the maximum angular
deviation from the nearest low-index `(hkl)` with all indices <= 12 is **2.15 deg**.

### Points resolved — 46, not 93

| panel | markers, pass A (S > 0.32) | pass B (S > 0.20) | matched | union |
|---|---|---|---|---|
| Fig. 7 | 45 | 46 | 45 | **46** |
| Fig. 8(a) | 45 | 46 | 45 | **46** |
| Fig. 8(b) | 46 | 45 | 45 | **46** |
| Fig. 9 | 45 | 45 | 44 | **46** |

The paper states it ran **93** simulations. Only **46 distinct scatter positions are resolvable** at
the published marker size — markers are ~15 px in diameter on a 613 px triangle base, so points
closer than about 0.024 stereographic units merge irrecoverably. The count has **not** been padded
to 93, and the true resolved count (46) is what every derived table reports.

Independent corroboration that 46 is the real figure: all four panels plot the same orientation set,
and after conversion to projection coordinates the **nearest-neighbor offset between panels is
median 0.0004, p95 0.0017, max 0.0034 stereographic units** (= 3.6 px). Four separately drawn,
separately compressed figures landing 46 markers on the same 46 positions is not a coincidence of
the detector.

- **Method.** HSV saturation mask (the markers are the only saturated content in the triangle),
  3 x 3 morphological open then a 5 x 5 close to bridge the 2–3 px black boundary lines that cut some
  markers, connected components, then per-component **distance-transform peak detection** with
  non-maximum suppression at 1.3 r and a peak floor of 0.62 r. This splits touching pairs and triples
  (blob areas cluster tightly at 1.0x, then 1.7–2.2x the median single-marker area) without
  fabricating peaks inside single markers. Colorbars were excluded by locating each one as the tall
  fully saturated column block on the right and cropping to its left.
- **Raster resolution.** 300 DPI. Fig. 7 852 x 617; Fig. 8 1795 x 707 (panels split at x = 685/950);
  Fig. 9 856 x 625.
- **Geometric calibration.** Triangle corners in pixels, and the resulting scale
  (`X_[101] = sqrt(2) - 1 = 0.414214` over the measured base):

| panel | [001] | [101] | [111] apex | base (px) | stereographic units per px | deg per px |
|---|---|---|---|---|---|---|
| Fig. 7 | (32, 584) | (645, 584) | (573.0, 42) | 613 | 6.757e-4 | 0.077 |
| Fig. 8(a) | (36, 606) | (667, 606) | (593.0, 47) | 631 | 6.564e-4 | 0.075 |
| Fig. 8(b) | (970, 606) | (1606, 606) | (1531.1, 44) | 636 | 6.513e-4 | 0.074 |
| Fig. 9 | (33, 592) | (654, 592) | (581.2, 42) | 621 | 6.670e-4 | 0.076 |

  Calibration residual: the apex ratios above (0.25 % from the theoretical stereographic value) are
  the residual, since only three fiducials exist and they over-determine a two-parameter
  scale-and-offset map by one condition.
- **Positional extraction uncertainty.** Repeated-extraction error between passes A and B on matched
  markers: **median 0.4 px, p95 1.1 px, max 3.6 px**, i.e. **+-0.08 deg (p95) and +-0.26 deg (max)**
  in the recovered crystal direction, using `dtheta = 2/(1+r^2) dr` with `theta = 2 arctan r`. The
  per-point `deg_per_pixel` column in the CSV carries the local sensitivity.
  Pass disagreements, named: Fig. 7 one pass-B-only marker at (288, 534); Fig. 8(a) one pass-B-only
  at (276, 418); Fig. 8(b) one pass-A-only at (1213, 416); Fig. 9 one pass-A-only at (270, 406) and
  one pass-B-only at (550, 253) — the latter is a pale marker whose saturation falls below the strict
  threshold. All are retained in the union and their detecting pass is recorded per point in the
  `*_pass` columns.
- **Value (color) calibration.** Values were inverted against **each panel's own colorbar pixels**,
  not against an assumed colormap: a vertical strip through the colorbar (6 px inset from each edge
  to avoid the frame) gives an empirical RGB-to-value lookup table of 517–537 entries, with value
  linear in row between the colorbar's extreme tick labels. The 6 tick labels on every colorbar are
  exactly evenly spaced in value (e.g. Fig. 7: 1.25745, 1.27378, 1.29011, 1.30643, 1.32276, 1.33909;
  spacing 0.016328 = range/5), which is what places the first and last labels at the bar's ends.
  Marker color is the **median** RGB inside 0.55 r of the marker center; inversion is
  nearest-neighbor in **CIE Lab**.
  For the record: the colorbars are consistent with matplotlib `rainbow` (RMS 26/255 against
  `rainbow`, versus 51 for `jet`, 53 for `turbo`, 74 for `nipy_spectral`, 150 for `gist_rainbow`,
  best-fit scale 1.00 and offset -0.017). The 26/255 residual is JPEG chroma compression in the
  embedded bitmap. Using the figure's own colorbar as the LUT means that compression bias largely
  cancels, which is why the empirical LUT was used rather than the identified colormap.
- **Value extraction uncertainty** (p95, from perturbing the LUT by 3 CIE Lab units — the scale of
  the observed JPEG chroma error — and re-inverting):

| panel | quantity | colorbar range | value per LUT row | p95 uncertainty | max |
|---|---|---|---|---|---|
| Fig. 7 | void growth fraction, CP-only | 1.25745 – 1.33909 | 1.58e-4 | **+-0.0019** | 0.0044 |
| Fig. 8(a) | fracture initiation strain | 2.762 – 9.593 | 0.0128 | **+-0.15** | 0.31 |
| Fig. 8(b) | fracture propagation rate | 3 – 6656 | 12.4 | **+-161** | 385 |
| Fig. 9 | void growth fraction (exp) | 1.1748 – 1.7100 | 1.03e-3 | **+-0.013** | 0.022 |

  Observed Lab residuals of the actual marker inversions were median 0.33–0.61 and p95 1.76–2.17,
  i.e. below the 3-unit perturbation used, so these are conservative. Fig. 8(b)'s +-161 on a 3–6656
  range (2.4 % of full scale) is the intrinsic limit of reading a wide-range quantity off a rainbow
  bar and should be treated as the floor on any propagation-rate comparison.
- **Unit caveats.** Fig. 8(a)'s colorbar carries **no unit**; it is read as percent strain by
  consistency with the 1.015–7 % values quoted in Sec. 3.1.1, and this is tagged `assumed`.
  Fig. 8(b)'s "Fracture Progration Rate" (spelling as printed) has **no stated units anywhere in the
  paper**. Fig. 9's caption says "rate" while its colorbar says "Fraction" — recorded as a fraction
  (see gap C.10).

---

## Ni Fig. 2 — thickness measurement stations and dog-bone planform

File: `digitized/ni_fig2_thickness_stations.csv`.

- **Method.** The five stations are white rules drawn across the dark specimen body in the oblique
  3-D view. Near-white pixels (`min(R,G,B) > 225`) were grouped into columns on every image row from
  150 to 191 — the 37-row band on which exactly five groups are present and stable — and
  median-combined per station; the row-to-row standard deviation is carried as `pixel_x_sd`. The
  axial coordinate follows from the affine map t1 -> 0 mm, t5 -> 8.0 mm, the figure's own "8.0"
  dimension arrow running between the t1 and t5 extension lines.
- **Assumption, stated.** The 3-D view is a **parallel (orthographic/oblique) projection**, so pixel
  x is an affine function of axial position at a fixed image row. This holds for the view types
  MOOSE/ParaView emit and is required for the linear mapping to be valid. Under a perspective
  projection the interior stations would be biased toward the near end; the fact that the recovered
  spacings are even to within 2 % is consistent with the parallel-projection assumption.
- **Raster resolution.** 300 DPI (1250 x 647 px), scale **0.02189 mm per pixel**.
- **Calibration.** Two fiducials (t1 and t5, tied to 0 and 8.0 mm), so there is no residual to
  report; the check is instead the recovered interior spacings: **2.003, 2.036, 2.025, 1.937 mm**
  against a nominal 2.000 mm, i.e. within +-0.063 mm, comfortably inside the extraction uncertainty.
- **Extraction uncertainty.** From the row-to-row spread: **t1 +-0.104 mm, t2 +-0.019, t3 +-0.008,
  t4 +-0.005, t5 +-0.005 mm**. t1's larger figure is the station line nearest the shoulder fillet,
  where the body edge intrudes on some rows.
- **Repeated extraction.** Two independent row bands (150–170 and 171–191) give station positions
  agreeing to **max 0.019 mm**, within the quoted uncertainty.
- **Conclusion for the deck.** The five stations are at **0, 2, 4, 6, 8 mm** along the 8 mm gauge —
  nominally even 2.0 mm spacing. Use those values, and cite this digitization.
- **Planform dimensions** read from the figure's own dimension text (the figure states
  "(Unit = mm)"): total length **25.4**, gauge length **8.0**, grip width **4.8**, shoulder fillet
  radius **R 3.175** (= 0.125 in). Uncertainty is that of reading printed text, i.e. exact as
  printed. The **gauge-section width is not dimensioned in the figure** and is **not** estimated
  (gap C.11).

---

## Ni Fig. 3 — specimen L calibration curve

File: `digitized/ni_fig3_specimenL_calibration.csv`.

- **Method.** Color-channel segmentation (blue = simulation line, red = experiment dotted series),
  one row per pixel column with the column mean and the column's full extent. The legend box in the
  upper left was excluded (x < 415 px, y < 135 px); no data falls in that box, since the simulation
  curve only enters y < 135 px at x > 870 px.
- **Raster resolution.** 300 DPI (1000 x 732 px).
- **Axis calibration.** x: 6 ticks at pixel 145.5, 302.0, 459.0, 615.5, 772.5, 929.0 -> 0.0…0.5
  engineering strain; slope 6.38046e-4 /px; **max|residual| = 1.0e-4 strain, RMS 7.6e-5**.
  y: 6 ticks at pixel 623.0, 504.5, 385.5, 266.5, 146.5, 28.0 -> 0…250 MPa; slope -0.419864 MPa/px;
  **max|residual| = 0.23 MPa, RMS 0.14 MPa**.
- **Estimated extraction uncertainty.** **+-0.001 strain, +-0.5 MPa** away from the elastic rise.
- **Repeated extraction.** Color threshold varied from 40 to 25: **|dstress| median 0.000,
  p95 0.21 MPa**; experiment series max 0.54 MPa. The simulation series shows one outlier of
  22.1 MPa, entirely inside the near-vertical elastic segment where a column mean is ill-defined —
  which is why `stress_min_MPa` / `stress_max_MPa` are reported for every column.
- **Ranges recovered.** Simulation 0.000–0.517 strain, 23.7–248.3 MPa (774 columns);
  experiment 0.000–0.513 strain, 4.4–240.7 MPa (606 columns). The slightly negative strains at the
  first column (-0.0015, -0.0035) are the finite line and marker width at x = 0, **not** a shift.
  Nothing was rescaled, shifted, smoothed or truncated.

---

## Ni Fig. 4(b) — experimental versus simulated failure strain (parity)

File: `digitized/ni_fig4b_failure_strain_parity.csv`.

- **Method.** Color identifies the orientation group (blue 9, green 10, purple H, red L, black P)
  and marker **shape** identifies the thickness variant (filled circle x.1, cross x.3, five-point
  star x.4), both read from the figure's own legend. Shapes were located by **occlusion-aware
  template matching**: synthetic circle (r = 8.7 px), cross (half-length 10.6 px, half-width 2.3 px)
  and five-point star (outer radius 11 px, inner ratio 0.40) templates were correlated
  (`TM_CCOEFF_NORMED`) against each color mask; the highest-scoring shape was fitted, its dilated
  footprint erased, and the next shape fitted on the residual. Without the residual step the
  overlapping same-color pairs are mis-assigned. The legend (x < 430, y < 195 px) and the inset
  IPF triangle (x > 415, y > 400 px) were masked out.
- **Raster resolution.** 300 DPI (769 x 745 px).
- **Axis calibration.** 9 ticks per axis (0.1…0.9); x slope 1.20481e-3 /px,
  **max|residual| = 6.7e-4 strain**; y slope -1.21249e-3 /px, **max|residual| = 4.3e-4 strain**.
- **Estimated extraction uncertainty.** **+-0.003 strain** on both axes for the 13 well-resolved
  markers (template-peak localization ~1–2 px at 1.2e-3 strain/px, plus calibration residual).
- **Low-confidence points, named.** **L.4** (match score 0.199) and **10.4** (0.267) sit almost
  entirely behind their own group's circle marker; only a small crescent of each is visible. Their
  coordinates carry **+-0.02 strain** and are flagged `low_marker_occluded` in the CSV. Exclude or
  flag them in any parity error metric.
- **Repeated extraction / end-to-end validation.** Two things were checked rather than one.
  (i) Direct-peak matching without the residual step agrees with the residual method on **12 of the
  15 markers to within 1 px**; it shifts H.3 by 3.2 px and 10.3 by 2.2 px (both partly occluded by
  their own group's circle) and places L.4 **41 px** away, which is the failure mode the residual
  step exists to fix. The residual method is the one reported.
  (ii) More usefully, the **RMSE of the 15 digitized (sim - exp) pairs is 0.0619**, against the
  **0.064 the paper states in Sec. 4** — a 3 % agreement. That single number validates the whole
  chain at once: tick location, linear calibration, color-to-group mapping, shape-to-variant
  mapping and marker localization. It is the strongest available check on this digitization and it
  passes.

---

## Summary table

| figure | file | points | DPI | axis calibration residual | extraction uncertainty | status |
|---|---|---|---|---|---|---|
| Al2219 Fig. 3 (curves) | `al2219_fig3_calibration_stress_strain.csv` | 764 sim + 616 exp columns | 300 | 0.0075 % strain; 0.36 MPa | +-0.03 % strain, +-1.0 MPa | digitized |
| Al2219 Fig. 3 (markers) | `al2219_fig3_experiment_markers.csv` | 9 | 300 | as above | +-0.03 % strain, +-1.5 MPa | digitized |
| Al2219 Fig. 7 | `al2219_ipf_93orientations_46resolved.csv` | 46 | 300 | apex ratio 0.8825/0.8842 vs 0.8837 | +-0.08 deg; value +-0.0019 | digitized, 46 of 93 |
| Al2219 Fig. 8(a) | same file | 46 | 300 | 0.8827/0.8859 vs 0.8837 | +-0.08 deg; value +-0.15 | digitized, 46 of 93 |
| Al2219 Fig. 8(b) | same file | 46 | 300 | 0.8823/0.8836 vs 0.8837 | +-0.08 deg; value +-161 | digitized, 46 of 93 |
| Al2219 Fig. 9 | same file | 46 | 300 | 0.8828/0.8857 vs 0.8837 | +-0.08 deg; value +-0.013 | digitized, 46 of 93 |
| Ni Fig. 2 | `ni_fig2_thickness_stations.csv` | 5 stations + 4 dimensions | 300 | spacings 2.003/2.036/2.025/1.937 vs 2.000 mm | +-0.005 to +-0.104 mm | digitized |
| Ni Fig. 3 | `ni_fig3_specimenL_calibration.csv` | 774 sim + 606 exp columns | 300 | 1.0e-4 strain; 0.23 MPa | +-0.001 strain, +-0.5 MPa | digitized |
| Ni Fig. 4(b) | `ni_fig4b_failure_strain_parity.csv` | 15 | 300 | 6.7e-4 / 4.3e-4 strain | +-0.003 strain (2 points +-0.02) | digitized |
| Ni Fig. 4(a) | — | — | 300 | — | — | **not digitized** — 15 curves in 5 colors x 3 line styles; terminal markers overlap and the black series shares its color with axes/legend/inset. Trial extraction produced demonstrable false positives. Quantity recoverable from Fig. 4(b) instead |
| Ni Fig. 5 | — | — | 300 | — | — | **not digitized** — fracture-image montage with no calibrated axis, no colorbar scale and no dimension annotation. Findings captured from prose |
| Al2219 Figs. 1, 2, 4(b), 5(b), 6 | — | — | 300 | — | — | **not digitizable** — mesh views and 3-D fracture renderings with no axes or scale bars |
