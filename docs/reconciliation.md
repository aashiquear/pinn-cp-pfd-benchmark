
# Reconciliation — independent extraction vs. companion data sheet

**Sources of record:** Rezwan, Aragon, Montes de Oca Zapiain, Lim, *Int. J. Plasticity* **191** (2025)
104372 (Al2219; cited below as `rezwan-2025-intragranular`) and Rezwan, Fitzgerald, Ruggles,
Gilliland, Carroll, Lim, *Comput. Mater. Sci.* **272** (2026) 114872 (single-crystal Ni;
`rezwan-2026-predicting`).

**Companion sheet status:** navigation aid only. Everything below was extracted from the two papers
first, without reference to the sheet, and then compared to it. Where the two disagree, the paper
governs. Where the paper itself is internally inconsistent, that is recorded as a source defect
rather than resolved silently.

**Headline:** the sheet is substantially correct. Every numeric value it reports from Tables 1–3 of
both papers matches this extraction exactly. The discrepancies are (a) three symbol errors, one of
which is consequential by six orders of magnitude, (b) four places where the sheet reports an
inference as a source value, (c) five source defects and ambiguities the sheet does not flag, and
(d) one case where the sheet's fallback plan is only 49 % achievable.

---

## 1. Arithmetic sanity check, recomputed independently

The sheet (A.2) asserts that the athermal Taylor term is near 43.6 MPa and that
`tau_c = 0.112 MPa` is a separate additive lattice-friction term. **Both are confirmed.**

Recomputed from this extraction's own Table 2 values, in the paper's own um–MPa system:

```
sqrt(rho0_SSD)          = sqrt(220.0 um^-2)                 = 14.832397 um^-1
b * sqrt(rho0_SSD)      = 0.000286 um * 14.832397 um^-1      = 0.00424207   (dimensionless)
mu * b * sqrt(rho0_SSD) = 25 700 MPa * 0.00424207            = 109.02108 MPa
alpha0 * mu * b * sqrt(rho0_SSD)
                        = 0.4 * 109.02108 MPa               =  43.60843 MPa
```

Equation (10) of the source is `tau_c^a = tau_c^0 + alpha0*mu*b*sqrt(rho)`, i.e. strictly additive.
Therefore

```
tau_c^a(t=0) = 0.112 + 43.60843 = 43.72043 MPa
Taylor share = 43.60843 / 43.72043 = 99.74 %
Peierls share = 0.112 / 43.72043 =  0.26 %
```

**Confirmed:** `tau_c = 0.112 MPa` is a small additive lattice-friction (Peierls) term, not the
initial slip resistance. A reimplementation that loads 0.112 MPa as the total initial CRSS yields a
yield stress ~400x too low. Note also that Table 2 *labels* this quantity "Critical resolved shear
stress, tau_c", while Eq. (10) identifies it as the constant Peierls stress `tau_c^0` — the loose
table label is the origin of the trap the sheet warns about.

**One addition the sheet does not make.** Equation (11) defines the density that enters Eq. (10) as
the latent-hardened sum `rho = sum_b chi^ab ( rho_SSD^b + |rho_eGND^b| + |rho_sGND^b| )`. If that sum
is evaluated literally over 12 FCC systems at the initial state with `chi = 1` on the diagonal and
`q = 1.4` off it, then

```
rho_eff = rho0_SSD * (1 + 11q) = 220.0 * 16.4 = 3608 um^-2
alpha0*mu*b*sqrt(rho_eff) = 176.60 MPa            (a factor sqrt(1+11q) = 4.050 larger)
```

So the 43.6 MPa figure corresponds to `rho = rho0_SSD` — the **self** term only, not the
latent-hardened sum. Which of the two the original implementation used is not stated in the paper.
**This must be resolved before any calibration run**, because it moves the initial CRSS between
43.7 and 176.7 MPa. Recorded as gap C.15.

For comparison, the same computation for Ni (mm–MPa system):
`0.4 * 26 000 MPa * 2.5e-7 mm * 1532.97 mm^-1 = 3.9857 MPa`, so `tau_c^m(t=0) = 4.2457 MPa` and the
Peierls share is 6.1 % — a materially different balance from Al2219's 0.26 %.

---

## 2. Discrepancy table

| # | Quantity | My value (from the paper) | Sheet value | Assessment | Resolution |
|---|---|---|---|---|---|
| 1 | Al2219 Allen–Cahn kinetic parameter | Table 3 gives **`eta` = 1e-3 s/um**; Eq. (31) defines `1/eta = L`, so the mobility is `L = 1e3 um/s` | "Kinetic coefficient, **L** = 1e-3 s/um" (A.3) | **Sheet symbol error, consequential.** The tabulated number is `eta`, not `L`; they are reciprocals | Load 1e-3 s/um as `eta`. If PINNACLE's damage kernel takes a mobility `L`, pass 1e3 um/s. Entering 1e-3 as `L` slows damage kinetics by 1e6 and the crack will not propagate |
| 2 | Al2219 SSD storage coefficient symbol | Al2219 Table 2 symbol is **`k_m`** (= 0.3432); the Ni Table 3 symbol for the same quantity is **`k_0`** (= 0.0295) while the Ni Eq. (10) writes `k_m` | "`k0` = 0.3432" (A.2) for Al2219 | Value correct; **symbol collision originates in the sources**, and the sheet propagates the Ni table's symbol onto the Al2219 paper | Use `k_m` as the equation symbol for both materials and record `k_0` as the Ni-table alias. Do not let the two papers' symbol choice imply two different quantities |
| 3 | Plastic-work prefactor | Both papers use **`zeta`** = 1.0 (Al Eq. 31, Ni Eq. 28) | "beta = 1.0" (B.2) | Sheet symbol error, harmless | Use `zeta` |
| 4 | Al2219 void position | **Not stated.** Sec. 2.4 says only "a cubic domain measuring 50 um, containing a spherical void with a radius of 5 um"; Fig. 1(b) is an uncoordinated cut view | "single spherical, intracrystalline, radius 5 um, **centered**" (A.4) | **Sheet reports an inference as a source value** | Tagged `unavailable` in `benchmark_dataset.csv`. Centring is the only symmetric choice consistent with the reported fracture patterns and is what this study should assume — but it must be declared as an assumption, not inherited. New gap C.8 |
| 5 | Al2219 slip family | Paper states FCC and refers throughout to {111} slip, but **never writes "{111}<110>, 12 systems"**. The Ni paper does state it explicitly | "Crystal structure: FCC. Slip: {111}<110>, 12 systems." tagged **[T]** (A.1) | Sheet's reading is almost certainly right but is **not** table-sourced for Al2219 | Tagged `assumed` for Al2219, `reported_text` for Ni |
| 6 | Al2219 latent-hardening matrix form | Al2219 gives `q = 1.4` and calls `chi^ab` "the latent hardening matrix" but **never writes the matrix**. The Ni paper Eq. (9) gives `chi = 1 (m=m')`, `q (m!=m')` | Not addressed | Sheet gap | Take the Ni form for both; tagged `assumed` for Al2219 |
| 7 | Al2219 damage-variable convention | **The source contradicts itself.** Sec. 2.3 first paragraph: the scalar "transitions between 1 (intact material) and 0 (fully broken material)". Sec. 2.3 second paragraph and every equation (degradation `[(1-c)^2 + k]`, Eq. 31): 0 = intact, 1 = damaged | Not addressed | **Source defect the sheet does not flag** | Use **0 = intact, 1 = fully damaged** — the equation-consistent convention, and the one the Ni paper states cleanly |
| 8 | SSD annihilation term | Al2219 Eq. (12): `rho_dot_SSD = (k_m*sqrt(rho) - 2*y_c*rho_SSD) * gamma_dot / b`. Ni Eq. (10) prints `- (2/y_c)*rho_SSD` instead | Not addressed | **The Ni form is dimensionally impossible.** With `y_c` a length, `(2/y_c)*rho_SSD` has units L^-3 while `k_m*sqrt(rho)` has L^-1. The Al2219 form gives L^-1 - L·L^-2 = L^-1, which divided by `b` and multiplied by `gamma_dot` yields L^-2 per unit time, as required | Treat Ni Eq. (10) as a typesetting error; **use the Al2219 form `2*y_c*rho_SSD` for both materials.** This is load-bearing for the Kocks–Mecking mapping (C.5): the annihilation coefficient is `2*y_c`, not `2/y_c` |
| 9 | Ni slip-rate law actually used | Ni Eq. (6) prints a **thresholded** power law `gamma_dot = gamma_dot_0 <(|tau|-tau_c)/tau_c>^(1/m_s) sign(tau)`, but the same section then states that the regularization had minimal effect and "the original power-law slip relation is retained in this study" | Not addressed (the sheet records no slip law) | **Source ambiguity the sheet does not flag.** The printed equation is not the equation used | Implement the **unthresholded** Al2219 Eq. (9) form `gamma_dot = gamma_dot_0 (|tau|/tau_c)^(1/m_s) sign(tau)` for both materials, and keep the thresholded form available as a sensitivity variant |
| 10 | Ni Table 3 provenance | Table 3 is captioned "taken from Ref. [1]" (Lim 2025 quantitative), but Sec. 3.5 says only the **elastic constants and lattice parameters** were adopted and "the remaining plasticity parameters were tuned in this work" by manual trial-and-error against specimen L | The sheet records the Sec. 3.5 statement correctly but lists the whole of Table 3 as [T] | **Source defect**: caption and body text disagree | Elastic constants → `reported_table` (from Lim 2025); the ten plasticity parameters → `calibrated`. Tagged accordingly |
| 11 | Ni parity numbers | The paper **states in the text** that the parity RMSE is **0.064** across all 15 specimens (Sec. 4), and defines the failure strain as the axial strain at which `c > 0.85` across the cross-section | "Numerical parity values: **[F]** (paper Fig. and dataset)" (B.9) | **Sheet under-reports what is recoverable from the text.** The aggregate error metric and the metric definition are both in prose; only the 15 individual coordinates are figure-only | RMSE 0.064 tagged `reported_text`. The 15 individual (exp, sim) pairs digitized from Fig. 4(b) — and their own RMSE recomputes to **0.0619**, 3 % from the reported 0.064, which validates the digitization chain end to end |
| 12 | Failure-metric thresholds | Al2219: fracture **initiation** at `c = 0.65` anywhere in the domain (Sec. 3.2). Ni: **failure strain** at `c > 0.85` **across the cross-section** (Sec. 4) | A.9 records 0.65; the 0.85 threshold is not recorded anywhere | **Sheet gap with a comparison trap.** One is a local initiation criterion, the other a through-thickness separation criterion | Both recorded. Any benchmark table reporting "failure strain" for both materials **must** state which definition applies per material; the two numbers are not interchangeable |
| 13 | Ni dog-bone planform | Fig. 2 prints, as dimension text with an explicit "(Unit = mm)": total length **25.4**, gauge length **8.0**, grip width **4.8**, shoulder fillet radius **R 3.175** (= 0.125 in). The **gauge-section width is not dimensioned** | "Dog-bone planform dimensions other than length and gage length: **[F]** (Fig. 2)" (B.5) | Sheet correct that they are figure-only; **two of the three are now recovered** | Grip width and fillet radius tagged `digitized`. Gauge-section width tagged `unavailable` — it cannot be recovered from the oblique 3-D view without assuming fillet tangency, and is **not** estimated here. New gap C.11 |
| 14 | Ni thickness station locations | Digitized from Fig. 2: **0, 2.00, 4.04, 6.06, 8.00 mm** from t1, uncertainty +-0.10 mm (t1) to +-0.005 mm (t5) | "[F] ... not given as coordinates in the text" (B.4) | Sheet correct; **now recovered** | Consistent with nominally **even 2.0 mm spacing** across the 8 mm gauge to within the extraction uncertainty. Use 0/2/4/6/8 mm and cite the digitization |
| 15 | The 93 orientations | Digitization of Figs. 7, 8(a), 8(b) and 9 resolves **46 distinct scatter positions**, identical in all four panels (max cross-panel offset 0.0034 stereographic units = 3.6 px). The other 47 are not separable at the published marker size (~15 px diameter on a 613 px triangle base) | A.8 presents IPF digitization as the fallback route to the 93 orientations without stating a recoverable count | **Sheet's fallback is only 49 % achievable, and recovers the wrong quantity anyway** | 46 loading-axis directions delivered as Miller indices and unit vectors, **not** padded to 93. The sheet's own critical limitation is confirmed and restated: an IPF point is 2 rotational DOF, a full orientation is 3, so **no** amount of digitization recovers the Euler angles. Refines gap C.1 |
| 16 | Fig. 9 quantity name | Caption: "Void growth rate observed during the experiment". The figure's own colorbar: "**Void Growth Fraction (Exp)**", range 1.1748–1.7100 (dimensionless) | "experimentally observed void growth **rate**" (A.8) | **Source defect**: caption and axis label disagree. The sheet follows the caption | Recorded as a **fraction** (dimensionless), because that is what the axis says and because the numeric range 1.17–1.71 is a ratio, not a rate. New gap C.10 |
| 17 | Fig. 8(a) units | Colorbar reads 2.762–9.593 with **no unit**. Sec. 3.1.1 quotes fracture strains of 1.015–7 % for the extreme orientations | Not addressed | Sheet gap | Read as **percent strain**, tagged `assumed` in the digitized CSV and in `benchmark_dataset.csv`. Fig. 8(b) "Fracture Progration Rate" (sic) has **no stated units at all** — recorded as unstated |
| 18 | Al2219 initial time step | Both papers: 0.045. **Neither states the time unit** | "Initial time step: 0.045" (A.6), no caveat | Sheet omits the caveat | Tagged with the caveat. New gap C.14 |
| 19 | Ni fracture angle | "approximately 45 degrees" (Abstract, Sec. 4). **The paper never states whether the angle is measured from the loading axis or from the transverse/normal direction** | B.9 records "approximately 45 degrees" without the convention caveat | Sheet gap | 45 deg recorded with an explicit note that the convention is undefined. Fix and document one convention before any angle comparison. New gap C.9 |
| 20 | Ni specimen 9.2 | The three tested group-9 specimens are **9.1, 9.3, 9.4** (Table 2). The paper's phrase "specimens (e.g., 9.1–9.3 and P.3–P.4)" therefore denotes 9.1 and 9.3 | B.9 quotes "9.1-9.3" verbatim | Not an error, but the range notation is misleading | Read as {9.1, 9.3} and {P.3, P.4} |
| 21 | Gc across the two materials | Al2219: 22.75 MPa-um = **22.75 J/m^2**. Ni: **125 kJ/m^2** = 1.25e5 J/m^2 = 125 MPa-mm. Ratio **5495** | Both values reported correctly and separately; the magnitude gap is not flagged | Sheet gap | Both are **effective calibration parameters**, not measured toughnesses — 22.75 J/m^2 is a brittle-ceramic-scale value for a ductile Al alloy. The gap is not necessarily an error (the phase-field driving force scales with `Gc/l`, and the Ni `l_0` is unreported), but **no Gc may be transferred between the two decks**, and any cross-material narrative must state this. New gap C.16 |
| 22 | Digitizable Ni figures | Fig. 2 (geometry) and Fig. 3 (calibration curve) digitized. Fig. 4(b) (parity) digitized. **Fig. 4(a) could not be reliably digitized** (below). **Fig. 5 carries no numerically digitizable axis** — it is a montage of experimental and simulated fracture images over a strain color field | C.7 lists "Figs. 3, 4, 5 of the Ni paper" as the digitization set | Sheet over-scopes | See "Figures not digitized" below |
| 23 | Al2219 polycrystal calibration domain | 50 x 50 x 100 um, **30 grains**, random orientations from a 3-D phase-field grain-growth simulation; same element type and size as the single-crystal model. **The 30 orientations themselves are unavailable** | A.2 records the method but not the domain dimensions or grain count | Sheet gap | Recorded. Any reproduction of the Fig. 3 calibration must generate its own random set and report the seed. New gap C.13 |

### Values checked and found in exact agreement

Al2219 `C11/C12/C44` (108 / 61.3 / 28.5 GPa); `b` 0.000286 um; `gamma_dot_0` 1e-5 1/s; `m` 0.12;
`mu` 25.7 GPa; `tau_c` 0.112 MPa; `q` 1.4; `alpha0` 0.4; `rho0_SSD` 220.0 um^-2; `k_m` 0.3432;
`y_c` 0.00268 um; PFD `k` 1e-3, `l` 1 um, `Gc` 22.75 MPa-um; domain 50 um cube, void radius 5 um;
mesh density 1 / 0.65; HEX8 + linear Lagrange; both loading definitions; all six solver settings;
the four extreme orientations including the overbar reading of `[1bar 4 9]` as `(-1, 4, 9)`;
the 0.65 initiation threshold; the void-growth-fraction definition (12 % strain).

Ni `C11/C12/C44` (253 000 / 152 000 / 124 000 MPa); `b` 2.5e-7 mm; `gamma_dot_0` 1e-5 1/s;
`m_s` 0.12; `mu` 26 000 MPa; `tau_c` 0.26 MPa; `q` 1.4; `alpha0` 0.4; `rho0_SSD` 2.35e6 mm^-2;
`k_0` 0.0295; `y_c` 0.79e-6 mm; `Gc` 125 kJ/m^2; all 15 x 5 = 75 thickness values of Table 2 to
four decimals; all 15 Bunge Euler triples and 5 Schmid factors of Table 1; 9880 HEX8 elements;
spline thickness interpolation; 25.4 mm length, 8 mm gauge, ~0.35 mm nominal thickness, 30–69 um
(2–18 %) variation; all boundary conditions; all solver settings; the processing route
(2.5 mm VAR 99.999 % Ni disk, 8 h at 1350 C in vacuum, EDM slicing into x.1–x.4, EBSD on x.2);
the DIC settings (6.66 Hz, subset 35 px, step 5 px, strain window 5).

The sheet's unit conversions in B.1 (`rho0_SSD` 2.35e6 mm^-2 = 2.35 um^-2; `b` 2.5e-7 mm =
0.00025 um; `y_c` 0.79e-6 mm = 0.00079 um) were recomputed and are all correct.

---

## Figures not digitized, and why

- **Ni Fig. 4(a)** — simulated engineering stress–strain curves for all 15 specimens, drawn as
  5 colors x 3 line styles with end-of-curve `x` markers. Color separates the orientation groups
  cleanly, but within a group the three line styles overlap heavily and the terminal `x` markers of
  different variants frequently coincide or fall on top of neighboring curves; the black (group P)
  series additionally shares its color with the axes, tick labels, legend glyphs and the inset
  triangle. Trial extraction returned spurious markers on curve crossings (e.g. a "group 10" point
  at 0.009 strain / 77 MPa, which is on the elastic knee, not a failure marker), so no numbers from
  this panel are reported. The quantity it encodes — simulated failure strain per specimen — is
  recovered instead from Fig. 4(b), which is unambiguous. Simulated *failure stress* per specimen is
  therefore **not** recovered; it is a simulation output the reference-FEM stage will regenerate.
- **Ni Fig. 5** — side-by-side experimental and simulated fracture-line montages over an
  `e_xx` strain field. There is no calibrated axis, no colorbar scale and no dimension annotation,
  so neither fracture location nor fracture angle can be measured from it to any stated accuracy.
  The qualitative findings it supports (45 deg inclined surfaces for groups 9 and P, distributed
  necking for groups 10 and H, failure at a local thickness minimum except for 9.4 and P.1) are
  captured from the prose instead.
- **Al2219 Figs. 1, 2, 4(b), 5(b), 6** — mesh views and 3-D fracture-surface renderings. No axes,
  no scale bars, no colorbars; nothing numerically digitizable. The strain values annotated in the
  Fig. 4 and Fig. 5 captions are captured as `reported_text`.
