# Ductile-like failure of single crystals in coupled crystal-plasticity phase-field models: literature review, assessment and test plan

Date: 2026-11-15. Scope: the fifteen single-crystal Ni specimens of (rezwan 2026 predicting)
as simulated with PINNACLE (`docs/reports/ni_single_crystal_pinnacle_benchmark_report.md`).
Question: which modifications of the common CP-PFD damage model are used for ductile-like
fracture, which of them could explain the early failure of orientation groups 9 and H, and
how to test them. Every reference below is in `bib/references.bib`; DOIs were verified
against Crossref on 2026-11-15 unless flagged. Section 9 lists what was not verified.

## 1. Summary

1. **A finding in our own setup comes before the literature.** PINNACLE applies the three
   Bunge angles in the opposite order to the paper. The maximum Schmid factors of the
   paper's Table 1 (0.44, 0.48, 0.44, 0.42, 0.46 for 9, 10, H, L, P) are reproduced to two
   decimals by the standard Bunge convention and not by the rotation PINNACLE builds
   (0.39, 0.48, 0.34, 0.40, 0.35). Passing the angles as ($\varphi_2$, $\Phi$, $\varphi_1$)
   reproduces the paper's values in PINNACLE. All Ni runs so far therefore simulated five
   different crystals from the tested ones. Groups H and P were run as hard, near-[111]
   crystals; in the paper H is a soft, near-[101] crystal. Section 2 gives the evidence.
   This must be corrected and the calibration repeated before any fracture-model change
   is judged.
2. **A constant prefactor on the plastic work cannot change the ranking of orientations.**
   With the elastic energy negligible, the phase-field equation depends on $\zeta / G_c$
   only. After $G_c$ is refitted on L.1, any constant $\zeta$ gives the same result for
   every specimen. The author's hypothesis is testable only in a form where the plastic
   contribution depends on the state: a threshold, a stored-energy fraction that varies
   with the dislocation density, or a slip-system-resolved measure (Section 6).
3. **The literature offers five families of modification** (Section 4): a plastic-work or
   plastic-strain threshold; degradation of the yield surface, not only of the elastic
   stiffness; a fracture toughness that degrades with plastic strain; stress-triaxiality
   dependence; and a stored-energy driving force normalized by the dislocation density.
   None of the CP-PFD papers found validates failure strain against a series of
   single-crystal tensile orientations. The benchmark paper is the only such data set
   found.
4. **Experiments on FCC single crystals** describe tensile failure as a plastic
   instability set by slip geometry: lattice rotation, geometric softening, necking, then
   localized shear on one or two slip systems (Section 5). The failure strain is set
   mainly by when the neck forms, which is a crystal-plasticity result, not a
   fracture-energy result.
5. **Recommended order of work:** correct the orientation convention; rerun crystal
   plasticity on L.1 and the $G_c$ fit; check 9.1 and H.1; only if they still fail early,
   test the threshold form and then yield degradation. Section 8 gives the plan.

## 2. The orientation convention

### 2.1 Evidence

The paper's Table 1 gives the Bunge angles and the maximum Schmid factor $M_s$ among the
twelve {111}<110> systems for the loading direction $x$. Table A compares them with the
two possible readings of the angles.

**Table A.** Maximum Schmid factor on the loading axis $x$.

| group | $\varphi_1$, $\Phi$, $\varphi_2$ (deg) | paper, Table 1 | standard Bunge | PINNACLE as run | PINNACLE with ($\varphi_2$, $\Phi$, $\varphi_1$) |
|---|---|---|---|---|---|
| 9 | 291.2, 47.9, 49.3 | 0.44 | 0.436 | 0.393 | 0.436 |
| 10 | 126.0, 45.0, 56.5 | 0.48 | 0.478 | 0.475 | 0.478 |
| H | 189.8, 50.6, 40.7 | 0.44 | 0.444 | 0.342 | 0.444 |
| L | 60.0, 23.5, 80.2 | 0.42 | 0.417 | 0.402 | 0.417 |
| P | 335.1, 45.4, 77.0 | 0.46 | 0.460 | 0.352 | 0.460 |

"Standard Bunge" is $g = Z(\varphi_2)\,X(\Phi)\,Z(\varphi_1)$ with passive matrices, mapping
sample to crystal axes, so the loading direction in the crystal is $g\,\mathbf{e}_x$.
PINNACLE (`materials/schmid.py: bunge_rotation`, `materials/cp_field_callback.py`,
`materials/elasticity.py`) builds $R = R_z(\varphi_2) R_x(\Phi) R_z(\varphi_1)$ with active
matrices and uses it as the crystal-to-sample rotation. The crystal-to-sample rotation of
the standard convention is $g^T = R_z(\varphi_1) R_x(\Phi) R_z(\varphi_2)$. The two differ
by the exchange of $\varphi_1$ and $\varphi_2$. The last column was computed with
PINNACLE's own material object (Schmid tensors `P`, component $xx$).

Three further observations agree with this.

- The inverse pole figure of the paper's Fig. 4 puts 9 and H at the [101] corner, L toward
  [111], 10 on the [001]-[101] side and P in the interior. The standard convention gives
  loading directions of (0.71, 0.69, 0.13) for 9, (0.72, 0.68, 0.13) for H,
  (0.70, 0.63, 0.35) for L, (0.80, 0.57, 0.17) for 10 and (0.82, 0.49, 0.30) for P, as
  drawn. PINNACLE as run puts H at (0.68, 0.54, 0.50) and P at (0.69, 0.53, 0.49), both
  near [111].
- The peak stresses of our fifteen-specimen runs rank exactly by the as-run Schmid
  factors: H 247 MPa (0.342), P 234 (0.352), 9 224 (0.393), L 223 (0.402), 10 158 (0.475).
- The paper states that the more ductile groups (H) show lower flow stress and the less
  ductile group L the higher. Our runs give H the highest flow stress of all.

### 2.2 Consequences

- Groups 9 and H are, in the tested crystals, softer than L, with two nearly equal slip
  systems (Schmid 0.44 and 0.42 to 0.43). A plastic-work criterion then lets them strain
  further than L. Our runs gave them the same or a higher stress than L and therefore the
  same or a smaller failure strain. The reading of report §7.4, that one parameter set
  cannot serve both classes, rests on runs with the wrong crystals and is open again.
- The crystal-plasticity set was fitted to the curve of L with a crystal that is not L.
  The three fitted values moved to 0.86, 0.83 and 0.59 of the paper's values, and a strain
  rate of $10^{-3}$ s$^{-1}$ was assumed to raise the stress. Both steps must be repeated.
- The surrogate's Schmid features use the same rotation, so the surrogate is consistent
  with its training runs, but its fifteen "real specimens" are not the real crystals.
- Parts I and II (Al2219 polycrystal plate) use the same rotation. With thirty grains the
  effect on the averaged response is expected to be small, but the grain orientations are
  not the ones intended. This needs its own check.

## 3. What the present model does

The model is AT2 with the history field $\mathcal{H} = \max[\chi\psi^+ + \zeta\langle W_p - W_c\rangle]$,
$\chi = \zeta = 1$, $W_c = 0$, and $(G_c/\ell_0)\,d - G_c \ell_0 \nabla^2 d = 2(1-d)\mathcal{H}$.
The degradation $(1-d)^2$ multiplies the elastic second Piola-Kirchhoff stress.

**Damage is homogeneous from the first step.** Without a threshold the uniform solution is
$d = 2\mathcal{H} / (G_c/\ell_0 + 2\mathcal{H})$ with $G_c/\ell_0 = 82$ MPa. Table B gives
the largest $d$ in our finite-element runs before any localization.

**Table B.** Maximum phase field and engineering stress in the runs of report §7.

| gauge strain | L.1: $d$ | L.1: stress (MPa) | 9.1: $d$ | H.1: $d$ | 10.1: $d$ |
|---|---|---|---|---|---|
| 0.10 | 0.21 | 134 | 0.23 | 0.25 | 0.18 |
| 0.20 | 0.41 | 181 | 0.44 | 0.46 | 0.34 |
| 0.30 | 0.56 | 209 | 0.59 | 0.60 | 0.48 |
| 0.40 | 0.69 | 222 | 0.72 | 0.73 | 0.59 |

**The stress hardly feels it.** At 40% strain L.1 carries 222 MPa with $d = 0.69$ against
229 MPa without damage. The reason is that only the elastic stiffness is degraded. Slip
needs the resolved shear stress to reach the undegraded slip resistance, so the elastic
strain grows by $1/g$ and the flow stress stays. The phase field acts as a counter of
plastic work that has no mechanical effect until $g$ is a few percent. Failure is then, to
a good approximation, the statement that the local plastic work reaches a fixed multiple
of $G_c/\ell_0$.

**Two consequences.**

- With $\psi^+ \ll W_p$ (here about 0.3 against 80 MPa) the equation contains $\zeta$ and
  $G_c$ only as $\zeta / G_c$. A constant $\zeta$ is the same as a different $G_c$.
- Two orientations with the same stress-strain curve reach the critical work at the same
  strain. In our runs L.1 and 9.1 had the same curve and failed at 0.464 and 0.444. No
  scalar function of the plastic work can separate them. A model can separate them only
  through a different flow stress, a different localization, or a quantity that is not
  the plastic work.

A measure that is invariant to any such rescaling is the ratio of failure strains to L.1.
The experiment gives 1.18 for 9.1 and 1.19 for H.1; our runs gave 0.96 and 0.92.

## 4. Modifications used in the literature

**Table C.** Modifications of the phase-field ductile fracture model. CP = applied with
crystal plasticity. A dagger marks an equation given in its commonly summarized form that
was not read in the source (Section 9).

| # | Modification | Form | New parameters | Physical content | Key references | CP |
|---|---|---|---|---|---|---|
| 1 | Weight on the plastic work | $\mathcal{H} = \beta_e \psi^+ + \beta_p W_p$ † | $\beta_p$ (our $\zeta$) | Share of the plastic work that drives the crack | (borden 2016 triaxiality), (miehe 2016 gradient), (rezwan 2026 predicting) | yes |
| 2 | Plastic-work threshold | $\beta_p \langle W_p - W_0 \rangle$ † | $W_0$ | No plastic damage before a critical work | (borden 2016 triaxiality) | yes, PINNACLE `plastic_work_threshold` |
| 3 | Critical plastic strain in the driving force | driving state built from $\langle \varepsilon_p / \varepsilon_{crit} - 1 \rangle$ † | $\varepsilon_{crit}$ | Damage starts at a critical equivalent plastic strain | (miehe 2015 ductile) | no |
| 4 | Plastic-strain-dependent degradation | $g = (1-d)^{2p}$, $p = \varepsilon_p / \varepsilon_{p,crit}$ | $\varepsilon_{p,crit}$ | Elastic energy drives; plastic strain raises the sensitivity | (ambati 2015 ductile), (ambati 2016 verification) | yes: (delorenzis 2016 crystal) with accumulated slip |
| 5 | Elastic energy only | $\mathcal{H} = \psi^+$ | none | Plasticity acts through the stress field alone | (duda 2015 gradient) | yes: (shanthraj 2016 viscoplastic), (diehl 2017 void) |
| 6 | Yield-surface (plastic) degradation | flow resistance multiplied by $g_p(d)$ | form of $g_p$ | Damaged material loses flow strength; plastic softening and localization inside the crack | (borden 2016 triaxiality), (alessi 2018 comparison), (vakili 2022 aluminum) | yes: (vakili 2022 aluminum) degrades the stored hardening energy $\varphi^2 \sum_\alpha \int g^\alpha d\gamma^\alpha$ |
| 7 | Toughness degrading with plastic strain | $G_c(\bar\varepsilon_p)$, e.g. $\left[1-(1-\beta)(1-e^{-\bar\varepsilon_p/\varepsilon_0})\right] G_c$ | $\beta$, $\varepsilon_0$ | Void growth and coalescence lower the resistance | (dittmann 2018 variational), (yin 2020 toughness), (hu 2021 coalescence), (han 2022 toughness) | not found |
| 8 | Stress-triaxiality dependence | threshold or driving term scaled by triaxiality † | triaxiality parameters | Void growth depends on the mean stress | (borden 2016 triaxiality), (abrarivajari 2022 multiaxial), (sur 2025 hydrostatic) | not found |
| 9 | Porous plasticity coupling | Gurson-type yield; phase field driven by void fraction | void parameters | Explicit void nucleation and growth | (aldakheel 2018 gurson) | not found |
| 10 | Stored-energy density | $G_s = \int \xi\, \boldsymbol\sigma : d\boldsymbol\varepsilon^p / \sqrt{\rho_{SSD} + \rho_{GND}}$, threshold $G_{crit}$ | $\xi$, $G_{crit}$ | Only the stored share of the work, per unit of dislocation mean free path, breaks bonds | (wilson 2018 driving), (lucarini 2023 cracking) | yes, fatigue |
| 11 | Defect energy | stored elastic plus defect (hardening) energy drives the crack | none new | Energy stored in the dislocation structure | (tu 2020 series), (maloth 2023 microstructures) | yes |
| 12 | Anisotropic or plane-resolved fracture | cleavage-plane or slip-plane damage variables | plane toughness | Cleavage and slip-plane decohesion | (shanthraj 2017 cleavage), (li 2024 cleavage) | yes |
| 13 | Elastic stage and length-independent strength | AT1 or cohesive-type degradation (Lorentz, PF-CZM) | strength | No damage below a threshold; strength independent of $\ell_0$ | (pham 2011 gradient), (lorentz 2011 cohesive), (wu 2017 unified), (hu 2021 coalescence) | not found |
| 14 | Gradient-enhanced plasticity | plastic length scale next to $\ell_0$ | plastic length | Damage zone stays inside the plastic zone; delays localization | (miehe 2016 gradient), (auth 2024 thermodynamic), (shishvan 2021 mechanism) | yes |

Reviews: (alessi 2018 comparison) compares the coupling choices; (marengo 2023 review)
surveys the small-strain formulations.

**Notes on the entries that matter here.**

- Items 1 and 2 are the model we run. (borden 2016 triaxiality) keeps the plastic weight
  at one and uses the threshold $W_0$ to delay plastic damage; this was read in a
  secondary source, not in the paper.
- Item 6 is the main structural difference between our model and the ductile models of
  the literature. Without it the phase field has no plastic softening, as Table B shows.
- Item 10 is the closest published form of the author's idea of a plastic contributor.
  With Taylor hardening, $\tau_c \approx \alpha \mu b \sqrt{\rho}$, the rate of $G_s$ is
  close to $\xi \alpha \mu b \sum_\alpha |\dot\gamma^\alpha|$: the criterion becomes a
  critical accumulated slip, not a critical work. It no longer penalizes an orientation
  for carrying a high stress. It was developed and validated for fatigue crack
  nucleation and short crack growth, with $\xi$ a constant; its use for monotonic ductile
  failure is not established in the sources found.
- (vakili 2022 aluminum) states the absence of ductile phase-field models for crystal
  plasticity at large deformation as its motivation.

## 5. Crystal plasticity applications and experimental evidence

**CP-PFD models.** (delorenzis 2016 crystal) couples single-crystal plasticity to the
model of item 4 and finds fracture initiating at the maximum accumulated slip.
(hernandezpadilla 2017 crystal) gives a two-dimensional model for a Ni-base single
crystal. (diehl 2017 void) shows crack initiation around a void depending on orientation
and pore shape. (xue 2025 orientation) simulates orientation-dependent crack growth in
ductile single crystals. (auth 2024 thermodynamic) adds gradient-enhanced hardening for
transgranular cracks in single grains. None of these compares a predicted failure strain
with tensile tests of several orientations. (rezwan 2026 predicting) does, with
$\zeta = 1$ and one $G_c$, and reports an RMSE of 0.064 with the largest scatter in the
high-ductility cases.

**Experiments on FCC single crystals.**

- (simoto 1965 copper): in copper crystals oriented for multiple slip, necking starts by
  general deformation with lattice rotation in the neck, followed by localized shear
  with void growth after a well-developed neck has formed.
- (jasienski 1985 instabilite): copper, aluminum and nickel crystals in single and
  multiple slip orientations, with stress-strain curves followed beyond the necking
  point; the shear bands are related to the activation of coplanar slip systems.
- (asaro 1979 geometrical) and (peirce 1982 localized): localization in ductile single
  crystals follows from geometrical softening by lattice rotation; the finite-element
  analysis reproduces necking and shear bands with crystal plasticity alone, without a
  damage model.
- (zhang 2007 necking): necking of a square copper single-crystal specimen simulated
  with crystal plasticity, sensitive to a small misalignment of the loading axis.
- (kocks 2003 hardening): the hardening stages of FCC crystals depend on the number of
  active systems; the hardening rate enters the Considère condition for necking.
- The benchmark paper itself: groups near [101] (10, H) were the more ductile, groups
  toward [111] (L, P) the less; the maximum Schmid factor alone does not rank the
  ductility; the authors attribute the remaining error to the onset and evolution of the
  plastic instability.

**Reading.** In thin single-crystal sheet the strain to failure is mainly the strain at
which the neck forms plus a short post-necking stage. The first is set by hardening rate
and lattice rotation, which the crystal-plasticity model predicts if the orientation and
the hardening law are right. The fracture model should then act late and locally. Our
present model does the opposite: damage is spread over the gage from the start (Table B).

## 6. Assessment of the hypothesis

The hypothesis is that a plastic contributor parameter controls the ductile-like behavior.

| Form | Can it change the failure-strain ratio of two orientations? | Reason |
|---|---|---|
| Constant $\zeta$ | No | Only $\zeta/G_c$ enters; the refit of $G_c$ on L.1 cancels it |
| Constant stored fraction $\xi$ of the work | No | Same as constant $\zeta$ |
| Threshold $W_c$ with $G_c$ refitted | Yes, weakly | Shifts the criterion from total work toward work after localization; orientations that neck late gain |
| Stored energy per mean free path (item 10) | Yes | Turns critical work into critical slip; removes the penalty of a high flow stress |
| Slip-resolved or orientation-dependent $\zeta$ | Yes | By construction; no physical basis found in the literature, and it is a fit per orientation |
| Yield degradation (item 6) | Yes | Couples damage to localization; changes the post-necking stage |

So the hypothesis is right in spirit and wrong in its simplest form. A plain sweep of
$\zeta$ at fixed $G_c$ repeats the existing $G_c$ search of report Table 5. The published
form closest to it is the stored-energy density of item 10.

Whether any of this is needed is not known until Section 2 is corrected. The paper's own
model, with $\zeta = 1$ and the right crystals, predicts H late, not early.

## 7. Recommendation

1. **Correct the Euler convention and repeat the calibration.** It is an error of fact,
   it is cheap to fix, and it changes the flow stress of four of the five groups.
2. **If 9 and H are still early, make the damage act late:** a plastic-work threshold
   with $G_c$ refitted (available now), then an AT1 or cohesive-type degradation with an
   elastic stage (item 13). The measure of success is the pre-localization $d$ of Table B
   falling below about 0.1 and the peak stress of L.1 reaching the experimental 239 MPa.
3. **If the post-necking stage is then wrong, add yield degradation** (item 6). This is a
   PINNACLE change to the crystal-plasticity update and needs an ADR.
4. **The stored-energy density** (item 10) is the form in which to test the author's
   idea, as a fourth step, because it changes the criterion from work to slip. It needs
   a new driving-force option in PINNACLE.

Toughness degradation, triaxiality and porous plasticity (items 7 to 9) address void
growth in polycrystals at high triaxiality. In a 0.35 mm sheet the triaxiality is low and
these are not the first candidates.

## 8. Test plan

Every step uses one or three specimens; no campaign is started without the author's word.
The invariant measure throughout is the failure-strain ratio to L.1 (experiment: 9.1 at
1.18, H.1 at 1.19, 10.1 at 1.17, P.1 at 1.01).

| Step | Test | Runs | Decides |
|---|---|---|---|
| T0 | Fix the convention: PINNACLE key or corrected rotation with a test against the paper's Table 1 Schmid factors (ADR); until then decks pass ($\varphi_2$, $\Phi$, $\varphi_1$) | none | Framework or deck fix, author's choice |
| T1 | Material-point curves of the five orientations in both conventions | seconds | Flow-stress ranking; expected L highest, 10 lowest |
| T2 | L.1 without damage: the paper's Table 3 set and the carried set, at $10^{-4}$ and $10^{-3}$ s$^{-1}$ | 4 runs of about 25 min | Whether the paper's set and a lower rate now match the curve of L; refit if not |
| T3 | $G_c$ on L.1 | about 6 runs | New single set |
| T4 | 9.1 and H.1 with the new set | 2 runs | Ratios against 1.18 and 1.19. Within 10%: go to the fifteen specimens and stop here |
| T5 | Threshold: $W_c$ at 25, 50, 75% of the plastic work of L.1 at its peak, $G_c$ refitted each time; L.1, 9.1, H.1 | about 9 runs plus fits | Whether late damage changes the ratios; Table B and the peak stress as checks |
| T6 | Constant $\zeta$ = 0.5 and 0.1 with $G_c$ scaled by the same factor, on L.1 | 2 runs | Confirms numerically that constant $\zeta$ is a rescaling (expected: identical failure strain) |
| T7 | Yield degradation in PINNACLE (ADR, tests), then L.1, 9.1, H.1 | after T5 | Post-necking behavior and fracture angle |
| T8 | Stored-energy driving force $G_s$ (ADR, tests), $\xi$ constant, $G_{crit}$ fitted on L.1; L.1, 9.1, H.1 | after T7 | The author's hypothesis in its state-dependent form |

Open inputs that would sharpen the tests: the experimental stress-strain curves of groups
9, 10, H and P. The paper shows only that of L (its Fig. 3); its Fig. 4(a) is simulated.
With them, T2 could be checked on all five orientations. The Al2219 parts need a separate
check of the convention after T0.

## 9. Not verified

- The equations marked with a dagger in Table C were not read in the primary papers. The
  threshold form of (borden 2016 triaxiality) comes from a secondary summary; the
  triaxiality terms were not read at all.
- (maloth 2023 microstructures) and (tu 2020 series): the defect-energy driving force is
  taken from abstracts.
- (vakili 2022 aluminum) is a preprint without a Crossref DOI; its Eq. (17) was read.
- (jasienski 1985 instabilite), (simoto 1965 copper), (asaro 1979 geometrical),
  (peirce 1982 localized), (zhang 2007 necking), (kocks 2003 hardening): statements come
  from abstracts or general knowledge of the papers, not from a full reading.
- No source was found on the effect of a regularization length larger than the sheet
  thickness ($\ell_0 = 0.6$ mm against 0.35 mm). The mesh check of report §5 covers the
  element size only.
- The strain rate of the tests is not stated in the paper.
- The values read from the paper's Fig. 4 (ranking of simulated flow stresses, position
  of the groups in the inverse pole figure) are qualitative readings.

## 10. References

Keys of `bib/references.bib`:

ambati-2015-ductile (10.1007/s00466-015-1151-4); ambati-2016-verification
(10.1007/s00466-015-1225-3); borden-2016-triaxiality (10.1016/j.cma.2016.09.005);
miehe-2015-ductile (10.1016/j.cma.2014.11.017); miehe-2016-gradient
(10.1016/j.ijplas.2016.04.011); alessi-2018-comparison (10.1007/978-3-319-60885-3_1);
dittmann-2018-variational (10.1016/j.cma.2018.07.029); yin-2020-toughness
(10.1016/j.cma.2020.113068); hu-2021-coalescence (10.1007/s00466-021-02033-1);
han-2022-toughness (10.1007/s00466-021-02087-1); aldakheel-2018-gurson
(10.1007/s00466-017-1530-0); marengo-2023-review (10.1016/j.euromechsol.2023.105058);
delorenzis-2016-crystal (10.1002/gamm.201610002); hernandezpadilla-2017-crystal
(10.1007/s00161-015-0471-0); shanthraj-2016-viscoplastic (10.1016/j.cma.2016.05.006);
shanthraj-2017-cleavage (10.1016/j.jmps.2016.10.012); diehl-2017-void
(10.1007/s11837-017-2308-8); vakili-2022-aluminum (arXiv:2204.13994, flagged);
wilson-2018-driving (10.1016/j.jmps.2018.07.005); lucarini-2023-cracking
(10.1016/j.ijfatigue.2023.107670); tu-2020-series (10.1016/j.engfracmech.2020.106970);
maloth-2023-microstructures (10.1615/intjmultcompeng.2022042164); li-2024-cleavage
(10.1016/j.ijplas.2024.103888); auth-2024-thermodynamic
(10.1016/j.euromechsol.2024.105418); xue-2025-orientation (10.1557/s43578-025-01614-7);
shishvan-2021-mechanism (10.1016/j.engfracmech.2021.107927); sur-2025-hydrostatic
(10.1007/s11012-025-01971-x); abrarivajari-2022-multiaxial (10.1016/j.cma.2022.115467);
duda-2015-gradient (10.1016/j.ijplas.2014.09.005); pham-2011-gradient
(10.1177/1056789510386852); lorentz-2011-cohesive (10.1016/j.crme.2010.10.010);
wu-2017-unified (10.1016/j.jmps.2017.03.015); peirce-1982-localized
(10.1016/0001-6160(82)90005-0); asaro-1979-geometrical (10.1016/0001-6160(79)90036-1);
simoto-1965-copper (10.1080/14786436508218874); jasienski-1985-instabilite
(10.1051/mattech/198573080091); zhang-2007-necking (10.1007/s11431-007-0038-9);
kocks-2003-hardening (10.1016/s0079-6425(02)00003-8); rezwan-2026-predicting.

Entries added to the bibliography in the same pass but not cited above:
rosakis-2000-partition (10.1016/s0022-5096(99)00048-4).
