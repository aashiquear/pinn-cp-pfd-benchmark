# 14. Part IV plan: phase-field crystallization benchmark

Started 2026-11-17. Guide papers: (granasy 2019 nucleation) (the review) and
(granasy 2005 spherulites). Aim: for every case of the two papers that PINNACLE can be made to
run, one standalone input deck, its expected outcome, the comparison with the published result,
and the report section with the mechanism. Route order per case: FE, then PINN+FEM, then PINN.
Every run keeps its field series; animations ship beside the figures (author, 2026-11-17).

This document is the plan and the capability record. Results go to
`docs/reports/phase_field_crystallization_pinnacle_benchmark_report.md`; reproduction steps to
`docs/reproduce_crystallization.md`.

## 1. Models of the two papers

| id | model | fields | source |
|---|---|---|---|
| M1 | single-field: F = int eps^2/2 \|grad phi\|^2 + w g(phi) + p(phi) Delta_f, Allen-Cahn dynamics, optional noise | phi | review Eq. 12-16 |
| M2 | M1 + wall: contact-angle boundary condition (Model I: n.grad phi = sqrt(2 Delta_f(phi))/eps cos(Psi); Model II: phi = const) | phi | review Eq. 17-19 |
| M3 | binary ideal solution (Warren-Boettinger): f = w(c) T g + p f_L + (1-p) f_S, conserved c | phi, c | review Eq. 20; spherulite paper Eq. 2.1-2.4 |
| M4 | M3 + orientation field: f_ori = H T [1-p(phi)] {x F0 + (1-x) F1}(xi0 \|grad theta\|), anisotropy s(vartheta, theta), noise | phi, c, theta | spherulite paper Eq. 2.1-2.6; review Eq. 27 |
| M5 | phase-field crystal (Swift-Hohenberg free energy, conserved dynamics) | psi | review Eq. 30, 36 |

## 2. Case matrix

Tier = accuracy tier of PROJECT_GUIDELINES §5 that the case can reach. "PINN" is the pure
network route; a case marked "no" has a structural reason given below the table.

| case | content | model | target | tier | FE | PINN+FEM | PINN |
|---|---|---|---|---|---|---|---|
| A1 | planar equilibrium interface: profile, gamma, 10-90 % width | M1 | closed form | 2 | yes | yes | yes |
| A2 | steady front velocity v = -6 M Delta_f delta (p = phi^2(3-2phi)); numerical v for the standard p | M1 | closed form | 2 | yes | yes | yes |
| A3 | critical nucleus in Ni: radial profiles, W*/kT, R* against undercooling | M1, spherical | review Fig. 3 (digitized), own shooting solution, CNT | 2, 3 | yes | yes | yes |
| A4 | effect of g and p on W* (standard, GL bcc, GL fcc) | M1, spherical | review Fig. 5, own shooting solution | 2, 3 | yes | - | - |
| A5 | sub- and supercritical seeds: a seed below R* melts, above it grows (the nucleus is the saddle) | M1, spherical, transient | A3 | 2 | yes | yes | yes |
| B1 | heterogeneous nucleus on a wall, Psi = 61.2 deg, Model I and II, W_het / W_cap against undercooling | M2, 2-D | review Fig. 6 | 2, 3 | yes | later | later |
| B2 | noise-induced nucleation and growth, transformed fraction, Avrami exponent p = 1 + d = 3 | M1 + noise, 2-D | JMAK theory; review Fig. 10 | 2 (statistical) | yes | no | no |
| C1 | orientational free energy f_ori against misorientation (k = 2, m = 3, n = 1/2, x = 0.2) | M4 material | spherulite Fig. 3 | 2, 3 | material check | - | - |
| C2 | single-crystal growth forms, isotropic and four-fold, S = 0.8 and 1.0 | M3 (+ M4 one grain) | spherulite Fig. 4 a, b, d, e | 3 (morphology), R(t) law | yes | no | no |
| C3 | polycrystalline dendrite and category 1 spherulite (reduced orientational mobility) | M4 | spherulite Fig. 4 c, f; Fig. 5 | 3 (morphology) | yes | no | no |
| C4 | needle to spherulite with supersaturation, 30 deg branching | M4, k = 2 | spherulite Fig. 6, 7 | 3 (morphology) | yes | no | no |
| C5 | transformation kinetics of impinging spherulites, p = 3.04 +/- 0.02 | M4 | spherulite Fig. 8 | 3 (statistical) | reduced size | no | no |
| C6 | heat-treatment sequences and the morphology gallery (Table II) | M4 | spherulite Fig. 9, 10 | 3 (morphology) | selected | no | no |
| D1 | PFC: phase diagram points, hexagonal nucleus | M5 | review Fig. 35, 37 | 3 | later | no | no |

Structural limits of the network routes:

- Noise-driven cases (B2, C2-C6). The published patterns are realizations of a stochastic
  equation: white noise is drawn per cell and per step. A network that maps (x, t) to the field
  has no random input and cannot represent one realization, and PINNACLE refuses
  `langevin_noise` on the PINN path by design. These cases are reported as FE-only, in the same
  way the coupled CP-PFD field solve is hybrid-only (PINNACLE ADR 0017).
- A4 repeats A3 with other functions; the network routes are shown on A3 only.

Sizes. The spherulite paper used an explicit finite-difference code on 500^2 to 5000^2 grids
(40 CPU-years in total). C2-C4 at 500^2 are within reach of an explicit scheme on one machine;
C5 (5000^2) and the 4000^2 run of Fig. 5 are reproduced at a reduced size with the statistics
stated.

## 3. PINNACLE capability against the cases (survey of 9548123, 2026-11-17)

Exists: `time_derivative + ac_bulk + ac_interface` on `pf_fem_transient` (bdf1, bdf2, convex
splitting; Newton with the exact block Jacobian), free-energy authoring (`parsed_free_energy`,
arithmetic of one variable), `langevin_noise`, `pf_contact_angle` (for phi in [-1, 1]),
`nucleation_inserter`, `feature_volume_fraction`, `pf_pinn` (strong / weak / energy forms, time
slabs, causal weights, hard transforms), split Cahn-Hilliard.

Gaps, in the order the cases need them:

| # | gap | needed by | kind |
|---|---|---|---|
| G1 | no field time series from `pf_fem_transient` (final `fields.vtk` only) | every case (animations) | output |
| G2 | no spherical or axisymmetric volume weight in the FE assembly or in `pf_pinn` | A3-A5 | new capability |
| G3 | `ac_bulk` on the PINN path accepts only three free-energy materials; a parsed (tilted) energy is refused | A1-A5 PINN | fix |
| G4 | `pinn_fem_hybrid` on an Allen-Cahn deck solves a steady Laplace problem (returned field is zero): no physical PINN+FEM route for phase-field | A1-A5 PINN+FEM | new capability |
| G5 | operator kernels (`time_derivative`, `ac_interface`, ...) assemble dense n x n matrices in a Python loop: about 2e4 nodes is the ceiling | B, C | performance |
| G6 | no steady nonlinear solve started from the IC (workaround: one backward-Euler step with dt = 1e10) | A1, A3, B1 | option |
| G7 | contact-angle BC only for phi in [-1, 1] and with a fixed (1 - phi^2) form; Models I and II for phi in [0, 1] | B1 | new BC form |
| G8 | noise amplitude has no dt or cell-volume scaling and cannot depend on phi | B2, C | extension |
| G9 | no free energy of two fields, no log terms, no field-dependent mobility: the binary model M3 | C | new kernels |
| G10 | no orientation field (theta, \|grad theta\|, anisotropy s(vartheta - theta)) | C | new kernels |
| G11 | no explicit scheme and no periodic boundaries in `pf_fem_transient` | C at 500^2 | new scheme |
| G12 | no sixth-order (PFC) operator | D | new kernels |

## 4. Stages

- Stage A (single-field model, cases A1-A5): gaps G1-G4, G6. FE, PINN+FEM and PINN on each.
- Stage B (walls and noise, B1-B2): gaps G5, G7, G8.
- Stage C (binary and orientation-field models, C1-C6): gaps G9-G11.
- Stage D (phase-field crystal): gap G12.

Each stage ends with its decks, figures, animations and report section, and with the PINNACLE
changes committed and pushed on `feat/pinn_benchmark` after their tests pass. A stage with a
parameter family (A3: nine undercoolings; C4: six supersaturations) runs one case first.

## 5. Reference solutions made independently of PINNACLE

`scripts/cryst_reference.py` (deck `decks/crystallization/reference/ni_nucleus_3d.i`): closed
forms of A1 and A2, and the shooting solution of the critical nucleus (A3, A4) in scaled
variables. Result (2026-11-17): with eps^2 and w proportional to T (gamma = gamma_eq T/T_f,
constant interface width) the standard model gives W*/kT = 99.2 / 67.7 / 47.9 / 34.8 / 25.8 /
19.5 at Delta T/T_f = 0.15 / 0.175 / 0.20 / 0.225 / 0.25 / 0.275. The 'standard PF' line of the
review's Fig. 3(c), digitized by `scripts/cryst_digitize_fig3c.py` (224 columns, half stroke
1.4 kT), reads 100.4 / 67.6 / 48.0 / 34.8 / 25.8 / 19.6: agreement within 1.1 % (tier 3 for the
reference). With eps^2 and w held at their melting-point values the barrier is 62 to 176 %
higher, so the published line uses the proportional scaling; the review does not say so for
the PF line (it does for the CNT line), and this is recorded as `derived`.

Model parameters of Ni (review Table 2, `reported_table`): gamma = 0.302 J/m^2 and
d_10-90 = 1 nm give eps^2 = 4.1234e-10 J/m, w = 1.59255e10 J/m^3, delta = 0.22756 nm;
H_f / v_m = 2.46182e9 J/m^3.

## 6. Status (2026-11-17)

Stage A: FE done for A1-A5; PINN and PINN+FEM done for A1 and A3 (0.20). Gaps closed in
PINNACLE `ddaeb91` (ADR 0037): G1, G2 (FE and stationary PINN; not `pf_pinn`), G3, G4 (stationary
states), G6 (kept as one bdf1 step with a large dt), and the conversion half of G5 (operators
stored sparse; their build is still dense). Open in stage A, awaiting the author: the network
routes over the undercooling family, the guess scan of the three routes, and the network
routes of the transient cases A2 and A5 (spherical term in `pf_pinn`, per-step predictor).

