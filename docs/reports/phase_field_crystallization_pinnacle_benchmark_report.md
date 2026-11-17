# Phase-field crystallization benchmark of PINNACLE (Part IV)

Status 2026-11-17: stage A (single-field model, cases A1 to A5) is complete on the FE route;
the stationary cases A1 and A3 are also done on the PINN and PINN+FEM routes. Stages B to D
(walls and noise, binary and orientation-field models, phase-field crystal) are planned in
`docs/14_phase_field_crystallization_plan.md` and have no results yet. Abstract and
introduction are written when all stages are done.

Guide papers: (granasy 2019 nucleation), cited below as "the review", and
(granasy 2005 spherulites). PINNACLE revision of every result: `ddaeb91` on
`feat/pinn_benchmark` (ADR 0037). Reproduction: `docs/reproduce_crystallization.md`.
Every deck named here is a complete input: `pinnacle run <deck>`.

## 1. The single-field model

The local state is one non-conserved order parameter, $\phi = 0$ in the liquid and
$\phi = 1$ in the crystal. The free energy of the review's Eq. 12, isotropic ($s = 1$), is

$$
F[\phi] = \int \Bigl\{ \frac{\epsilon^2}{2}\,|\nabla\phi|^2 + w\,g(\phi) + p(\phi)\,\Delta f \Bigr\}\,dV ,
$$

with the double well $g$, the interpolation function $p$ ($p(0) = 0$, $p(1) = 1$) and the
driving force $\Delta f = f_s - f_l < 0$ in the undercooled liquid. The standard choice is

$$
g = \tfrac14\,\phi^2(1-\phi)^2 , \qquad p = \phi^3(10 - 15\phi + 6\phi^2) .
$$

Relaxational (Allen-Cahn) dynamics, the review's Eq. 16 without noise:

$$
\dot\phi = -M_\phi\,\frac{\delta F}{\delta\phi}
         = M_\phi \bigl\{ \epsilon^2\nabla^2\phi - w\,g'(\phi) - p'(\phi)\,\Delta f \bigr\} .
$$

In PINNACLE the three terms are three kernels on one variable: `time_derivative`,
`ac_interface` (gradient coefficient $\epsilon^2$, mobility $M_\phi$) and `ac_bulk` bound to a
`parsed_free_energy` material that holds $w\,g + p\,\Delta f$ as an expression.

### 1.1 Closed forms used as targets

Planar interface at the melting point ($\Delta f = 0$). The first integral of the stationary
equation, $\tfrac12\epsilon^2\phi_x^2 = w\,g(\phi)$, gives for the quartic well

$$
\phi(x) = \frac{1}{1 + e^{x/\delta}}, \qquad \delta = \sqrt{2\epsilon^2/w}, \qquad
\gamma_{SL} = \int \epsilon^2\phi_x^2\,dx = \frac{\epsilon^2}{6\,\delta}, \qquad
d_{10\text{-}90} = 2\ln 9\;\delta .
$$

For a general well, $\gamma_{SL} = \epsilon\sqrt{2w}\int_0^1\!\sqrt g\,d\phi$ and
$d_{10\text{-}90} = (\epsilon/\sqrt{2w})\int_{0.1}^{0.9} d\phi/\sqrt g$; these two relations fix
$\epsilon^2$ and $w$ from the interface free energy and thickness for any $g$ of the review's
Table 1.

Steady front. With $p = \phi^2(3-2\phi)$ the travelling profile is the equilibrium logistic
and the velocity is

$$
v = -6\,M_\phi\,\Delta f\,\delta \;=\; -12\,M_\phi\,\Delta f\,\sqrt{\epsilon^2/(2w)} ,
$$

the expression quoted in the review after Eq. 16. For the standard $p$ there is no closed
form; the reference is the travelling-wave solution $\phi(x - vt)$ of
$\epsilon^2\phi'' + (v/M_\phi)\,\phi' = w g' + p'\Delta f$ found by shooting on $v$.

Critical nucleus. The stationary state with spherical symmetry and liquid far away solves the
review's Eq. 15,

$$
\phi'' + \frac{2}{r}\,\phi' = \frac{1}{\epsilon^2}\bigl\{ w\,g'(\phi) + p'(\phi)\,\Delta f \bigr\},
\qquad \phi'(0) = 0, \quad \phi(\infty) = 0 ,
$$

and its excess free energy is the nucleation barrier $W^* = F[\phi^*]$. It is a saddle point
of $F$. The classical droplet with a sharp interface gives $R^* = 2\gamma/|\Delta f|$ and
$W^*_{CNT} = \tfrac{16\pi}{3}\gamma^3/\Delta f^2$.

### 1.2 Material: nickel

Values of the review's Table 2 (`reported_table`): $\gamma_{SL,eq} = 0.302$ J m$^{-2}$,
$d_{10\text{-}90} = 1.0$ nm, $\Delta H_f = 17.29$ kJ mol$^{-1}$, $T_f = 1728.3$ K,
$\rho = 8357$ kg m$^{-3}$, molar mass 58.6934 g mol$^{-1}$; Turnbull's driving force
$\Delta f = \Delta H_f (T - T_f)/(T_f v_m)$.

Derived (`derived`): $v_m = 7.0233\times10^{-6}$ m$^3$ mol$^{-1}$,
$\Delta H_f/v_m = 2.46182$ GPa, and at $T_f$: $\epsilon^2 = 0.412338$ aJ nm$^{-1}$,
$w = 15.92548$ GPa, $\delta = 0.227560$ nm.

Temperature dependence. The review states that $\epsilon^2 \propto T$ and $w \propto T$ give
$\gamma_{SL} \propto T$ at constant thickness, and uses $\gamma = \gamma_{eq}T/T_f$ for the
droplet-model line of its Fig. 3; it does not say which scaling the 'standard PF' line uses.
This was settled by computation (section 3.1): the published line is reproduced with the
proportional scaling and missed by 62 to 176 % with constant $\epsilon^2$ and $w$. Every Ni
deck uses the proportional scaling (`derived`).

The decks use nm, aJ and GPa (1 aJ nm$^{-3}$ = 1 GPa) and $M_\phi = 1$, so time is in units
of $1/(M_\phi\,\mathrm{GPa})$.

## 2. Routes and what PINNACLE needed

| route | executioner | how a stationary state is obtained |
|---|---|---|
| FE | `pf_fem_transient` | one backward-Euler step with $\Delta t = 10^{12}$: Newton's method on the weak stationary equation, started from the initial condition |
| PINN | `pf_pinn_stationary` | a network $\phi_\theta(\mathbf x)$ is fitted to the initial condition, then trained on the squared stationary residual |
| PINN+FEM | `pf_pinn_fem_hybrid` | the briefly trained network is the start of the FE Newton solve, whose result is reported |

Spherical symmetry is a property of the mesh: `[Mesh] coordinate_system = spherical` puts the
volume element $4\pi r^2$ into every quadrature weight, so the weak form
$(\epsilon^2\phi', \psi')_{4\pi r^2} + (f', \psi)_{4\pi r^2} = 0$ is the radial equation above
with the symmetry condition at $r = 0$ natural, and the energy monitor returns $W^*$. In the
strong form the interface kernel adds $(2/r)\,\phi'$; the stationary network reads $r^2$, so
$\phi'(0) = 0$ holds for any weights.

Why a separate stationary network. A saddle of $F$ is unstable under the gradient flow, so a
time-dependent PINN cannot hold it. Every stationary state is a zero, hence a minimum, of
the loss $\langle(\delta F/\delta\phi)^2\rangle$. The liquid $\phi \equiv 0$ is a zero as well;
the loss is therefore divided by
$\min\bigl(1, \langle\phi^2\rangle/\langle\phi_{guess}^2\rangle\bigr)^p$ (`deflation_power`),
which leaves the nucleus a minimum and removes the liquid from the minima.

These capabilities did not exist at PINNACLE `9548123` and were added for this part
(ADR 0037, commit `ddaeb91`): the curvilinear volume element, the field time series of
`pf_fem_transient`, parsed free energies on the PINN path, the two executioners, and sparse
storage of the linear operators. The earlier `pinn_fem_hybrid` on an Allen-Cahn deck returns
a zero field and is not used.

## 3. Results

Accuracy tiers follow `PROJECT_GUIDELINES.md` section 5. The reference solutions of
`scripts/cryst_reference.py` are independent of PINNACLE (closed forms; SciPy DOP853 shooting,
relative tolerance $10^{-11}$): agreement with them is tier 2. Agreement with the review's
published curve is tier 3. No experiment enters this stage.

### 3.1 A3: critical nucleus in undercooled Ni (review Fig. 3)

Decks `decks/crystallization/single_field/a3_ni_nucleus_x<ΔT/T_f>_fem.i`, nine reduced
undercoolings from 0.10 to 0.30; radial line of 12 to 23 nm (ten equimolar radii),
$h = 0.005$ nm.

![A3](../../figures/crystallization/a3_ni_nucleus.png)

*Critical nucleus in Ni. Left: radial profiles (lines: reference; circles: PINNACLE FE).
Middle: barrier against undercooling, with the digitized 'standard PF' line of the review's
Fig. 3(c). Right: equimolar radius against the classical radius.*

| $\Delta T/T_f$ | $W^*/kT$ FE | reference | published (digitized) | droplet model | $\phi(0)$ | $R_e$ (Å) |
|---|---|---|---|---|---|---|
| 0.100 | 255.26 | 255.26 | outside the axis | 258.49 | 0.9989 | 22.57 |
| 0.150 | 99.23 | 99.22 | 100.4 | 102.47 | 0.9756 | 14.76 |
| 0.175 | 67.68 | 67.68 | 67.6 | 70.92 | 0.9452 | 12.67 |
| 0.200 | 47.87 | 47.86 | 48.0 | 51.06 | 0.9035 | 11.18 |
| 0.225 | 34.77 | 34.77 | 34.8 | 37.86 | 0.8551 | 10.07 |
| 0.250 | 25.80 | 25.79 | 25.8 | 28.72 | 0.8037 | 9.22 |
| 0.275 | 19.46 | 19.46 | 19.6 | 22.18 | 0.7519 | 8.54 |
| 0.300 | 14.90 | 14.90 | 15.7 (read at 0.295, the end of the line) | 17.37 | 0.7014 | 7.99 |

- FE against the reference (tier 2): $W^*$ within $2.2\times10^{-5}$ relative on all nine
  runs, $\phi(0)$ within $4\times10^{-5}$, $R_e$ within $2.2\times10^{-5}$; 3 to 5 Newton
  iterations, 0.8 to 1.1 s per run.
- Reference against the published line (tier 3): within 1.1 % between 0.15 and 0.275. The
  digitized stroke is 2.7 $kT$ thick (half thickness 1.4 $kT$, the extraction uncertainty;
  `data/digitized/cryst_review_fig3c_standard_pf.csv`, method in
  `scripts/cryst_digitize_fig3c.py`).
- The diffuse nucleus lies 1 to 14 % below the droplet model, and its centre is not bulk
  crystal at large undercooling ($\phi(0) = 0.70$ at 0.30), as in the review's Fig. 3(b).

The three routes at $\Delta T/T_f = 0.20$ (decks `a3_ni_nucleus_x0.200_{fem,pinn,hybrid}.i`;
reference $W^* = 0.913705$ aJ $= 47.865\,kT$):

![A3 routes](../../figures/crystallization/a3_routes_x0.200.png)

| route | initial guess | $W^*$ (aJ), seed 0 | relative error over 5 seeds: mean ± std, worst | wall time |
|---|---|---|---|---|
| FE | droplet fitted to the reference | 0.913726 | $+2.2\times10^{-5}$ (deterministic) | 0.8 s |
| PINN | classical droplet, $R = 2\gamma/\lvert\Delta f\rvert$ | 0.913756 | $(+6.6 \pm 5.8)\times10^{-5}$, $1.7\times10^{-4}$ | 8 to 11 s |
| PINN+FEM | classical droplet | 0.913726 | $+2.2\times10^{-5}$, the FE value on every seed | 9 to 14 s |

- The PINN+FEM result equals the FE result to all printed digits because the corrector is
  the same Newton solve; the network's field was within $1.0$ to $7.1\times10^{-4}$ (relative
  L2) of it and the corrector took 3 iterations on every seed.
- The network routes do not beat the FE route here in accuracy or time. Their measured
  property is the range of starts: with deflation the stationary network reached the nucleus
  from droplets of 0.36 to 1.56 classical radii (PINNACLE test deck, scaled units); without
  deflation the smallest start was lost (relative residual 0.13, reported by the run).
- Not yet measured on the Ni decks: the same scan for the plain FE Newton solve. One check
  was made: from the classical droplet (amplitude 1, radius 0.981 nm) the FE solve of the
  0.20 case converges in 5 iterations, so at this undercooling the FE route does not need
  the fitted guess. The GL-fcc cases of section 3.2 did need it.

### 3.2 A4: other double-well and interpolation functions (review Fig. 5)

Decks `a4_ni_nucleus_<well>_<interpolation>_x<...>_fem.i`: the GL-fcc well
$g = \phi^2(1-\phi^2)^2/6$ and the interpolations $p_{bcc} = \phi^3(4-3\phi)$,
$p_{fcc} = \phi^4(3-2\phi^2)$ of the review's Table 1, each with $\epsilon^2$ and $w$ matched to
the same $\gamma_{SL}$ and thickness.

![A4](../../figures/crystallization/a4_functions.png)

*Lines: reference. Circles: PINNACLE FE.*

| $\Delta T/T_f$ | standard | $g_{fcc}$ & $p_{PF}$ | GL $g_{bcc}$ & $p_{bcc}$ | GL $g_{fcc}$ & $p_{fcc}$ |
|---|---|---|---|---|
| 0.15 | 99.2 | 85.4 | 147.3 | 175.0 |
| 0.20 | 47.9 | 37.6 | 81.3 | 102.2 |
| 0.25 | 25.8 | 17.6 | 50.3 | 67.0 |
| 0.30 | 14.9 | 8.6 | 33.2 | 47.2 |

All twelve FE runs are within $2.2\times10^{-5}$ of the reference (tier 2). The ordering and
the sensitivity are those of the review's Fig. 5: the Ginzburg-Landau pairs give barriers 1.5 to 3.2 times that of the standard pair, and
changing $g$ alone (first two columns) lowers it by 14 to 43 %. The
published curves of Fig. 5 were not digitized; the comparison with that figure is
qualitative (ordering and magnitude at 0.15 to 0.30).

### 3.3 A1: planar interface

Decks `a1_planar_interface_{fem,pinn,hybrid}.i`. FE: an interface twice too wide relaxes over
200 steps of $\Delta t = 0.01$ (relaxation time $\delta^2/(M_\phi\epsilon^2) = 0.126$).

![A1](../../figures/crystallization/a1_planar_interface.png)

| route | $\gamma_{SL}$ (J m$^{-2}$) | profile, max deviation from the logistic | 10-90 % width (nm) | note |
|---|---|---|---|---|
| target | 0.302000 | | 1.0000 | inputs of the parameter derivation |
| FE (time march) | 0.302005 | $3.2\times10^{-6}$ | 1.0000 | interface stays at $x = 0$ |
| PINN | 0.302005 (5 seeds: $+1.7\times10^{-5} \pm 1.5\times10^{-6}$) | $1.5\times10^{-4}$ | 1.0001 | interface at $x = 0.081$ nm (seed 0) |
| PINN+FEM | 0.302005 on every seed | $4.7\times10^{-6}$ | 1.0001 | interface at $x = 0.047$ nm; 4 to 19 corrector iterations |

Limit of the stationary routes on this case: the stationary equation does not fix the
position of a planar interface (translation invariance), and the two Dirichlet ends pin it
only exponentially weakly. Both stationary routes return the correct profile and energy at a
position shifted by 0.05 to 0.08 nm from the guess; the hybrid's corrector needs up to 19
iterations because its matrix is nearly singular along that shift, and the reported
prediction distance (0.07 to 2.7 %) is mostly the shift. The time march has no such freedom.

### 3.4 A2: steady crystallization front

Decks `a2_front_velocity_{tanh,standard}_fem.i`, $\Delta T/T_f = 0.20$, 40 nm, $h = 0.01$ nm,
BDF2 with $\Delta t = 0.02$ to $t = 35$; velocity from a linear fit of the solid length
$\int\phi\,dx$ for $t > 5$.

![A2](../../figures/crystallization/a2_front_velocity.png)

| $p$ | $v/M_\phi$ FE (nm GPa) | reference | relative error |
|---|---|---|---|
| $\phi^2(3-2\phi)$ | 0.672180 | 0.672253 (closed form) | $-1.1\times10^{-4}$ |
| $\phi^3(10-15\phi+6\phi^2)$ | 0.644397 | 0.644490 (travelling wave) | $-1.4\times10^{-4}$ |

The free energy falls at the rate $v\,\Delta f$ to ten digits in both runs, which is the
energy balance of a steady front. Tier 2.

### 3.5 A5: seeds below and above the critical size

Decks `a5_seed_{sub,super}critical_fem.i`: the tanh droplet fitted to the nucleus of case A3
with its radius multiplied by 0.9 and by 1.1, marched in time with the spherical volume
element.

![A5](../../figures/crystallization/a5_seeds.png)

The smaller seed melts (solid volume 4.2 nm$^3$ to zero by $t \approx 2.6$), the larger one
grows without bound; both start below the barrier ($F = 0.890$ and $0.887$ aJ against
$W^* = 0.9137$ aJ) and their free energy only decreases. This is the saddle-point character
of the nucleus shown by the dynamics, and the reason a time march cannot produce case A3.
Consistency check of A3 by an independent route in the same code (tier 1).

### 3.6 Animations

`figures/crystallization/a1_planar_interface.gif`, `a2_front_velocity_standard.gif`,
`a5_seed_subcritical.gif`, `a5_seed_supercritical.gif`: the profile over time from the
`fields_step*.vtk` series every run writes (`snapshots_every`).

## 4. Limits of stage A

- Network routes cover the stationary cases A1 and A3 at one undercooling. A2 and A5 are
  time dependent; they need `pf_pinn`, which does not yet carry the spherical term (A5), and
  a per-step predictor for the hybrid. Not run.
- The PINN and PINN+FEM routes were not run over the nine undercoolings or the other
  functions (one case first; the family awaits the author).
- The range of starts of the three routes (guess scan) is measured for the stationary PINN
  on the scaled test deck only.
- Neural results are reported over 5 seeds; the cost of the development runs behind
  `pf_pinn_stationary` (about 15 exploratory trainings of 5 to 15 s) is not in the wall
  times of the tables.
- Hardware of every number: Intel Core i9-9980HK, 4 threads, float64, numpy backend for FE
  and torch 2.13.0 (CPU) for the network routes; no GPU on this machine.
