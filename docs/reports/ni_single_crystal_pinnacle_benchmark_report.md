# Crystal plasticity and phase-field ductile fracture of single-crystal Ni tensile specimens with PINNACLE: finite-element calibration and the PINN+FEM, PINN and surrogate routes

Technical report, draft of 2026-10-23. Status: methodology (§2),
numerical method (§3), the crystal-plasticity parameter calibration (§4),
the mesh sensitivity (§5) and the two phase-field fracture calibrations,
without (§6) and with (§7) a plastic-work threshold, written for the
production mesh h020; the fifteen specimens and the neural and surrogate
routes follow as they are confirmed by the author; abstract and introduction are held until the
study is closed. Every number here traces to a run in `docs/run_ledger.md`
(rows `ni-*`) and to a deck under `decks/ni/`; the lab notebook is
`docs/12_ni_study.md`, the plan `docs/11_ni_plan.md`. Figures of §4 are
produced by `scripts/ni_calibration_figures.py` from
`decks/ni/figures/cal_figures.i`, those of §5 by `scripts/ni_mesh_figures.py`
and `scripts/ni_fields.py` from `decks/ni/figures/mesh_figures_h020.i` and
`decks/ni/figures/fields_mesh_h020.i`, those of §6 and §7 by
`scripts/ni_pfd_curve_figures.py`, `scripts/ni_pfd_search_figures.py` and
`scripts/ni_pfd_shape_figures.py` from `decks/ni/figures/pfd_*.i`.

## Abstract

*(to be written after all studies are complete)*

## 1. Introduction

*(to be written after all studies are complete)*

## 2. Crystal-plasticity finite-element methodology

The constitutive model is the one of the Al2219 report (Part I there),
written again here with the Ni conventions; the specimens are single
crystals, so there is one orientation per specimen and no grain structure.

### 2.1 Kinematics

The deformation gradient at each material point is split multiplicatively into
an elastic and a plastic part,

$$
\mathbf{F} = \mathbf{F}^e \mathbf{F}^p , \qquad
\mathbf{L}^p = \dot{\mathbf{F}}^p \mathbf{F}^{p-1} = \sum_{\alpha=1}^{12} \dot{\gamma}^\alpha \, \mathbf{s}^\alpha \otimes \mathbf{n}^\alpha ,
\tag{1}
$$

where $\mathbf{s}^\alpha$ and $\mathbf{n}^\alpha$ are the slip direction and
slip-plane normal of system $\alpha$ in the sample frame and
$\dot{\gamma}^\alpha$ its slip rate. The twelve $\{111\}\langle 110\rangle$
systems of the FCC lattice are used. The crystal frame of a specimen is set by
the Bunge triple $(\varphi_1, \Phi, \varphi_2)$ of its orientation group,
which rotates both the slip systems and the elastic stiffness.

### 2.2 Elasticity

The elastic response is cubic and hyperelastic in the Green-Lagrange strain of
the elastic part,

$$
\mathbf{E}^e = \tfrac{1}{2}\left(\mathbf{F}^{eT}\mathbf{F}^e - \mathbf{I}\right), \qquad
\mathbf{S} = \mathbb{C} : \mathbf{E}^e , \qquad
\boldsymbol{\sigma} = \frac{1}{\det \mathbf{F}^e}\, \mathbf{F}^e \mathbf{S} \mathbf{F}^{eT} ,
\tag{2}
$$

with $\mathbb{C}$ the cubic stiffness ($C_{11}$, $C_{12}$, $C_{44}$) rotated
into the sample frame, $\mathbf{S}$ the second Piola-Kirchhoff stress and
$\boldsymbol{\sigma}$ the Cauchy stress. The resolved shear stress on system
$\alpha$ is

$$
\tau^\alpha = \mathbf{S} : \mathbf{M}^\alpha , \qquad
\mathbf{M}^\alpha = \tfrac{1}{2}\left(\mathbf{s}^\alpha \otimes \mathbf{n}^\alpha + \mathbf{n}^\alpha \otimes \mathbf{s}^\alpha\right).
\tag{3}
$$

### 2.3 Flow rule

Slip follows a rate-dependent power law (Hutchinson form),

$$
\dot{\gamma}^\alpha = \dot{\gamma}_0 \left| \frac{\tau^\alpha}{\tau_c^\alpha} \right|^{n} \operatorname{sign}(\tau^\alpha), \qquad n = 1/m_s ,
\tag{4}
$$

with reference slip rate $\dot{\gamma}_0$, rate exponent $n$ and the rate
sensitivity $m_s$ of (rezwan 2026 predicting). No kinematic backstress is
used. With $\dot{\gamma}_0 = 10^{-5}$ s$^{-1}$ and $n = 8.33$, a slip rate of
0.2 s$^{-1}$ (single slip at the test's gage rate, §2.6) needs
$|\tau| = 3.3\,\tau_c$: at the strain rate of these tests the flow stress is
mostly rate stress (§4.6 gives its effect on the calibrated set).

### 2.4 Slip resistance and dislocation-density evolution

The slip resistance is a Taylor relation on the forest density, with the
latent-hardened interaction sum,

$$
\tau_c^\alpha = \tau_{c0} + \alpha_0\, \mu\, b \sqrt{\sum_{\beta=1}^{12} a^{\alpha\beta} \rho^\beta}, \qquad
a^{\alpha\beta} = \begin{cases} 1 & \alpha = \beta \\ q & \alpha \ne \beta \end{cases}
\tag{5}
$$

where $\rho^\beta$ is the statistically stored dislocation density of system
$\beta$, $\mu$ the shear modulus, $b$ the Burgers vector, $\alpha_0$ the
Taylor prefactor, $q$ the latent-to-self hardening ratio and $\tau_{c0}$ the
lattice friction (the paper's Peierls stress). With all twelve densities
equal, the interaction sum multiplies the density by $1 + 11q = 16.4$; the
calibrated parameters of §4 are specific to this convention, as in the
Al2219 report.

Each density evolves by a Kocks-Mecking storage-recovery law written per
system,

$$
\dot{\rho}^\alpha = \frac{1}{b}\left( k_m \sqrt{\rho^\alpha} - 2 y_c\, \rho^\alpha \right) \left| \dot{\gamma}^\alpha \right| ,
\qquad \rho^\alpha(0) = \rho_0 ,
\tag{6}
$$

with the storage coefficient $k_m$ (PINNACLE's `k1`) and the annihilation
distance $y_c$ (PINNACLE's `annihilation_diameter`, recovery coefficient
$k_2 = 2 y_c$). The law saturates at

$$
\rho_{sat} = \left( \frac{k_m}{2 y_c} \right)^2 , \qquad
\tau_{c,sat} = \tau_{c0} + \alpha_0 \mu b \sqrt{(1 + 11 q)\, \rho_{sat}} ,
\tag{7}
$$

so $k_m$ sets the hardening rate after yield and, together with $y_c$, the
height of the plateau; $\rho_0$ sets the initial resistance through Eq. (5);
$\tau_{c0}$ is a uniform offset of the resistance. §4.5 quantifies each of
these on the specimen.

No geometrically necessary densities or damage are active in the runs of
this part of the report: the specimen is the undamaged crystal-plasticity
baseline that the phase-field fracture stage builds on.

### 2.5 Material parameters

Table 1 lists every parameter of Eqs. (1)-(7), its value in the runs of §4,
its provenance tag (PROJECT_GUIDELINES §4) and its role in the study. Units
are millimeter, megapascal and second throughout, the paper's system, so
Table 3 of (rezwan 2026 predicting) enters the decks unchanged. PINNACLE's
law has no built-in length unit: the same specimen run in micrometers and in
millimeters gives the same stress (`docs/12_ni_study.md` §1).

**Table 1.** Crystal-plasticity parameters for single-crystal Ni.

| symbol | meaning | value | unit | provenance | role |
|---|---|---|---|---|---|
| $C_{11}$ | cubic stiffness | 253.0 | GPa | reported_table ((rezwan 2026 predicting), Table 3, from (lim 2025 quantitative)) | fixed |
| $C_{12}$ | cubic stiffness | 152.0 | GPa | reported_table | fixed |
| $C_{44}$ | cubic stiffness | 124.0 | GPa | reported_table | fixed |
| $\mu$ | shear modulus in Eq. (5) | 26.0 | GPa | reported_table | fixed |
| $b$ | Burgers vector | $2.5 \times 10^{-7}$ | mm | reported_table | fixed |
| $q$ | latent-hardening ratio | 1.4 | - | reported_table | fixed |
| $\alpha_0$ | Taylor prefactor | 0.4 | - | reported_table | fixed |
| $\dot{\gamma}_0$ | reference slip rate | $10^{-5}$ | s$^{-1}$ | reported_table | fixed |
| $m_s$ ($n = 1/m_s$) | rate sensitivity | 0.12 (8.333) | - | reported_table | fixed |
| $\rho_0$ | initial SSD density per system | $7.256 \times 10^5$ | mm$^{-2}$ | calibrated (§4); paper $2.35 \times 10^6$ | free |
| $k_m$ | storage coefficient | 0.01147 | - | calibrated (§4); paper 0.0295 | free |
| $y_c$ | annihilation distance | $1.542 \times 10^{-6}$ | mm | calibrated (§4); paper $0.79 \times 10^{-6}$ | free |
| $\tau_{c0}$ | lattice friction | 0.244 | MPa | calibrated (§4); paper 0.26 | free |
| $\dot{\varepsilon}$ | applied bar strain rate | $5 \times 10^{-2}$ | s$^{-1}$ | assumed (the paper does not state the test rate; author's choice of 2026-10-15) | fixed |

The fixed set is the author's decision: the elastic constants, $b$, $\mu$,
$q$, $\alpha_0$, $\dot{\gamma}_0$ and $m_s$ stay at the paper's Ni values;
the four free parameters are the ones the paper itself calibrated on
specimen L, and the ones calibrated for Al2219.

### 2.6 Single-crystal dog-bone specimen model

**Specimens.** (rezwan 2026 predicting) tests fifteen single-crystal Ni
dog-bone coupons: five orientation groups (9, 10, H, L, P; Table 1 of the
paper, `data/extracted/ni_table1_orientations.csv`) with three thickness
variants each (x.1, x.3, x.4), whose gage thickness was measured at five
stations 2 mm apart (Table 2 of the paper,
`data/extracted/ni_table2_thickness.csv`). The calibration of §4 uses
specimen L.1: orientation L, Bunge $(60.0°, 23.5°, 80.2°)$, maximum Schmid
factor 0.42, gage thickness 0.3399 / 0.3441 / 0.3540 / 0.3625 / 0.3715 mm at
the five stations (mean 0.3544 mm).

**Geometry and mesh.** The planform is that of Fig. 2 of the paper: 25.4 mm
long, 4.8 mm wide grips, circular fillets of radius 3.175 mm tangent to a
1.85 x 8.0 mm gage (`data/extracted/ni_fig2_planform.csv`; the gage width is
measured on the plan view with the grip width as scale, $\pm 0.05$ mm, since
the paper does not dimension it). The thickness follows the five measured
stations by a natural cubic spline inside the gage, blends linearly across
each fillet to a common 0.35 mm, and is constant in the grips (author's
decision: identical ends, the measured variation in the gage). The mesh is a
single mapped hexahedral block built by the `dogbone` generator added to
hexaform for this study (`scripts/ni_specimen_mesh.py decks/ni/mesh/h020.i`).
The production level is h020: 0.2 mm cells in the gage, 0.4 mm in the
grips, one HEX8 element through the thickness (pseudo-3D, as in the Al2219
plates), 88 x 9 x 1 divisions, 792 elements, 1,780 nodes, 5,340 displacement
degrees of freedom and 6,336 quadrature points; minimum scaled Jacobian
0.58, no inverted element (Fig. 1). The 8 mm gage carries 40 x 9 cells. The
mesh study (h040, h020, h010 and h020 with two elements through the
thickness) is reported in §5 with the calibrated set; the coarser h040 is 10
to 20% stiffer than h020 between 10 and 45% gauge strain, which is the reason
for the production level.
Node sets `gage-` and `gage+` sit on the first and last station planes and
carry the gauge strain.

**Figure 1.** The h020 mesh of specimen L.1 (blocks gage, shoulder and grip;
the five thickness stations dotted), its planform and thickness profile, and
the thickness profiles of the fifteen specimens.

![mesh](../../figures/ni/mesh_h020.png)

**Orientation.** One orientation per specimen, the Bunge triple of its
group, given to the crystal block of the material (`euler_angles`); the
sample frame is $x$ = loading direction, $y$ = transverse, $z$ = thickness,
as in the paper.

**Boundary conditions and loading.** The grip end $x = 0$ is clamped
($\mathbf{u} = 0$ on node set `x-`); the grip end $x = 25.4$ mm is clamped
in $y$ and $z$ and displaced along $+x$ by a linear ramp to 30% of the bar
length (7.62 mm) in 6.0 s, i.e. a bar strain rate of
$5 \times 10^{-2}$ s$^{-1}$ (1.27 mm s$^{-1}$). All other faces are
traction-free. Because the gage is the narrowest and thinnest section, it
carries most of the elongation: the gauge strain runs 1.8 times the bar
strain up to 20% gauge strain and 2.0 times by the end, so the gage strain
rate is about 0.09 to 0.10 s$^{-1}$.

**Stress and strain measures.** The reported strain is the gauge strain,
the FE counterpart of the paper's DIC gage average,

$$
\varepsilon = \frac{\langle u_x \rangle_{gage+} - \langle u_x \rangle_{gage-}}{\langle X \rangle_{gage+} - \langle X \rangle_{gage-}} ,
\tag{8}
$$

the difference of the mean axial displacements of the two station planes
over their initial distance (8 mm; PINNACLE `gauge_node_sets`, ADR 0032). The
reported stress is the engineering stress, the reaction force on the pulled
grip over the mean initial gage section,

$$
\sigma = \frac{R_x}{w_g\, \bar{t}} , \qquad \bar{t} = \tfrac{1}{5}\sum_{k=1}^{5} t_k ,
\tag{9}
$$

with $w_g = 1.85$ mm and $\bar{t}$ the mean of the five station thicknesses
(0.6556 mm$^2$ for L.1), as the paper defines its engineering stress.

**Reference data.** The ground truth is the experimental engineering
stress-strain curve of specimen L in Fig. 3 of (rezwan 2026 predicting),
digitized at the figure's native 300 DPI (axis calibration residual
$10^{-4}$ strain and 0.23 MPa; repeated-extraction error 0.21 MPa;
`data/digitized/ni_fig3_specimenL_calibration.csv`). Thirteen markers, from
0.5% to 50% gauge strain, are the median stress of the digitized columns
within $\pm 0.4$% strain of each target, with the spread of those columns
as their uncertainty (0.1 to 1.7 MPa, and 25 MPa at the 0.5% marker, which
sits on the elastic-plastic knee; `data/digitized/ni_fig3_experiment_markers.csv`).
The experimental curve ends at 51.3% strain. The paper does not say which of
L.1, L.3 or L.4 Fig. 3 shows (its Fig. 4(b) gives failure strains of 46.1 to
48.1% for the three); the curve is taken as the hardening reference and L.1
as the calibration geometry.

## 3. Numerical method (PINNACLE)

All runs use PINNACLE (branch `feat/pinn_benchmark`, revision d91ccdb or
later; ddc76c3 or later for every run with a load cutback), torch backend on the CPU, double precision, two threads per run on
the compute host, with a 16 GB resident-memory cap per run (§3.4). A run is fully described by one specimen deck
(`decks/ni/cp/*.i`), rendered by `scripts/ni_specimen.py` into a PINNACLE
deck from the template `decks/ni/templates/cp_specimen_template.i`; the
calibration iterations are study decks (`decks/ni/cp/cal_gn*.i`) run by
`scripts/ni_study.py`.

### 3.1 Finite-element discretization

Trilinear HEX8 elements with 2 x 2 x 2 Gauss quadrature carry the displacement
field. The finite-strain problem is posed on the displacement gradient
$\mathbf{H} = \nabla \mathbf{u}$ at each quadrature point, from which
$\mathbf{F} = \mathbf{I} + \mathbf{H}$; the material returns the Cauchy stress
of Eq. (2) and the assembler forms the residual of the balance of linear
momentum (`stress_divergence_cp` kernel, finite-strain mode) with Dirichlet
conditions on node sets. The reaction of Eq. (9) is the sum of the internal
force over the constrained degrees of freedom of the pulled face.

### 3.2 Load stepping

The quasi-static load program (`cp_quasi_static` executioner) divides the
6.0 s ramp to 30% bar strain into 300 equal steps of $\Delta t = 0.02$ s,
i.e. $10^{-3}$ bar strain and about $2 \times 10^{-3}$ gauge strain per
step. A step that fails to converge is halved and retried; a run has a
budget of ten halvings in all (PINNACLE `max_cutbacks`, counted over the
run, not per step). The time increment of the constitutive update follows
the applied load increment, so a halved step keeps the strain rate and the
curve continues past a cutback unchanged. A curve ends where the budget is
spent: the steps the executioner takes after that have not converged and are
not reported. The calibrated run converges at every step to 59% gauge strain
(305 steps, its five cutbacks all beyond 41%), past the last marker (50%)
and the end of the experimental curve (51.3%).

### 3.3 Newton solve and tangent

Each load step is solved by Newton's method on the displacement vector to an
absolute residual norm of $10^{-3}$ N over the free degrees of freedom (about
$6 \times 10^{-6}$ of the 160 N load at 50% strain), with no relative
criterion, at most 50 iterations, with a backtracking line search. The
algorithmic tangent $\partial \boldsymbol{\sigma} / \partial \mathbf{H}$ is
formed at every quadrature point by one-sided finite differences of the
integrated stress of §3.4 (nine perturbations of $\mathbf{H}$, perturbation
$10^{-6}$) and refreshed at every Newton iteration, then assembled sparse
and solved with a sparse LU factorization (`linear_method = direct`): the
finite-strain crystal-plasticity tangent is non-symmetric. The calibrated
run takes 19 Newton iterations per step on average (15 below 40% gauge
strain; 11 for the first step from the virgin state, 43 at most, near the
end).

### 3.4 Constitutive integration

The state $\{\mathbf{F}^p, \rho^\alpha\}$ at each quadrature point is
integrated over the load step inside every residual evaluation, by a
backward-Euler update of the slip increments (PINNACLE ADR 0033,
`state_update = implicit`): for the trial deformation $\mathbf{F}_{n+1}$ of
a Newton iterate, the slip increments $\Delta\gamma^\alpha$ solve

$$
\Delta\gamma^\alpha = \Delta t\, \dot{\gamma}^\alpha\!\left(\tau^\alpha(\mathbf{F}_{n+1}, \Delta\boldsymbol{\gamma}),\ \tau_c^\alpha(\rho_{n+1}(\Delta\boldsymbol{\gamma}))\right), \qquad
\mathbf{F}^p_{n+1} = \Big(\mathbf{I} + \sum_\alpha \Delta\gamma^\alpha\, \mathbf{s}^\alpha \otimes \mathbf{n}^\alpha\Big)\mathbf{F}^p_n ,
\tag{10}
$$

with $\rho_{n+1}$ advanced by Eq. (6) with $|\Delta\gamma^\alpha|$. A local
Newton per point solves Eq. (10) to $10^{-10}$ in the slip increments, with
the semi-analytic Jacobian
$\mathbf{I} + \Delta t\,(n|\dot{\gamma}|/|\tau|)\,\mathbf{K} + \Delta t\,(n\dot{\gamma}/\tau_c)\,\mathbf{h}\operatorname{sign}(\Delta\gamma)$
($\mathbf{K}$ the elastic Schmid coupling, $\mathbf{h}$ the Taylor-hardening
sensitivity), the material's exact residual at the trial end state,
residual backtracking, and a damped first phase for trials far from
equilibrium. When it fails, the step is taken as 2 to 64 backward-Euler
sub-increments along the straight deformation path, and only then by
path-following explicit sub-stepping. The stress of the step is the stress
of the integrated end state, so the slip of the step enters the equilibrium
of the same step; the state is committed when the step converges. The last
fallback is rare (a trial far from equilibrium after repeated cutbacks), but
once reached it has driven a run's memory to 16 GB within seconds in some
runs (the cause is not located); every run therefore carries a
resident-memory cap (16 GB, a cgroup limit), at which the run is stopped.
The steps it has already converged stand.

This replaces the explicit update of the Al2219 report, in which the stress
inside Newton was evaluated with the state frozen at the last converged step.
That scheme is first order in the load increment: on the Ni dog-bone with the
paper's parameters it gives 507 / 189 / 133 / 111 / 121 MPa at 3% gauge
strain with 16 / 64 / 128 / 256 / 1,024 steps to 2% bar strain, against a
converged 126 MPa (material point under true uniaxial stress, 3,000 explicit
steps), and would need about 30,000 steps on the calibration path. The
backward-Euler update reproduces a 1,024-step explicit path within 1% on a
single element at four steps of $10^{-3}$, and the dog-bone reference within 4% at
3% gauge strain and 2% at 4.75% with 64 steps (`docs/12_ni_study.md` §1.1;
tier 2, independently verified numerical agreement). On the Al2219 plate of
Part II the two schemes give the same calibrated fracture toughness within
1.5% (`docs/10_pfd_study.md` §10).

### 3.5 Output

The executioner records, for every accepted load step, the time, the
reaction on the pulled face, the gauge strain of Eq. (8), the Newton count
and the cutbacks; `scripts/ni_specimen.py` turns them into the engineering
curve of Eq. (9) (`results.csv`). Field output (displacement, stress,
strain, equivalent plastic strain, slip activity and mean density) is
written every 20 steps.

### 3.6 Reproducibility

Each run is reproducible from its specimen deck, the mesh deck and the data
files it names; the commands are in `docs/reproduce_ni_specimen.md`, and the
run ledger gives the deck, scripts, framework revision, host, cores and wall
time of every run. The calibrated h020 run (305 steps to 30% bar strain)
took 82 s per load step on average on two threads (68 s below 40% gauge
strain), 6.9 h in all, with five runs side by side.

## 4. Results: parameter calibration

### 4.1 Objective

The free parameters $\rho_0$, $k_m$, $y_c$ and $\tau_{c0}$ are those that
minimize the root-mean-square deviation between the engineering stress of
specimen L.1 on the h020 mesh, interpolated at the thirteen marker gauge
strains from 0.5 to 50% (§2.6), and the marker stresses, in this framework,
on this mesh and at the assumed strain rate. The published parameter table
is not a target. In PINNACLE's law (the latent-hardened sum of Eq. (5), the
per-system law of Eq. (6)) at the assumed rate, the Table 3 set of
(rezwan 2026 predicting) hardens about four times too fast: 942 MPa at 66%
gauge strain against the experiment's 240 MPa at 51%, an RMSE of 373 MPa
(on the preliminary mesh h040, whose excess stiffness of about 20% (§5) does
not account for the factor). The calibrated values are therefore reported as
such (provenance `calibrated`).

### 4.2 Procedure

1. Each iteration ran the center and one multiplicative perturbation per
   parameter ($\rho_0$, $k_m$, $y_c$ x 1.2; $\tau_{c0}$ x 1.5) as five
   concurrent runs of the specimen to 30% bar strain (study decks
   `decks/ni/cp/h020/cal_gn*.i`, driver `scripts/ni_study.py`). The
   finite-difference Jacobian $\partial \sigma_i / \partial \ln p_j$ at the
   marker strains gave a damped Gauss-Newton step in log parameters, limited
   to $\pm 0.7$ in $\ln p$ per iteration (`scripts/ni_gauss_newton.py`,
   which writes the next center and study deck).
2. The starting point was the set of a preliminary calibration on the
   coarser mesh h040. On h020 that set gives a curve 12 to 18% below the
   markers from 10% strain on (RMSE 29.7 MPa), the mesh effect of §5.
3. The iteration stopped when the next step predicted a gain below 1 MPa.
   Two iterations were needed (Fig. 2): the first center (RMSE 29.7 MPa)
   stepped to the final set with a predicted 3.6 MPa; the second center
   measured 4.1 MPa, and its own step (-0.7% on $\rho_0$, +5.9% on $k_m$,
   +12.5% on $y_c$, -3.1% on $\tau_{c0}$) predicted 3.4 MPa. That step moves
   $k_m$ and $y_c$ together, along the one direction the markers do not
   resolve (§4.3), and was not taken.
4. An iteration took 6 to 7 h. The center run cut back five times beyond
   41% gauge strain and the $\tau_{c0}$ run four times, all recovered. In the
   $k_m$ run every step beyond 41% needed one cutback, so it spent its budget
   of ten (§3.2) at 43.0% gauge strain; its curve, and the $k_m$ column of
   the Jacobian, end there. The step not taken is the only quantity that
   used that column beyond 40%. The run later reached the memory cap of §3.4
   and was stopped; its converged steps are unaffected.

**Figure 2.** Gauss-Newton convergence on h020: the center curve of each
iteration against the experimental markers (left) and the RMSE over the
thirteen markers, with the linear prediction of the previous step (right).

![convergence](../../figures/ni/cal/cal_convergence.png)

### 4.3 Calibrated parameter set

**Table 2.** Calibrated set on h020 (deck `decks/ni/cp/h020/calibrated_L1_h020.i`,
run `results/ni_study/cal_h020/gn2/base`), against the paper's values, with
the linearized one-standard-deviation uncertainty from the residual at the
eleven markers from 0.5 to 40% (where all four perturbation curves exist)
and the FE Jacobian, for the four-parameter fit and with $\tau_{c0}$ held.

| parameter | value | unit | paper | ratio | relative 1$\sigma$, four free | relative 1$\sigma$, $\tau_{c0}$ held |
|---|---|---|---|---|---|---|
| $\rho_0$ | $7.256 \times 10^5$ | mm$^{-2}$ | $2.35 \times 10^6$ | 0.31 | $\pm 400$% | $\pm 6$% |
| $k_m$ | 0.01147 | - | 0.0295 | 0.39 | $\pm 81$% | $\pm 20$% |
| $y_c$ | $1.542 \times 10^{-6}$ | mm | $0.79 \times 10^{-6}$ | 1.95 | $\pm 56$% | $\pm 36$% |
| $\tau_{c0}$ | 0.244 | MPa | 0.26 | 0.94 | $\pm 6{,}200$% | held |

The uncertainties say how tightly the markers pin each value, not how far
the curve would move:

- **$\tau_{c0}$ is not identifiable here.** At 0.24 MPa it is small against
  the Taylor term of Eq. (5), so it moves the curve by only 3 MPa per e-fold,
  uniformly along the strain, and its column of the Jacobian is nearly
  parallel to that of $\rho_0$ (correlation -1.0). That near-degeneracy
  inflates all four uncertainties of the four-parameter fit. It stays close
  to the paper's value, where it started.
- **$\rho_0$ is well determined once $\tau_{c0}$ is held** ($\pm 6$%): it
  alone sets the level at yield.
- **$k_m$ and $y_c$ are determined only as a combination.** With
  $\tau_{c0}$ held their correlation is 1.00: the markers fix
  $k_m / \sqrt{y_c}$ to $\pm 2$%, not the two separately, because the curve
  does not reach the plateau of Eq. (7) within 50% strain. The step not
  taken in iteration 2 moved along this direction.

**Figure 3.** Calibrated specimen against the experimental markers.

![calibrated curve](../../figures/ni/cal/cal_curve.png)

Characteristics of the calibrated curve: yield within the first load step
(76 MPa at 0.24% gauge strain; the load steps do not resolve the elastic
branch); engineering stress 87 / 104 / 134 / 179 / 211 / 233 / 248 MPa at
2 / 5 / 10 / 20 / 30 / 40 / 50% gauge strain; a hardening rate
$d\sigma/d\varepsilon$ of 0.55 to 0.59 GPa from 2 to 10% strain, falling to
0.37 GPa at 20%, 0.18 GPa at 40% and 0.12 GPa at 50%.

### 4.4 Error analysis

**Table 3.** Calibrated specimen at the experimental markers (engineering
stress against gauge strain).

| strain (%) | experiment (MPa) | FE (MPa) | difference (MPa) | relative error |
|---|---|---|---|---|
| 0.5 | 85.0 | 78.8 | -6.2 | -7.3% |
| 1 | 85.2 | 81.5 | -3.7 | -4.3% |
| 2 | 87.7 | 86.9 | -0.8 | -1.0% |
| 5 | 98.4 | 104.0 | +5.6 | +5.7% |
| 10 | 129.3 | 134.0 | +4.7 | +3.6% |
| 15 | 154.9 | 159.3 | +4.4 | +2.9% |
| 20 | 180.1 | 179.4 | -0.7 | -0.4% |
| 25 | 195.4 | 196.4 | +1.0 | +0.5% |
| 30 | 211.8 | 210.9 | -0.9 | -0.4% |
| 35 | 222.9 | 222.9 | +0.0 | +0.0% |
| 40 | 232.8 | 232.8 | +0.0 | +0.0% |
| 45 | 239.1 | 240.9 | +1.9 | +0.8% |
| 50 | 238.4 | 247.7 | +9.2 | +3.9% |

Root-mean-square deviation over the thirteen markers: 4.1 MPa; over the
ten markers from 2 to 45% strain: 2.8 MPa, with a largest deviation of
5.6 MPa (5.7%, at 5%). From 2% to 45% strain the specimen follows the
experimental curve within a few MPa, through the change from the near-linear
hardening of the first 15% to the flattening beyond 30%. This is a tier-4
agreement (PROJECT_GUIDELINES §5: agreement with experiment).

The deviations are not random. The FE curve is 4 to 6 MPa low at 0.5 and 1%
and 4 to 6 MPa high from 5 to 15%: the experiment hardens slowly between 2
and 5% (0.36 GPa) and faster between 5 and 10% (0.62 GPa), and the FE curve,
with a nearly constant rate of 0.55 to 0.59 GPa over that range, cannot
follow that delayed onset. No step of the four parameters removes the pattern
(the linear prediction of the next step is 3.4 MPa). The 0.5% marker sits on
the elastic-plastic knee, where the digitized columns spread by 25 MPa, so
its deviation is within the data's own uncertainty; at 1% the spread is
0.2 MPa.

The other miss is the 50% marker, 3.9% high. The experiment stops hardening
between 45 and 50% (239.1 and 238.4 MPa) as it approaches failure, and
crystal plasticity alone keeps hardening (0.12 GPa at 50%). The end of each
experimental curve is the onset of damage, which this part does not model;
the phase-field fracture stage calibrates its toughness on the failure strain
of Fig. 4(b) of (rezwan 2026 predicting) (46.3% for L.1), not on this marker.

**Figure 4.** Deviation of the calibrated specimen from the experimental
markers.

![error](../../figures/ni/cal/cal_error.png)

### 4.5 Effect of each parameter

Each parameter was perturbed one at a time from the calibrated set while the
others were held (Table 4, Fig. 5, Fig. 6). The Jacobian of the final
iteration is the same information as a sensitivity
$\partial \sigma / \partial \ln p$, i.e. the stress change per e-fold change
of the parameter, along the curve (Fig. 7). The $k_m$ curve ends at 43.0%
gauge strain, where its run had spent its cutback budget (§4.2); its
sensitivity is drawn to the 40% marker, and its RMSE is taken over the
eleven markers to 40%.

**Table 4.** One-at-a-time perturbations of the calibrated set on h020
(`decks/ni/cp/h020/cal_gn2.i`): engineering stress at the marker strains and
the sensitivity $\partial\sigma/\partial\ln p$.

| case | 2% | 10% | 20% | 30% | 40% | 50% | RMSE (MPa) | $\partial\sigma/\partial\ln p$ at 2 / 10 / 30% (MPa) |
|---|---|---|---|---|---|---|---|---|
| calibrated | 87 | 134 | 179 | 211 | 233 | 248 | 4.1 | |
| $\rho_0$ x 1.2 = $8.707 \times 10^5$ mm$^{-2}$ | 94 | 141 | 185 | 216 | 237 | 252 | 7.9 | 39 / 39 / 27 |
| $k_m$ x 1.2 = 0.01376 | 89 | 145 | 201 | 240 | 267 | (ends at 43.0%) | 21.4 (to 40%) | 10 / 59 / 161 |
| $y_c$ x 1.2 = $1.851 \times 10^{-6}$ mm | 86 | 129 | 169 | 196 | 214 | 226 | 11.7 | -3 / -26 / -83 |
| $\tau_{c0}$ x 1.5 = 0.366 MPa | 88 | 135 | 181 | 212 | 234 | 249 | 4.6 | 3 / 3 / 3 |

What each parameter does to the curve, in this framework and with this law:

- **$\rho_0$ sets the yield level.** Raising the initial density raises
  $\tau_c$ through Eq. (5) from the first slip onward, so the whole curve
  moves up, by 39 MPa per e-fold from yield to 10% strain and less beyond
  (27 MPa at 30%, 23 MPa at 50%), as the stored density overtakes the
  initial one. It is the parameter that places the knee of the curve.
- **$k_m$ sets the hardening rate.** Its effect is small at yield (10 MPa per
  e-fold at 2%) and grows with strain (59 MPa at 10%, 161 MPa at 30%, 187 MPa
  at 40%): a larger storage coefficient makes the curve climb faster after
  yield and raises the plateau of Eq. (7). It is the parameter that sets the
  slope of the plastic branch.
- **$y_c$ sets the saturation.** Recovery is negligible at low density and
  grows with it, so a larger $y_c$ lowers the late curve (-83 MPa per e-fold
  at 30%, -121 MPa at 50%) and leaves yield unchanged. It bends the curve
  toward its plateau and, with $k_m$, fixes the plateau height.
- **$\tau_{c0}$ is a small uniform offset.** It shifts the curve by about
  3 MPa per e-fold at every strain, because it adds to $\tau_c$
  independently of the density and is small against the Taylor term. On
  this specimen, at this rate, it is immaterial to the curve, unlike in the
  Al2219 plate, where it was interchangeable with a small change of
  $\rho_0$.

The order of the calibration follows from these shapes: level ($\rho_0$)
first, slope ($k_m$) second, curvature ($y_c$) third; $\tau_{c0}$ is carried
along.

**Figure 5.** One parameter perturbed at a time from the calibrated set (left
to right, top to bottom: $\rho_0$, $k_m$, $y_c$, $\tau_{c0}$). The $k_m$
curve ends at 43.0% gauge strain (cross), where its run spent its cutback
budget.

![one at a time panels](../../figures/ni/cal/cal_oat_panels.png)

**Figure 6.** All four perturbations on one axis.

![one at a time overlay](../../figures/ni/cal/cal_oat_overlay.png)

**Figure 7.** Sensitivity $\partial\sigma/\partial\ln p$ of the engineering
stress to each parameter along the curve, from the finite-difference Jacobian
of the final iteration ($k_m$ to the 40% marker).

![sensitivity](../../figures/ni/cal/cal_sensitivity.png)

### 4.6 Scope of the calibrated set

The set of Table 2 holds for: this framework's constitutive conventions (the
latent-hardened sum of Eq. (5), the per-system law of Eq. (6), the rate law
of Eq. (4) with the fixed $\dot{\gamma}_0$ and $m_s$); the assumed bar
strain rate of $5 \times 10^{-2}$ s$^{-1}$ (with $m_s = 0.12$ a factor 3 in
rate moves the flow stress by about 14%); the h020 mesh (0.2 mm elements;
the mesh dependence is quantified in §5); the geometry of specimen L.1 and
the orientation of group L, while the paper does not say which of L.1, L.3
or L.4 its Fig. 3 curve shows (§2.6); and the hardening branch up to the
onset of damage, which the phase-field stage adds on top of this set.

## 5. Mesh sensitivity

### 5.1 Design

The calibrated set of Table 2 was run on four meshes of specimen L.1 with
the same orientation, loading and solver criteria, changing only the
discretization: h040 (0.4 mm gage cells, the coarse reference), h020 (the
production mesh, which is the calibrated run itself), h010 (0.1 mm) and h020
with two elements through the thickness (Table 5; study deck
`decks/ni/cp/h020/mesh_L1.i`, driver `scripts/ni_study.py`). Each run follows
the 300-step path of §3.2 and is reported to its last converged step: h040
and h020 reach the end of the path, h010 and the two-layer mesh spend their
cutback budget at 47.1% and 44.0% gauge strain, past the last marker their
comparison uses (40% and 45%).

**Table 5.** Meshes of the sensitivity study and the engineering stress
(MPa) at six gauge strains; the last column is the change from h020 at
10 / 20 / 30 / 40% strain.

| mesh | elements | gage cell (mm) | through thickness | 2% | 10% | 20% | 30% | 40% | 45% | curve ends at | change from h020 | wall time (threads) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h040 | 220 | 0.4 | 1 | 92.0 | 147.2 | 203.7 | 244.7 | 275.8 | 288.5 | 57.3% | +9.9 / +13.5 / +16.0 / +18.5% | 1.7 h (2) |
| h020 | 792 | 0.2 | 1 | 86.9 | 134.0 | 179.4 | 210.9 | 232.8 | 240.9 | 59.0% | | 6.9 h (2) |
| h010 | 3,168 | 0.1 | 1 | 85.3 | 129.7 | 172.2 | 199.9 | 217.5 | 223.6 | 47.1% | -3.2 / -4.0 / -5.2 / -6.6% | 23 h (8) |
| h020, two layers | 1,584 | 0.2 | 2 | 86.7 | 132.3 | 177.9 | 209.1 | 230.7 | - | 44.0% | -1.3 / -0.8 / -0.9 / -0.9% | 8.4 h (4) |

### 5.2 Results

**Figure 8.** Stress-strain curves of the four meshes with the calibrated
set, against the experimental markers.

![mesh curves](../../figures/ni/mesh_h020/mesh_curves.png)

**Figure 9.** Engineering stress at fixed gauge strains against the gage
cell size (open symbols: h020 with two elements through the thickness).

![mesh convergence](../../figures/ni/mesh_h020/mesh_convergence.png)

The response softens monotonically with refinement at every strain (Figs. 8
and 9), and the change per halving of the cell shrinks: from h040 to h020
the curve falls by 9 to 16% between 10 and 30% gauge strain, from h020 to
h010 by 3 to 5%. The ratio of successive changes is 0.30 to 0.37 at every
strain from 2 to 45% (apparent order 1.5 to 1.8), so the sequence converges.
Richardson extrapolation from the three in-plane levels places h020 above the
limit by 2.7% at 2% strain, 4.9% at 10%, 6.0% at 20%, 8.4% at 30% and 11.4%
at 40%, and h010 by 0.9 to 4.0% over the same range. The second element
through the thickness changes the curve by less than 1.3% at every strain:
one layer resolves the thickness direction of these 0.33 to 0.40 mm thick
coupons.

**Figure 10.** Equivalent plastic strain on the deformed gage (plan view)
for h040, h020 and h010 at 8, 16 and 24% bar strain (one color scale per
column; the panel titles give the gauge strain and the peak value).

![plastic strain fields](../../figures/ni/mesh_h020/mesh_eqps.png)

At the specimen scale the three meshes deform alike: the same bar strain
gives the same gauge strain (15.7 to 16.2%, 29.8 to 30.1% and 45.1 to 46.7%)
and the gage takes the same S-shaped lateral bend by 24% bar strain, the
single-crystal shear of the L orientation. The finer mesh resolves diagonal
bands of concentrated slip where the gage meets the fillets, which h040
smears and h020 only begins to show; its peak equivalent plastic strain is
correspondingly higher (0.25 / 0.36 / 0.51 against 0.20 to 0.21 / 0.32 /
0.48). The stiffer response of the coarse meshes is consistent with these
unresolved bands; as in the Al2219 plate, the constitutive law has no
internal length, so localized slip stays mesh-dependent.

### 5.3 Mesh selection

The author's choice is h020 (Table 5): its curve is within 5% of the
extrapolated limit to about 10% gauge strain and within 11% to 40%, against
h010 at four times the elements, about four times the time per load step on
four times the threads (348 s on eight threads against 82 s on two) and 11 GB
of memory per run against 2.4 to 6 GB. The calibrated set of Table 2 is
therefore an h020 set: it absorbs the remaining discretization error, so it
reproduces the experiment on h020, and on a finer mesh the same parameters
give a curve up to 7% lower at 40% strain (h010). Every later Ni stage runs
on h020 with this set. h040 is the reference that shows the convergence and
is not used further.

## 6. Phase-field ductile fracture: calibration without a plastic-work threshold

### 6.1 Fracture model

The calibrated crystal plasticity of §4 on the h020 mesh is coupled to the
AT2 phase field of the Al2219 report (Part II there, Eq. (20)): the damage
$d$ obeys

$$
\frac{G_c}{\ell_0}\, d - G_c \ell_0 \nabla^2 d = 2 (1-d)\, \mathcal{H}, \qquad
\mathcal{H}(t) = \max_{\tau \le t} \left[ \chi\, \psi^+(\tau) + \zeta \left\langle W_p(\tau) - W_c \right\rangle \right] ,
\tag{10}
$$

with $\psi^+$ the positive part of the elastic energy density (Miehe's
spectral split of the elastic Green-Lagrange strain with the Voigt average
of the cubic constants), $\chi$ a prefactor on it ($\chi = 1$ in this
section, the form of the Al2219 report), $W_p =
\int \sum_\alpha \tau^\alpha \dot{\gamma}^\alpha \, dt$ the accumulated
plastic work density per unit reference volume (the slip of Eq. (4) against
the resolved shear stress of Eq. (3)), $\zeta = 1$, $W_c$ the plastic-work
threshold and $\eta = 0$ (quasi-static). The second Piola-Kirchhoff stress
of Eq. (2) is degraded by $g(d) = (1-d)^2$ before the resolved shear
stresses of Eq. (3) are formed, so a damaged point keeps its flow stress
and loses its load-carrying capacity only as $g \to 0$. Mechanics and damage are solved in a staggered loop of up to three
alternations per load step (tolerance $10^{-3}$ on the damage increment)
along the 300-step path of §3.2, with $d = 0$ held on the grip faces. This
section reports the condition $W_c = 0$, $\chi = 1$, in which the plastic
work drives damage from the first load step; §7 reports the thresholded
condition.

**Failure criterion.** Failure is the first load step at which the nodes
with $d > 0.85$ cut the specimen in two: the intact nodes are joined along
element edges and the step counts as failed when no path of intact nodes
joins the two grip faces (PINNACLE `stop_separation`; the run ends there).
The failure strain is the gauge strain of that step, or of the last
converged step when the section loses its load in the same step. This is
the paper's criterion, $d > 0.85$ across the section (rezwan 2026
predicting); its $d \ge 0.8$ at any node is 0.5 to 0.8% strain earlier on
every calibrated run. The calibration target is the experimental failure
strain of specimen L.1, 0.463 (the experiment axis of the paper's Fig. 4(b)),
accepted within $\pm 0.004$ (two load steps). The DIC gauge strain is
reported throughout, as in §4; the experiment is shown up to its failure
strain.

### 6.2 Calibration of $G_c$ at each length scale

The length scale $\ell_0$ is a numerical parameter tied to the mesh, not a
material property, so it was fixed by the author's rule for the mesh at
hand: on h020 (gage cell $h = 0.2$ mm) $\ell_0$ is moved until the response
stops changing, and $G_c$ is the material parameter calibrated there. $G_c$
was first searched at $\ell_0 = 0.4$, 0.6 and 0.8 mm ($\ell_0/h = 2$, 3, 4)
by regula falsi in $\ln G_c$ on the failure strain (driver
`scripts/ni_pfd_calibrate.py`, decks `decks/ni/pfd/cal_L1_h020*.i`; six to
seven runs per length scale, about 5 to 6.6 h each on two threads).

**Figure 11.** The $G_c$ searches: failure gauge strain against $G_c$ at each
length scale (left) and against $G_c/\ell_0$ (right); the band is the
accepted tolerance around the target. Open squares are runs that ended
without the $d > 0.85$ zone spanning the section (one lost its load at 44.0%
and its mechanics stopped on the next step; one was killed at 47.7%).

![G_c searches](../../figures/ni/pfd/search_L1_cal.png)

**Table 6.** Calibrated $G_c$ at each length scale ($W_c = 0$, h020), with
the failure strain, the position of the fracture zone along the gage
(0 = the gage end nearer the clamped grip; the thinnest station is at 0.20
and the experiment's crack at 0.20), the width of the $d > 0.85$ zone on the
top face at failure, and the peak engineering stress.

| $\ell_0$ (mm) | $\ell_0/h$ | $G_c$ (MPa mm) | $G_c/\ell_0$ (MPa) | failure strain | zone center (of the gage) | zone width (mm) | peak stress (MPa) at strain |
|---|---|---|---|---|---|---|---|
| 0.4 | 2 | 21.1 | 52.8 | 0.461 | 0.17 | 2.2 | 226.0 at 44.8% |
| 0.6 | 3 | 30.95 | 51.6 | 0.460 | 0.21 | 2.8 | 225.4 at 44.4% |
| 0.8 | 4 | 40.75 | 50.9 | 0.463 | 0.25 | 3.2 | 225.4 at 44.6% |
| 1.2 | 6 | 61.1 | 50.9 | 0.472 | 0.32 | 3.4 | 226.8 at 45.7% |

The failure strain depends on $G_c$ and $\ell_0$ only through their ratio
(Fig. 11, right): at equal $G_c/\ell_0$ the three searched length scales
fail within 0.01 to 0.03 of each other, and the calibrated ratio is 51 to
53 MPa at every $\ell_0$. In an unnotched coupon the damage has to
nucleate, and the AT2 model nucleates at a driving force proportional to
$G_c/\ell_0$; the length scale therefore sets the material strength as
much as $G_c$ does, and a sweep of $\ell_0$ at fixed $G_c$ could not level
off. The plateau check was accordingly run with $G_c/\ell_0$ held: the
$\ell_0/h = 6$ point (1.2 mm, $G_c = 61.1$) is the last row of Table 6.

### 6.3 Length-scale plateau and the carried pair

**Figure 12.** Stress-strain curves (top) and maximum phase field (bottom)
of the four calibrated pairs of Table 6 against the experiment, with the
calibrated crystal plasticity without damage (dashed) and the target
(dotted). Crosses mark failure.

![plateau curves](../../figures/ni/pfd/curves_L1_plateau.png)

**Figure 13.** Axial Green-Lagrange strain $E_{xx}$ on the deformed top face
at failure for the four pairs, with the $d > 0.85$ zone outlined in black,
against the experiment's DIC panel of L.1 (axial strain at the onset of
failure, the paper's Fig. 5; only the experiment is shown).

![plateau shapes](../../figures/ni/pfd/shape_L1_plateau.png)

**Figure 14.** The phase field itself at failure for the same runs, with
the $d = 0.85$ contour in cyan.

![plateau damage](../../figures/ni/pfd/shape_L1_plateau_damage.png)

At fixed $G_c/\ell_0$ the stress-strain curve and the history of the maximum
phase field do not depend on $\ell_0$ at all: the four curves of Fig. 12
agree within 0.3% of stress at every strain from 5 to 43%, and the failure
strain is flat over $\ell_0/h = 2$ to 4 (0.460 to 0.463) and rises to 0.472
at $\ell_0/h = 6$. What the length scale changes is the fracture zone
(Figs. 13 and 14): the $d > 0.85$ band widens from 2.2 to 3.4 mm and its
center drifts from 0.17 to 0.32 of the gage, away from the thinnest station
and toward the middle, and at $\ell_0/h = 6$ the band is wide and shallow
enough that the section separates with a maximum $d$ of only 0.88. By the
author's rule, the first two consecutive values of $\ell_0/h$ that barely
change the result are 2 and 3 (curve, failure strain and zone center within
0.04 of the gage), and the larger is kept: the carried pair of the $W_c = 0$
condition is

$$
\ell_0 = 0.6 \ \mathrm{mm} \ (\ell_0/h = 3), \qquad G_c = 30.95 \ \mathrm{MPa\,mm} \ (31 \ \mathrm{kJ/m^2}),
$$

failure at 0.460, the zone at 0.21 of the gage. $\ell_0 = 0.6$ mm is fixed
for every later phase-field stage of this study.

### 6.4 What the $W_c = 0$ condition reproduces, and what it does not

**Figure 15.** Every run of the three searches ($\ell_0 = 0.4$, 0.6 and
0.8 mm), stress-strain curves above and maximum phase field below, for
reference.

![all searches](../../figures/ni/pfd/curves_L1_all.png)

The calibrated pair reproduces the failure strain and the fracture location
of L.1: 0.460 against 0.463, and the zone at 0.21 of the gage against the
experiment's crack at 0.20, with flanks within a few degrees of the normal
to the axis. It does not reproduce the stress-strain curve or the shape of
the fracture. With $W_c = 0$ the plastic work drives damage from the first
load step, so the phase field grows from the start and is diffuse: the
maximum $d$ is 0.29 at 10% gauge strain, 0.51 at 20%, 0.65 at 30% and 0.75
at 40% (Fig. 12, bottom; the same numbers at every $\ell_0$), the whole
gage is above $d = 0.6$ at failure, and the grips reach $d = 0.1$ to 0.4.
The degradation of this diffuse field lowers the stress below the
crystal-plasticity curve by 0.9% at 10%, 2.3% at 20%, 3.8% at 30% and 5.2%
at 40% gauge strain, so the curve sits 2 to 6% below the experiment it was
calibrated to without damage, peaks at 225 MPa against the experiment's
239, and leaves the experiment from about 20% strain. The load is lost
only in the last one or two steps, from 216 MPa to 157 MPa, when $g(d)$
reaches zero across the section; the experiment holds its stress to
failure. The fracture zone is a band 2.8 mm wide, nearly normal to the
axis, where the experiment and the paper's simulation show a narrow crack
inclined at about 55 to 60 degrees (Fig. 13).

The growth of the maximum phase field is the signature of what is missing.
It rises steeply at first and flattens from 20 to 30% strain (Fig. 12,
bottom; Fig. 15), because the plastic work accrues at a nearly constant
rate while the AT2 response to it saturates, $d \approx 2\mathcal{H}\ell_0/G_c /
(1 + 2\mathcal{H}\ell_0/G_c)$; a crack that nucleates late and runs would
instead show $d$ near zero for most of the path and an exponential rise just
before failure. In the Al2219 polycrystal of the Part II report the same
$W_c = 0$ coupling needed no threshold, because the plastic work
concentrates at the grain boundaries and the hole and the rest of the plate
stays undamaged; in a single crystal there is no neighboring grain to
compete with, the plastic work is nearly uniform along the gage, and all of
it is counted toward damage. This is the reason for the plastic-work
threshold $W_c$ of §7, which holds the phase field at zero until the plastic
work at a point crosses it, and for the prefactors of Eq. (10) on the two
energies. The $W_c = 0$ condition is kept as the reference against which
§7 is read.

### 6.5 Cost

A coupled run takes 82 s per load step on two threads for the mechanics
(19 Newton iterations) plus the three damage solves, 5 to 6.6 h to failure
at 46% gauge strain, within a 16 GB resident-memory cap; the three
searches took 27 h of wall time side by side, and the four calibrated
runs are reported from their own snapshots (every second step, and the
failure step).

## 7. Phase-field ductile fracture: calibration with a plastic-work threshold

### 7.1 Design

The $W_c = 0$ condition of §6 reproduces the failure strain but not the
path to it: the phase field grows from the first load step, the curve sits
2 to 6% under the experiment from 20% strain on, and the fracture zone is a
band normal to the axis (§6.4). The threshold stage asks the model to
follow the experimental curve for as long as possible before fracture and
then to fail at 0.463. The calibrated crystal plasticity without damage
matches the experiment within 0.8% to the 45% marker, so the plastic work
density at that strain is the natural scale of $W_c$: in the $W_c = 0$
check run the gage-peak plastic work is 61 MPa at 38%, 73.5 at 44% and 76
at 46.3%. A threshold above 76 MPa cannot fail at the target, so $W_c$ was
bracketed at 55, 65 and 75 MPa (damage onset near 36, 41 and 46% strain)
with $G_c$ searched to the failure strain at each value, at the fixed
$\ell_0 = 0.6$ mm of §6.3 (decks `decks/ni/pfd/wc*_L1_h020.i`, driver
`scripts/ni_pfd_calibrate.py`, one search per threshold side by side).

**Failure definition.** The separation criterion of §6.1 is kept, with one
addition for a sharp crack: when the cracking step's mechanics stop
converging before the $d > 0.85$ zone spans the section, the run counts
as failed at its last converged step if the maximum $d$ has reached 0.85
and the load has fallen by more than 3% from its peak (the section is
losing its load through the crack). This fires on one run of this section,
the carried one.

### 7.2 The elastic driving energy as a barrier ($\chi = 1$)

The first round kept the driving force of §6 ($\chi = 1$). Every $G_c$
small enough to fail by 0.463 cracked the fillets at 4 to 14% strain,
before any plastic work had reached $W_c$ (Table 7, Fig. 16): with the
threshold on the plastic work only, the elastic energy $\psi^+$ drives
damage from the first step, and the AT2 model has a tensile strength
$\sigma_c = \sqrt{27 E G_c / (256\,\ell_0)}$ (tanne 2018 crack) that falls
to 100 to 375 MPa at the $G_c/\ell_0$ of 0.5 to 7 MPa these runs needed,
below the stress the fillet transition carries (the check on the coupled
cube in `docs/12_ni_study.md` §6.4 confirms the coupler passes the stored
elastic energy itself; the fillet energy is real). At the one pair that
nearly reached the target ($W_c$ 65, $G_c$ 8) the damage localized at the
fillet, not in the gage. The threshold did what it was meant to do where
the fillets held: at $W_c$ 55 and $G_c$ 16 the curve lies on the
experiment within 1% to 45% (peak 241 MPa against 239), $d$ stays below
0.1 to 30% and rises from a knee at 37%; the run had not failed by 48%.

**Table 7.** First round, $\chi = 1$, $\ell_0 = 0.6$ mm. Onset is the strain
at which the maximum $d$ passes 0.1; location as in Table 6.

| $W_c$ (MPa) | $G_c$ (MPa mm) | $G_c/\ell_0$ (MPa) | end of run | onset | peak stress (MPa) | where |
|---|---|---|---|---|---|---|
| 75 | 0.3 | 0.5 | separated at 0.037 | 0.3% | 94 | fillet (-0.10 of the gage) |
| 75 | 1.2 | 2.0 | mechanics lost at 0.092, $d$ 0.67 | 2% | 127 | fillet |
| 65 | 2 | 3.3 | separated at 0.138 | 6% | 148 | fillet (-0.16) |
| 65 | 8 | 13 | mechanics lost at 0.458, $d$ 0.998 (28 nodes above 0.85) | 17% | 233 | fillet (-0.10) |
| 55 | 4 | 6.7 | mechanics lost at 0.122, $d$ 0.15 | - | 145 | - |
| 55 | 16 | 27 | no failure by 0.480, $d$ 0.72 | 30% | 241 | gage, diffuse |

**Figure 16.** First round ($\chi = 1$): stress-strain curves (top) and
maximum phase field (bottom) at each threshold, with the $W_c = 0$ pair of
§6.3 in gray.

![W_c round 1 curves](../../figures/ni/pfd/curves_L1_wc_round1.png)

The elastic energy therefore picks the fillet and the plastic work the
gage, and the two cannot be balanced by $G_c$ alone. PINNACLE gained a
prefactor $\chi$ on the elastic energy in Eq. (10) (the author's MOOSE
CP-PFD carries such prefactors on both energies); $\chi = 0$ is the purely
ductile driving force, in which $W_c$ alone sets the onset and the plastic
work alone picks the site, and $\chi = 1$ keeps every earlier deck as it
was.

### 7.3 Calibrated $G_c(W_c)$ with the purely ductile driving force ($\chi = 0$)

**Table 8.** Calibrated pairs at $\chi = 0$, $\ell_0 = 0.6$ mm, h020. Onset
as in Table 7; stress ratios against the experiment markers at 20 / 30 /
40 / 45%; the fracture zone from the top face at failure; the flank angle
from the principal direction of the $d > 0.85$ nodes (for the narrow
$W_c = 75$ zone the crack line of Fig. 19 is the better descriptor).

| $W_c$ (MPa) | $G_c$ (MPa mm) | failure strain | mode | onset | stress / experiment | peak (MPa) at strain | zone center | zone width (mm) | flank |
|---|---|---|---|---|---|---|---|---|---|
| 55 | 3.255 | 0.458 | separation | 35.8% | 0.996 / 0.996 / 0.986 / 0.984 | 234.4 at 45.8% | 0.39 | 2.5 | 55 deg |
| 65 | 1.384 | 0.459 | separation | 40.7% | 0.996 / 0.996 / 1.000 / 0.978 | 234.9 at 43.1% | 0.41 | 2.3 | 58 deg |
| 75 | 0.1 | 0.464 | load drop (5.9%) | 45.7% | 0.996 / 0.996 / 1.000 / 1.008 | 241.5 at 45.5% | 0.38 | 0.6 | inclined (line) |

**Figure 17.** Stress-strain curves (top) and maximum phase field (bottom)
of the three calibrated pairs against the experiment, with the $W_c = 0$
pair of §6.3 in gray. Crosses mark failure.

![W_c calibrated curves](../../figures/ni/pfd/curves_L1_wc_cal.png)

**Figure 18.** Every run of the three $\chi = 0$ searches, one panel per
threshold, for reference.

![W_c chi 0 searches](../../figures/ni/pfd/curves_L1_wc_chi0.png)

**Figure 19.** Axial Green-Lagrange strain $E_{xx}$ on the deformed top face
at failure for the three calibrated pairs, with the $d > 0.85$ zone
outlined in black and the crack line (for each node row across the gage
width, the node of largest $d$; colored by that $d$) against the
experiment's DIC panel.

![W_c calibrated shapes](../../figures/ni/pfd/shape_L1_wc_cal.png)

**Figure 20.** The phase field at failure for the same runs, with the
$d = 0.85$ contour in cyan and the crack line.

![W_c calibrated damage](../../figures/ni/pfd/shape_L1_wc_cal_damage.png)

With the purely ductile driving force every threshold reproduces the
experiment to failure: the three curves are indistinguishable from each
other and from the markers within 1 to 2% from 20 to 45% strain, peak at
234 to 242 MPa against the experiment's 239, and the phase field is exactly
zero until the plastic work crosses $W_c$ (36, 41 and 46% strain) and then
rises steeply to failure within 10, 5 and 0.7% strain of its onset (Fig.
17). The damage sits in the gage and the crack is inclined at 48 to 60
degrees across the width, at 0.38 to 0.41 of the gage (Figs. 19 and 20);
the experiment's crack is at 0.20 and 55 to 60 degrees. The three
thresholds differ in the toughness they need: $G_c$ falls thirty-fold from
55 to 75 MPa, because the plastic work that accrues beyond $W_c$ before
0.463 falls from about 17 to about 1 MPa, and the AT2 field reaches 0.85
only when that excess is several times $G_c/\ell_0$. The fracture zone
narrows accordingly, from 2.5 mm at $W_c = 55$ to 0.6 mm (three cells) at
75.

### 7.4 How and why the two calibrations differ

**Table 9.** The carried pair of each condition.

| | $W_c = 0$, $\chi = 1$ (§6) | $W_c = 75$ MPa, $\chi = 0$ (§7) |
|---|---|---|
| $G_c$ (MPa mm) at $\ell_0 = 0.6$ mm | 30.95 | 0.1 |
| what $G_c/\ell_0$ sets | the nucleation strength: failure strain follows $G_c/\ell_0$ alone (Fig. 11) | the energy that completes the cut once $W_p > W_c$ |
| damage onset (max $d$ > 0.1) | 3% | 45.7% |
| max $d$ at 20 / 40% strain | 0.51 / 0.75 | 0 / 0 |
| stress / experiment at 20 / 40% | 0.973 / 0.948 | 0.996 / 1.000 |
| peak stress (MPa) | 225 | 241.5 |
| failure strain | 0.460 | 0.464 |
| load at failure | lost in one to two steps, 216 to 157 MPa | 227 MPa, 6% below the peak, when the mechanics loses the cracking step |
| fracture zone at failure | band 2.3 mm wide, flanks within 3 degrees of the normal to the axis, at 0.21 of the gage; whole gage above $d$ = 0.6 | line 0.6 mm wide, inclined across the width, at 0.38 of the gage; the rest of the specimen at $d$ = 0 |
| mesh cost to failure | 231 steps, 6,600 Newton iterations | 243 steps, 4,600 Newton iterations |

The two conditions are different models of the same material, not two
values of one parameter. With $W_c = 0$ the whole plastic work drives
damage, so the phase field is a slowly saturating function of the
accumulated work everywhere in the gage (§6.4), the degradation of that
diffuse field is what lowers the curve, and $G_c/\ell_0$ acts as a strength
that decides when the saturating field reaches 0.85 across the thinnest
section; the toughness is large (31 kJ/m$^2$) because it must hold back
the work of the entire path. With $W_c = 75$ MPa nothing is damaged until
the work at a point exceeds the threshold, which in this single crystal
happens first at one row of the gage at 45.7%; the field there goes from 0
to 0.85 within 0.7% strain and the toughness needed to let it do so is
0.1 MPa mm, three hundred times smaller, because only the excess work
beyond the threshold, about 1 MPa, is available to pay for the crack. The
threshold replaces the dissipation that a polycrystal supplies through its
neighboring grains (the Al2219 plate needed none): in a single crystal the
plastic work is nearly uniform along the gage and all of it would count,
so the single-crystal coupling needs $W_c$ to concentrate the damage in
strain as well as in space. The location is the one feature neither
condition reproduces: the plastic work of the model peaks at 0.38 to 0.41
of the gage where the experiment cracked at 0.20, near the thinnest
measured station.

### 7.5 The elastic share at the carried pair

At the carried pair ($W_c$ 75, $G_c$ 0.1) the elastic prefactor was raised
to 0.3, 0.6 and 0.8 to see whether a growing elastic share moves the crack
(`decks/ni/pfd/wc75_gc0.1_chi*_L1_h020.i`). It does not move it along the
gage; it moves it to the fillet outright: $\chi = 0.8$ separated at 1.5%
strain at -0.04 of the gage, $\chi = 0.6$ at 2.1% at -0.06, and $\chi =
0.3$ lost its mechanics at 4.0% with $d = 0.95$ there, as the AT2 strength
of 33 to 54 MPa at $G_c/\ell_0 = 0.17$ MPa predicts: the larger the
elastic share, the earlier the fillet cracks. At the $G_c$ a late ductile
onset needs, the driving force has to be purely ductile.

**Figure 21.** The elastic share at the carried pair: curves and maximum
phase field for $\chi$ = 0, 0.3, 0.6 and 0.8.

![chi sweep](../../figures/ni/pfd/curves_L1_chi_sweep.png)

### 7.6 Carried set

Every later phase-field stage of this study runs with

$$
\chi = 0, \qquad W_c = 75 \ \mathrm{MPa}, \qquad G_c = 0.1 \ \mathrm{MPa\,mm}, \qquad \ell_0 = 0.6 \ \mathrm{mm},
$$

on h020 with the crystal-plasticity set of Table 2 and the failure
definition of §7.1. The $W_c = 0$ pair of §6.3 is kept as the reference
condition.

### 7.7 Cost

The $\chi = 0$ runs take 5 to 10 h each to failure on two threads; the three
searches needed fourteen runs (Fig. 18). Four runs of §6 and §7 were killed
by the 16 GB memory cap when the crystal-plasticity update of one point fell
back to sub-stepping; the cause was an autograd graph growing by one node
per sub-step under the torch backend (the assembler's displacement carries
grad), fixed in PINNACLE on 2026-10-24 (ADR 0033 addendum) by integrating
the state on detached deformations. The runs reported here were unaffected
in their results (a killed run was rerun or superseded; `docs/run_ledger.md`,
rows `ni-pfd-*`).
