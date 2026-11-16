# Crystal plasticity and phase-field ductile fracture of single-crystal Ni tensile specimens with PINNACLE: finite-element calibration and the PINN+FEM, PINN and surrogate routes

Technical report, restarted 2026-11-15 and revised 2026-11-16. The
finite-element route is in §2 to §7 (model, total-Lagrangian method, crystal plasticity on L.1,
mesh check, one $G_c$, fifteen specimens); the PINN+FEM, PINN and surrogate routes are in §8.
Revision of 2026-11-16: PINNACLE read the three Euler angles in exchanged order, so every
earlier run simulated crystals other than the measured ones. The angles now follow Bunge's
convention (PINNACLE ADR 0036), which reproduces the Schmid factors of the paper's Table 1,
and every run of this report was repeated. The version with the exchanged angles, which read
groups 9 and H as a separate ductile-like class, is withdrawn; its figures are kept in
`archive/figures/ni_exchanged_angles/` and its results under `results/ni_study/tl/`.
The first attempt (draft of 2026-10-23) is withdrawn: PINNACLE's finite-strain force was
the Cauchy stress over the undeformed section, so the reported engineering stress was the
true stress and no specimen could neck, and the assumed strain rate was not quasi-static
(`.claude/FINDINGS.md`, 2026-11-15). That draft and its figures are kept unchanged in
`archive/docs/ni_report_first_attempt_2026-11-15.md` and
`archive/figures/ni_first_attempt/`; none of its numbers is used here. Every number of
this report traces to a run in `docs/run_ledger.md` and to a deck under `decks/ni/`.

## Abstract

The fifteen single-crystal nickel tensile specimens of (rezwan 2026 predicting)
are simulated with PINNACLE: dislocation-density crystal plasticity coupled to
an AT2 phase field driven by the elastic energy and the plastic work, in a
total-Lagrangian finite-strain balance. One parameter set serves every run.
Three crystal-plasticity values are fitted to the stress-strain curve of
specimen L at an assumed gauge strain rate of $10^{-3}$ s$^{-1}$ (RMSE 4.1 MPa
over thirteen markers), and one fracture energy, $G_c = 48.67$ MPa mm at
$\ell_0 = 0.6$ mm, is fitted to the failure strain of specimen L.1. The Euler
angles follow Bunge's convention, which reproduces the Schmid factors the
paper reports. On the finite-element route the RMSE in failure strain over the
fifteen specimens is 0.039, and eleven specimens are within 5%. The largest
errors are P.3 and P.4, predicted 14 to 17% late, and 9.3 and 9.4, predicted
11 to 12% early.

Three network routes are measured against that reference. A network proposal
for each Newton solve, trained on the finite-element residual (PINN+FEM),
reproduces the finite-element result on all fifteen specimens to six digits
and costs 25% more time, because a linear extrapolation already meets the
residual target on 88% of the load steps. Accepting proposals without a Newton
solve (PINN) shifts the failure strain of L.1 by 0.7% and removes 52% of the
Newton iterations; the accepted fields are mostly the extrapolation. A
surrogate trained on thirty synthetic specimens predicts the finite-element
failure strain of the fifteen real specimens, none of them training data, with
an RMS difference of 3.1%, in under one second per specimen, after 15.0 hours
of training runs.

## 1. Introduction

Rezwan and co-workers tested fifteen single-crystal nickel dog-bone specimens,
five crystal orientations with three thickness profiles each, and predicted
their strain to failure with crystal plasticity coupled to a phase-field
fracture model, one parameter set for all (rezwan 2026 predicting). The data
set is small, well documented and discriminating: the orientation changes the
failure strain from 0.38 to 0.63, and the thickness profile moves the fracture
along the gage. It is used here as the single-crystal benchmark of PINNACLE.

The report has two aims. The first is a finite-element reference made with
PINNACLE alone: the same constitutive and fracture model, fitted with as few
free quantities as the paper uses, and compared with the paper's experiments
only. The paper's own simulation results are not used as targets. The second is
a measurement of what neural networks add on this problem, on three routes that
PINNACLE offers: a network that proposes the start of each Newton solve
(PINN+FEM), the same network with its proposal accepted when the residual is
small enough (PINN), and a surrogate that maps orientation and thickness
directly to the stress-strain curve and the failure descriptors.

§2 states the model and §3 the numerical method. §4 fixes the crystal
plasticity on specimen L.1, §5 checks the mesh, and §6 fits the one fracture
energy. §7 gives the fifteen specimens on the finite-element route and reads
the remaining errors. §8 gives the three network routes and
compares them in accuracy and time. Every run comes from a deck under
`decks/ni/` and has a row in `docs/run_ledger.md`.

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
which rotates both the slip systems and the elastic stiffness. The
crystal-to-sample rotation is $\mathbf{R} = \mathbf{R}_z(\varphi_1)\,
\mathbf{R}_x(\Phi)\, \mathbf{R}_z(\varphi_2)$, the transpose of Bunge's
orientation matrix (PINNACLE ADR 0036).

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
used.

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
equal, the interaction sum multiplies the density by $1 + 11q = 16.4$.

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
$\tau_{c0}$ is a uniform offset of the resistance.

No geometrically necessary densities or damage are active in the runs of
this part of the report: the specimen is the undamaged crystal-plasticity
baseline that the phase-field fracture stage builds on.

### 2.5 Material parameters

Table 1 lists every parameter of Eqs. (1)-(7), its value,
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
| $\rho_0$ | initial SSD density per system | $2.35 \times 10^6$ | mm$^{-2}$ | reported_table; refitted in §4 (Table 2) | free |
| $k_m$ | storage coefficient | 0.0295 | - | reported_table; refitted in §4 (Table 2) | free |
| $y_c$ | annihilation distance | $0.79 \times 10^{-6}$ | mm | reported_table; refitted in §4 (Table 2) | free |
| $\tau_{c0}$ | lattice friction | 0.26 | MPa | reported_table | fixed |
| $\dot{\varepsilon}$ | gauge strain rate | $10^{-3}$ | s$^{-1}$ | assumed (the paper does not state the test rate; §4) | fixed |

Table 3 of (rezwan 2026 predicting) is the starting point. Two things are
set in §4 on specimen L.1: the strain rate of the test, which the paper does
not give (with the rate sensitivity $m_s = 0.12$ a decade of rate moves the
flow stress by 32%), and the three hardening parameters $\rho_0$, $k_m$ and
$y_c$, refitted by 14 to 41% because the framework and the rate differ from
the paper's.

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
mesh level is checked in §5.
Node sets `gage-` and `gage+` sit on the first and last station planes and
carry the gauge strain.

**Figure 1.** The h020 mesh of specimen L.1 (blocks gage, shoulder and grip;
the five thickness stations dotted), its planform and thickness profile, and
the thickness profiles of the fifteen specimens.

![mesh](../../figures/ni/mesh_h020.png)

**Orientation.** One orientation per specimen, the Bunge triple of its
group, given to the crystal block of the material (`euler_angles`); the
sample frame is $x$ = loading direction, $y$ = transverse, $z$ = thickness,
as in the paper. The convention is checked on the paper's Table 1: the maximum
Schmid factors on the loading axis are 0.436, 0.478, 0.444, 0.417 and 0.460
for groups 9, 10, H, L and P, against the published 0.44, 0.48, 0.44, 0.42 and
0.46.

**Boundary conditions and loading.** The grip end $x = 0$ is clamped
($\mathbf{u} = 0$ on node set `x-`); the grip end $x = 25.4$ mm is clamped
in $y$ and $z$ and displaced along $+x$ by a linear ramp. All other faces are
traction-free. The ramp time sets the strain rate (§4).

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

All runs use PINNACLE (branch `feat/pinn_benchmark`, revision 9548123 with the
total-Lagrangian balance of ADR 0035 and the Euler convention of ADR 0036), torch backend on the CPU, double
precision, four threads per run on the compute host. A run is fully described
by one specimen deck (`decks/ni/cp/bunge/*.i`, `decks/ni/pfd/bunge/*.i`,
`decks/ni/final/bunge/*.i`), rendered by `scripts/ni_specimen.py` into a PINNACLE
deck from the templates under `decks/ni/templates/`.

### 3.1 Finite-strain balance

Trilinear HEX8 elements with 2 x 2 x 2 Gauss quadrature carry the displacement
field. The balance of linear momentum is written on the reference
configuration (total Lagrangian). At each quadrature point the displacement
gradient $H_{iJ} = \partial u_i / \partial X_J$ gives
$\mathbf{F} = \mathbf{I} + \mathbf{H}$, the material returns the Cauchy stress
of Eq. (2), and the internal force is the reference-volume integral of the
first Piola-Kirchhoff stress,

$$
\mathbf{P} = J\, \boldsymbol{\sigma}\, \mathbf{F}^{-T}, \qquad
R_{n i} = \sum_q w_q\, \frac{\partial N_n}{\partial X_J}\, P_{iJ} ,
\tag{10}
$$

with $J = \det \mathbf{F}$. The reaction of Eq. (9) is therefore the true force
on the pulled grip: it carries the thinning of the section, so the engineering
stress passes through a maximum and the specimen can neck. The tangent is the
exact derivative of Eq. (10),

$$
\frac{\partial P_{iJ}}{\partial F_{kL}} = J \left( F^{-1}_{Lk}\, (\boldsymbol{\sigma}\mathbf{F}^{-T})_{iJ}
 + \frac{\partial \sigma_{im}}{\partial F_{kL}} F^{-1}_{Jm}
 - \sigma_{im} F^{-1}_{Jk} F^{-1}_{Lm} \right) ,
\tag{11}
$$

with the material tangent $\partial \boldsymbol{\sigma} / \partial \mathbf{F}$
from the constitutive update (§3.3); the first and last terms are the geometric
stiffness. Four unit tests fix the formulation: the layout of the gradient, a
force that follows the deformed section, objectivity
($\mathbf{R}(\mathbf{Q}\mathbf{F}) = \mathbf{Q}\mathbf{R}(\mathbf{F})$), and the
tangent against central differences of the residual (relative error below
$10^{-5}$ with active plastic flow).

The earlier balance of PINNACLE integrated the Cauchy stress with the reference
gradients and weights and is kept as `kinematics = legacy` for the decks of
Parts I and II. On the Ni specimen it reported the true stress as the
engineering stress: at 59% gauge strain 257 MPa where the force over the
initial section is 159 MPa. Every run of this report sets
`kinematics = total_lagrangian`.

### 3.2 Load stepping and Newton solve

The pulled grip moves by $10^{-3}$ of the bar length per load step, about
$2 \times 10^{-3}$ gauge strain. The ramp to 30% bar strain takes 600 s (300
steps) in the crystal-plasticity runs and the ramp to 40% takes 800 s (400
steps) in the fracture runs, which is a gauge strain rate of about
$10^{-3}$ s$^{-1}$ (§4). A step that fails to converge is halved and retried,
with a budget of ten halvings per run. Each step is solved by Newton's method
with a backtracking line search to an absolute residual of 0.05 N
($5 \times 10^{-4}$ of the 100 to 160 N reaction) or a relative residual of
$10^{-6}$. The tangent of Eq. (11) is assembled sparse and factorized directly;
it is not symmetric. A linear extrapolation of the last two converged steps
starts each solve. The fracture runs take 3.1 Newton iterations per load step
on average.

### 3.3 Constitutive integration

The state $\{\mathbf{F}^p, \rho^\alpha\}$ at each quadrature point is
integrated over the load step inside every residual evaluation by a
backward-Euler update of the slip increments (PINNACLE ADR 0033,
`state_update = implicit`): for the trial deformation $\mathbf{F}_{n+1}$ the
increments $\Delta\gamma^\alpha$ solve

$$
\Delta\gamma^\alpha = \Delta t\, \dot{\gamma}^\alpha\!\left(\tau^\alpha(\mathbf{F}_{n+1}, \Delta\boldsymbol{\gamma}),\ \tau_c^\alpha(\rho_{n+1}(\Delta\boldsymbol{\gamma}))\right), \qquad
\mathbf{F}^p_{n+1} = \Big(\mathbf{I} + \sum_\alpha \Delta\gamma^\alpha\, \mathbf{s}^\alpha \otimes \mathbf{n}^\alpha\Big)\mathbf{F}^p_n ,
\tag{12}
$$

with $\rho_{n+1}$ advanced by Eq. (6). A local Newton per point solves Eq. (12)
to $10^{-10}$; the stress of the step is the stress of the integrated end
state, and its consistent tangent follows from the implicit-function theorem on
Eq. (12). The state is committed when the load step converges.

### 3.4 Phase-field solve

The damage equation (§6.1) is linear in $d$ for a given history field and is
solved once per alternation with a sparse direct solve. Mechanics and damage
alternate up to three times per load step (tolerance $10^{-3}$ on the damage
increment), with $d = 0$ held on the grip faces.

### 3.5 Reproducibility

Each run is reproducible from its specimen deck, the mesh deck and the data
files it names. The run ledger (`docs/run_ledger.md`, rows `ni-bunge-*`) gives the
deck, scripts, framework revision, host, cores and wall time of every run of
this report.

## 4. Crystal plasticity on specimen L.1

**Objective.** The crystal-plasticity set is fixed on the experimental curve of
specimen L (Fig. 3 of the paper, thirteen markers from 0.5% to 50% gauge
strain) with specimen L.1 on the h020 mesh and no damage.

**Strain rate.** The paper does not state the test rate. With the rate
sensitivity $m_s = 0.12$ the flow stress depends on it: the paper's Table 3
gives 137 / 161 / 163 MPa at 20 / 40 / 50% gauge strain at $10^{-4}$ s$^{-1}$,
155 / 183 / 185 MPa at $3 \times 10^{-4}$ s$^{-1}$ and 178 / 210 / 212 MPa at
$10^{-3}$ s$^{-1}$, against 180 / 233 / 238 MPa in the experiment (RMSE over
the thirteen markers 49.7, 33.2 and 14.3 MPa). A gauge
strain rate of $10^{-3}$ s$^{-1}$, a common quasi-static test rate, is assumed
from here on (provenance: assumed).

**Calibration.** At that rate the paper's set has an RMSE of 14.3 MPa over the
thirteen markers: it follows the experiment up to 25% strain and saturates too
low, at 212 MPa against 239 MPa. One damped Gauss-Newton step on $\rho_0$,
$k_m$ and $y_c$, from one-at-a-time perturbations of a factor 1.2, gives the
carried set (Table 2) with an RMSE of 4.1 MPa (Fig. 2, Table 3). The step
lowers the recovery terms and leaves the initial density at the paper's value
within 3%. The lattice friction is held at the paper's value. A second
candidate with a weaker damping was run and is worse (5.7 MPa); the search
stops there.

**Table 2.** Crystal-plasticity parameters fitted on L.1 (gauge strain rate
$10^{-3}$ s$^{-1}$); all other parameters as Table 1.

| parameter | paper, Table 3 | carried set | ratio |
|---|---|---|---|
| $\rho_0$ (mm$^{-2}$) | $2.35 \times 10^6$ | $2.428 \times 10^6$ | 1.03 |
| $k_m$ | 0.0295 | 0.02555 | 0.87 |
| $y_c$ (mm) | $7.9 \times 10^{-7}$ | $4.91 \times 10^{-7}$ | 0.62 |
| RMSE, 13 markers (MPa) | 14.3 | 4.1 | |

**Table 3.** Engineering stress (MPa) of L.1 without damage against the
experiment.

| gauge strain | experiment | paper's set | carried set |
|---|---|---|---|
| 1% | 85 | 80 | 81 |
| 5% | 98 | 106 | 105 |
| 10% | 129 | 136 | 134 |
| 20% | 180 | 178 | 181 |
| 30% | 212 | 200 | 213 |
| 40% | 233 | 209 | 229 |
| 50% | 238 | 212 | 236 |

**Figure 2.** L.1 without damage: the paper's set and the carried set against
the experimental markers.

![cal](../../figures/ni/cal/cal_curve.png)

The largest remaining error is at yield, where the carried set gives 78 to
81 MPa against 85 MPa below 1% strain. Without damage the engineering stress of
the carried set reaches its maximum, 237 MPa, at 59% gauge strain: the
geometric softening of the section alone does not end the test at the observed
46%. That is the role of the fracture model of §6.

## 5. Mesh check

The single parameter set of §4 and §6 is run unchanged on the h010 mesh of
L.1: 0.1 mm gage cells, 3,168 elements, four times the 792 elements of h020,
six cells per length scale $\ell_0$ instead of three (Table 4, Fig. 3).

**Table 4.** Specimen L.1 with crystal plasticity and fracture on two meshes,
same parameters.

| mesh | elements | failure strain | peak stress (MPa) | fracture location (of the gage) | wall time |
|---|---|---|---|---|---|
| h020 | 792 | 0.4610 | 224.1 | 0.29 | 27 min, 4 threads |
| h010 | 3,168 | 0.4671 | 221.0 | 0.23 | 2 h 26 min, 4 threads |
| difference | | +1.3% | -1.4% | -0.06 | |

**Figure 3.** L.1 on the h020 and h010 meshes: engineering stress and maximum
phase field against the gauge strain.

![mesh](../../figures/ni/pfd/curves_L1_mesh.png)

The failure strain and the peak stress change by 1.3% and 1.4%, inside the
5% acceptance of the study, and the fracture stays in the same quarter of the
gage. h020 is kept for every run that follows.

## 6. Phase-field fracture: one $G_c$ for all specimens

### 6.1 Fracture model

The fracture model is that of (rezwan 2026 predicting): an AT2 phase field
driven by the elastic energy and the full plastic work, with no threshold. The
damage $d$ obeys

$$
\frac{G_c}{\ell_0}\, d - G_c \ell_0 \nabla^2 d = 2 (1-d)\, \mathcal{H}, \qquad
\mathcal{H}(t) = \max_{\tau \le t} \left[ \psi^+(\tau) + \zeta\, W_p(\tau) \right] ,
\tag{13}
$$

with $\psi^+$ the positive part of the elastic energy density, $W_p = \int
\sum_\alpha \tau^\alpha \dot{\gamma}^\alpha \, dt$ the accumulated plastic work
density and $\zeta = 1$. The second Piola-Kirchhoff stress of Eq. (2) is
degraded by $g(d) = (1-d)^2$. The length scale is $\ell_0 = 0.6$ mm, three
gage cells of h020 (§5 checks six). $G_c$ is the one fitted quantity of the
fracture model, and one value serves all fifteen specimens, as in the paper.

**Failure criterion.** Failure is the first load step at which the nodes with
$d > 0.85$ cut the specimen in two (PINNACLE `stop_separation`); the failure
strain is the gauge strain of that step. This is the paper's criterion. The
calibration target is the experimental failure strain of L.1, 0.463 (the
experiment axis of the paper's Fig. 4(b)), accepted within $\pm 0.004$, two
load steps.

### 6.2 Calibration on L.1

The search runs in $\ln G_c$: a first pair brackets the target, then the
interval is closed by regula falsi (`scripts/ni_pfd_calibrate.py`,
`decks/ni/pfd/bunge/cal_L1_h020.i`). Three runs give Table 5 and Fig. 4.

**Table 5.** $G_c$ search on L.1 ($\ell_0 = 0.6$ mm).

| $G_c$ (MPa mm) | failure strain | peak stress (MPa) |
|---|---|---|
| 42.5 | 0.414 | 220 |
| **48.67** | **0.461** | **224** |
| 57.5 | 0.523 | 228 |

**Figure 4.** L.1 with crystal plasticity and fracture: engineering stress and
maximum phase field for three values of $G_c$, against the experiment up to its
failure strain.

![curves](../../figures/ni/pfd/curves_L1_cal.png)

The calibrated value is $G_c = 48.67$ MPa mm, with failure at 0.461 against
0.463. The paper reports $G_c = 125$ MPa mm without its length scale; in AT2 a
failure strain fixes $G_c$ only together with $\ell_0$, and the two frameworks
differ in the crystal-plasticity set and in the rate, so the values are not
expected to coincide. The failure strain is steep in $G_c$ near the target:
1% in $G_c$ moves it by about 0.004.

**Fracture location.** The simulated fracture crosses the gage at 0.29 of its
length, at 76° to the loading axis; the experiment fails at 0.20, on the same
side, on an inclined plane (Fig. 5).

**Figure 5.** L.1 at failure: DIC axial strain of the experiment (top) and the
simulated phase field (bottom). Dashed: the $d = 0.85$ outline. Solid line:
the ridge of maximum $d$ across the width.

![shape](../../figures/ni/pfd/shape_L1_cal_damage.png)

**What the model does not reproduce.** Without a threshold the plastic work
drives damage from the first load step, so $d$ grows over the whole gage
before it localizes (maximum $d$ 0.62 at 35% strain). The degraded stiffness
lowers the curve: the peak stress is 224 MPa against 239 MPa in the
experiment at 45% strain, 6% low, where the undamaged curve of §4 is within
3%. A refit of the crystal-plasticity set with damage active would remove this
offset; it is not done here.

## 7. The fifteen specimens on the finite-element route

### 7.1 Design

The fifteen specimens (five orientation groups, three thickness variants each)
are run with the single set: Table 1 with the three fitted values of Table 2,
$G_c = 48.67$ MPa mm, $\ell_0 = 0.6$ mm, each on its own h020 mesh with its
measured thickness profile (`decks/ni/final/bunge/fem_15_h020.i`). Only L.1
entered the calibration. Every run ends by separation; a run takes 19 to 33
minutes on four threads, 6.7 hours of summed run time for the fifteen.
The results are compared with the paper's experiments only: the failure
strains of its Fig. 4(b) and the fracture locations read from its Fig. 5
(reading uncertainty $\pm 0.08$ of the gage).

### 7.2 Strain to failure

**Table 6.** Experimental and predicted strain to failure and fracture
location (fraction of the gage length) of the fifteen specimens.

| specimen | failure strain, experiment | failure strain, FEM | error (%) | location, experiment | location, FEM | peak stress (MPa) |
|---|---|---|---|---|---|---|
| 9.1 | 0.546 | 0.529 | -3.3 | 0.15 | 0.31 | 199 |
| 9.3 | 0.616 | 0.547 | -11.2 | 0.85 | 0.65 | 204 |
| 9.4 | 0.626 | 0.548 | -12.4 | 0.60 | 0.38 | 200 |
| 10.1 | 0.540 | 0.515 | -4.6 | 0.20 | 0.33 | 173 |
| 10.3 | 0.528 | 0.513 | -3.0 | 0.30 | 0.35 | 174 |
| 10.4 | 0.530 | 0.485 | -8.4 | 0.25 | 0.30 | 171 |
| H.1 | 0.549 | 0.555 | +1.0 | 0.75 | 0.56 | 184 |
| H.3 | 0.566 | 0.558 | -1.3 | 0.75 | 0.53 | 183 |
| H.4 | 0.568 | 0.560 | -1.4 | 0.80 | 0.58 | 183 |
| L.1 | 0.463 | 0.461 | -0.5 | 0.20 | 0.29 | 224 |
| L.3 | 0.481 | 0.466 | -3.3 | 0.65 | 0.67 | 225 |
| L.4 | 0.461 | 0.449 | -2.5 | 0.75 | 0.75 | 220 |
| P.1 | 0.466 | 0.482 | +3.4 | 0.20 | 0.50 | 186 |
| P.3 | 0.378 | 0.443 | +17.4 | 0.75 | 0.79 | 184 |
| P.4 | 0.406 | 0.463 | +14.2 | 0.78 | 0.71 | 184 |

**Figure 6.** Experimental against predicted strain to failure of the fifteen
specimens, with the $\pm 10$% band.

![parity](../../figures/ni/final/fem_parity.png)

Eleven of the fifteen specimens are predicted within 5% and thirteen within
12.5%. Over all fifteen the RMSE is 0.039 in strain, 7.9% relative; the paper
reports 0.064 for its own simulations. The model follows the orientation: the
flow stress ranks L, 9, P, H, 10 from 224 to 173 MPa at the peak, and the
softer crystals strain further, as in the experiment. Groups H and L are
within 3.3%, group 10 within 8.4%. Specimen 9.1 is within 3.3%, while 9.3 and
9.4 fail 11 to 12% early, at 0.55 against 0.62 to 0.63.

The largest errors are P.3 and P.4, predicted 17% and 14% late. The model
gives the three specimens of group P failure strains of 0.44 to 0.48, where
the experiment spreads them from 0.38 to 0.47 with the thickness profile: P.3
and P.4 carry the steepest tapers of the set, $-0.073$ and $-0.054$ mm over
the gage. The model responds to the taper in the right direction and too
weakly.

### 7.3 Fracture location

**Figure 7.** Experimental against predicted fracture location along the gage.

![location](../../figures/ni/final/fem_location.png)

Thirteen of the fifteen fractures are on the experiment's side of the gage
center and six are within the $\pm 0.08$ reading uncertainty; the RMSE is
0.16 of the gage length. The predicted locations are drawn toward
the center compared with the experiment, most clearly in group H, where the
experiment breaks at the fillet (0.75 to 0.80) and the model at 0.53 to 0.58
(Fig. 8).

**Figure 8.** The fifteen specimens at failure, in the order L, P, 10, 9, H. Each
pair shows the DIC axial strain of the experiment at the onset of failure (top;
Fig. 5 of the paper) and the simulated axial strain $E_{xx}$ on the experiment
panel's own color scale (bottom), with the $d = 0.85$ outline (dashed) and the
ridge of maximum $d$ across the width (solid line).

![L.1](../../figures/ni/final/fem_shape_L.1.png)

![L.3](../../figures/ni/final/fem_shape_L.3.png)

![L.4](../../figures/ni/final/fem_shape_L.4.png)

![P.1](../../figures/ni/final/fem_shape_P.1.png)

![P.3](../../figures/ni/final/fem_shape_P.3.png)

![P.4](../../figures/ni/final/fem_shape_P.4.png)

![10.1](../../figures/ni/final/fem_shape_10.1.png)

![10.3](../../figures/ni/final/fem_shape_10.3.png)

![10.4](../../figures/ni/final/fem_shape_10.4.png)

![9.1](../../figures/ni/final/fem_shape_9.1.png)

![9.3](../../figures/ni/final/fem_shape_9.3.png)

![9.4](../../figures/ni/final/fem_shape_9.4.png)

![H.1](../../figures/ni/final/fem_shape_H.1.png)

![H.3](../../figures/ni/final/fem_shape_H.3.png)

![H.4](../../figures/ni/final/fem_shape_H.4.png)

### 7.4 Reading of the result

One crystal-plasticity and phase-field parameter set, fitted on one specimen,
predicts the strain to failure of all five orientation groups with an RMSE of
0.039. The mechanism is direct. The driving force of Eq. (13) is the
accumulated plastic work, close to the area under the stress-strain curve, and
the gage fails when that work reaches a level set by $G_c / \ell_0$. A crystal
with a lower flow stress reaches that level at a larger strain. Groups H and
10, near the [101] pole, carry 173 to 184 MPa at the peak and fail at 0.49 to
0.56; group L, toward [111], carries 224 MPa and fails at 0.45 to 0.47. This
is the ordering of the experiment, and it follows from the orientation alone.

Two limits remain. The fracture location is predicted closer to the gage
center than observed, and the fracture plane is steeper: the simulated
fractures cross the gage at 67° to 90° to the loading axis, where groups 9, L
and P break on planes near 45°. And the sensitivity to the thickness taper is
too weak in group P. Both concern the localization stage, in which the phase
field of this model is spread over the gage before it concentrates (§6.2).

A first version of this study applied the Euler angles with $\varphi_1$ and
$\varphi_2$ exchanged. It simulated groups H and P as hard crystals near
[111], predicted groups 9 and H 19 to 26% early, and read that as two classes
of orientation that one parameter set cannot serve. That reading is withdrawn;
its figures are kept in `archive/figures/ni_exchanged_angles/`.

## 8. PINN+FEM, PINN and surrogate routes

### 8.1 The two network-assisted solver routes

Both routes solve the same discrete problem as §7: the same mesh, the same
total-Lagrangian residual of Eq. (10), the same constitutive update and the
same phase-field solve. They differ from the finite-element route only in how
the mechanics solve of a load step starts and ends.

A network $\mathcal{N}_\theta$ (three hidden layers of 64 units, tanh) with
inputs the reference coordinates and the load factor $\lambda$ corrects the
linear extrapolation of the last two accepted fields,

$$
\mathbf{u}^{\ast}(\mathbf{X}, \lambda) = \mathbf{u}^{\mathrm{lin}}(\mathbf{X}, \lambda)
 + s\, \mathcal{N}_\theta(\mathbf{X}, \lambda), \tag{14}
$$

where $s$ is the largest nodal change of the last accepted increment. The
output layer starts at zero, so the first proposal is the extrapolation. The
network is trained at every load step on the assembled residual itself, with
the reaction as scale,

$$
\mathcal{L}(\theta) = \frac{\lVert \mathbf{R}_f(\mathbf{u}^{\ast}) \rVert^2}
 {\lVert \mathbf{R}_d(\mathbf{u}^{\ast}) \rVert^2}, \tag{15}
$$

where $\mathbf{R}_f$ is the residual on the free degrees of freedom and
$\mathbf{R}_d$ its part on the prescribed ones. The gradient passes through the
finite-element assembly with the transpose of the consistent tangent. Training
uses L-BFGS, five epochs at the first fit and four afterward, and stops as soon
as $\sqrt{\mathcal{L}} \le 0.05$; when the extrapolation already meets this
target no epoch is run.

- **PINN+FEM**: the proposal $\mathbf{u}^{\ast}$ is the starting field of the
  Newton solve, which runs to the tolerance of §3.2
  (`decks/ni/final/bunge/hybrid_base_h020.i`).
- **PINN**: the proposal is accepted as the solution of the mechanics solve,
  with no Newton iteration, when $\lVert \mathbf{R}_f \rVert \le 0.05\,
  \lVert \mathbf{R}_d \rVert$; otherwise Newton corrects it
  (`decks/ni/final/bunge/pinn_base_h020.i`).

The finite-element route of §7 starts every mechanics solve from the same
linear extrapolation without a network.

### 8.2 PINN+FEM on the fifteen specimens

The fifteen specimens were run with the single set of §7
(`decks/ni/final/bunge/hybrid_15_h020.i`). All fifteen end by separation.

| Specimen | Experiment | FEM | PINN+FEM | Location, FEM | Location, PINN+FEM | Newton iterations, FEM | Newton iterations, PINN+FEM | Load steps with training |
|---|---|---|---|---|---|---|---|---|
| L.1 | 0.463 | 0.4610 | 0.4610 | 0.29 | 0.29 | 704 | 704 | 26 of 222 |
| L.3 | 0.481 | 0.4656 | 0.4656 | 0.67 | 0.67 | 679 | 679 | 27 of 214 |
| L.4 | 0.461 | 0.4492 | 0.4492 | 0.75 | 0.75 | 695 | 695 | 25 of 218 |
| P.1 | 0.466 | 0.4817 | 0.4817 | 0.50 | 0.50 | 675 | 675 | 29 of 213 |
| P.3 | 0.378 | 0.4434 | 0.4434 | 0.79 | 0.79 | 637 | 637 | 28 of 198 |
| P.4 | 0.406 | 0.4629 | 0.4629 | 0.71 | 0.71 | 665 | 665 | 26 of 209 |
| 10.1 | 0.540 | 0.5153 | 0.5153 | 0.33 | 0.33 | 727 | 727 | 34 of 231 |
| 10.3 | 0.528 | 0.5128 | 0.5128 | 0.35 | 0.35 | 696 | 696 | 34 of 221 |
| 10.4 | 0.530 | 0.4851 | 0.4851 | 0.30 | 0.30 | 687 | 687 | 34 of 217 |
| 9.1 | 0.546 | 0.5287 | 0.5287 | 0.31 | 0.31 | 615 | 690 | 28 of 241 |
| 9.3 | 0.616 | 0.5471 | 0.5471 | 0.65 | 0.65 | 624 | 635 | 30 of 247 |
| 9.4 | 0.626 | 0.5485 | 0.5485 | 0.38 | 0.38 | 716 | 782 | 23 of 265 |
| H.1 | 0.549 | 0.5547 | 0.5547 | 0.56 | 0.56 | 712 | 798 | 29 of 251 |
| H.3 | 0.566 | 0.5583 | 0.5583 | 0.53 | 0.53 | 733 | 799 | 23 of 255 |
| H.4 | 0.568 | 0.5604 | 0.5604 | 0.58 | 0.58 | 755 | 814 | 19 of 261 |

**Table 7.** Strain to failure, fracture location (fraction of the gage) and
Newton iterations of the fifteen specimens on the finite-element and PINN+FEM
routes. The last column counts the load steps on which the network was trained,
out of the steps with a proposal.

![Parity of the two routes](../../figures/ni/final/routes_parity.png)

**Figure 9.** Experimental against predicted strain to failure of the fifteen
specimens on the finite-element route (left) and the PINN+FEM route (right).

The PINN+FEM route reproduces the finite-element route: the failure strains
agree to six digits and the fracture locations are the same. The error against
the experiment is therefore that of §7, a root-mean-square error of 0.039 in
strain over the fifteen specimens. The per-specimen curves and fields
(`figures/ni/final/hybrid_curve_*.png`, `hybrid_shape_*.png`) are
indistinguishable from those of Figure 8.

The reason is in the last column of Table 7. Over the fifteen runs the linear
extrapolation alone meets the 5 % residual target on 3,048 of 3,463 load steps,
so the network is trained on 415 steps only. On those steps four L-BFGS epochs
lower the residual ratio by 0.1 % at the median and by 20 % at the ninetieth
percentile. The Newton solve then starts from a field that is, for practical
purposes, the extrapolation the finite-element route also starts from. With
0.1 % of bar strain per step the load path is smooth enough that the
extrapolation leaves the network little to correct. No training trial was
refused by the constitutive update in this campaign.

The route does not save Newton iterations. Groups L, P and 10 need the same
iterations on both routes, 6,165 in total. Groups 9 and H need 4,518 on the
PINN+FEM route against 4,155, so the fifteen runs sum to 10,683 against
10,320. The additional iterations fall on load steps without training, where
the first mechanics solve of the step takes two iterations in place of one from
a starting field that differs only in the last digits. The cause was not
investigated; a residual that lies at the convergence tolerance after the first
iteration would show this pattern.

The route costs more than it saves. The fifteen runs take 8.4 hours of summed
run time against 6.7 hours for the finite-element route, on the same host with
seven runs at a time on four threads each; the proposals account for 1.5 hours
of that. The PINN+FEM campaign shared the host with the PINN run of §8.3, so
the difference of 25 % is indicative.

### 8.3 PINN on specimen L.1

The PINN route was run on L.1, the calibration specimen
(`decks/ni/final/bunge/pinn_base_h020.i`).

| Route | Strain to failure | Peak stress (MPa) | Fracture location | Newton iterations | Mechanics solves without Newton | Run time (min) |
|---|---|---|---|---|---|---|
| Experiment | 0.463 | 239 | 0.20 | | | |
| FEM | 0.4610 | 224.1 | 0.29 | 704 | 0 of 462 | 27 |
| PINN+FEM | 0.4610 | 224.1 | 0.29 | 704 | 0 of 462 | 32 |
| PINN | 0.4579 | 224.1 | 0.29 | 336 | 195 of 454 | 31 |

**Table 8.** Specimen L.1 on the three solver routes with the single set. The
run times of the first two rows are from the fifteen-specimen campaigns (seven
runs at a time); the PINN run shared the host with the seven runs of the
PINN+FEM campaign. They are not a like-for-like timing.

![L.1 on the three routes](../../figures/ni/final/routes_curves_L1.png)

**Figure 10.** L.1 on the three solver routes: engineering stress and maximum
phase field against gauge strain, with the experiment up to its failure strain.

![L.1 fields on the three routes](../../figures/ni/final/routes_shape_L1.png)

**Figure 11.** L.1 at failure on the three routes: DIC axial strain of the
experiment (top) and the simulated axial strain with the $d = 0.85$ outline
(dashed) and the ridge line of the phase field.

The PINN route fails at 45.8 % gauge strain against 46.1 % on the
finite-element route, a difference of 0.7 %, with the same peak stress and the
same fracture location. It accepts 195 of 454 mechanics solves without a Newton
iteration and needs 336 Newton iterations in place of 704. The run times of
Table 8 do not show this saving, because the three runs did not have the same
share of the host.

What is accepted must be stated precisely. The accepted fields are those whose
out-of-balance force is within 5 % of the reaction, and on this load path that
is in most cases the linear extrapolation itself: the network was trained on 26
of the 221 load steps. The saving of the PINN route on this specimen is thus
the saving of accepting a verified extrapolation at a looser equilibrium
tolerance than the Newton solve, not the gain of a trained field. The
accumulated effect of that looser tolerance is the 0.7 % shift of the failure
strain. A stricter residual target or longer training would test the network
itself; neither was run here.

### 8.4 Surrogate

**Design.** The surrogate maps a specimen, its orientation and its thickness
profile, to the result of the finite-element route without a solve. Its
training data are thirty synthetic specimens
(`decks/ni/final/bunge/surrogate_design_h020.i`): a Latin hypercube over the
three Euler angles, uniform over the rotation group, and over the thickness,
entered as a mean of 0.320 to 0.400 mm and a linear taper of $-0.08$ to
$+0.07$ mm over the five gage stations. These ranges contain the fifteen
measured profiles. Each synthetic specimen is meshed at h020 and run on the
finite-element route with the single set of §7. All thirty end by separation,
at failure strains of 0.415 to 0.558, peak stresses of 150 to 257 MPa and
fracture locations of 0.23 to 0.81 of the gage. The fifteen real specimens are
never training data.

**Network.** The surrogate is PINNACLE's `parametric_surrogate`
(`decks/ni/final/bunge/surrogate/surrogate_h020.i`). Its inputs are the twelve
Schmid factors of the loading axis, sorted, and the five station thicknesses.
A shared trunk of two layers of 32 units feeds four heads: the engineering
stress-strain curve, learned as the ratio to the material-point curve of the
same orientation with the same crystal-plasticity set, and the failure strain,
the fracture location and the peak stress. Five networks of 3,908 parameters
each, trained from different seeds on the same 24 specimens, form an ensemble;
their spread flags a query above 0.1. Six of the thirty specimens are held
out.

**Held-out synthetic specimens.** Figure 12 shows the six held-out curves and
Figure 13 the three scalars. The curve error is 5.8% on average and 9.6% at
most. The failure strain has an RMS error of 2.2%, with five specimens within
2.3% and one at 4.4%; the peak stress has 3.2%, and the fracture location 0.040
of the gage. The training specimens are reproduced within 0.2% on the scalar
heads, so with 24 specimens the networks interpolate, and the held-out error
is the measure of the surrogate.

![Held-out curves](../../figures/ni/final/surrogate_holdout_curves.png)

**Figure 12.** The six held-out synthetic specimens: engineering stress against
gauge strain from the finite-element run and from the surrogate, with the
material-point curve the surrogate scales. Dotted lines mark the two failure
strains.

![Held-out scalars](../../figures/ni/final/surrogate_holdout_parity.png)

**Figure 13.** Failure strain, fracture location and peak stress of the thirty
synthetic specimens: finite-element run against surrogate.

**The fifteen real specimens.** Table 9 and Figure 14 give the surrogate's
predictions for the real specimens against the finite-element route and the
experiment.

| specimen | failure strain, experiment | FEM | surrogate | surrogate against FEM (%) | location, experiment | FEM | surrogate | peak stress, FEM (MPa) | surrogate (MPa) | ensemble spread |
|---|---|---|---|---|---|---|---|---|---|---|
| 9.1 | 0.546 | 0.529 | 0.522 | -1.2 | 0.15 | 0.31 | 0.32 | 199 | 209 | 0.027 |
| 9.3 | 0.616 | 0.547 | 0.520 | -4.9 | 0.85 | 0.65 | 0.61 | 204 | 211 | 0.016 |
| 9.4 | 0.626 | 0.548 | 0.544 | -0.8 | 0.60 | 0.38 | 0.52 | 200 | 202 | 0.031 |
| 10.1 | 0.540 | 0.515 | 0.533 | +3.4 | 0.20 | 0.33 | 0.33 | 173 | 165 | 0.027 |
| 10.3 | 0.528 | 0.513 | 0.536 | +4.5 | 0.30 | 0.35 | 0.33 | 174 | 168 | 0.038 |
| 10.4 | 0.530 | 0.485 | 0.525 | +8.2 | 0.25 | 0.30 | 0.20 | 171 | 161 | 0.043 |
| H.1 | 0.549 | 0.555 | 0.549 | -1.1 | 0.75 | 0.56 | 0.64 | 184 | 190 | 0.021 |
| H.3 | 0.566 | 0.558 | 0.551 | -1.2 | 0.75 | 0.53 | 0.58 | 183 | 191 | 0.025 |
| H.4 | 0.568 | 0.560 | 0.556 | -0.8 | 0.80 | 0.58 | 0.61 | 183 | 189 | 0.022 |
| L.1 | 0.463 | 0.461 | 0.454 | -1.6 | 0.20 | 0.29 | 0.36 | 224 | 232 | 0.036 |
| L.3 | 0.481 | 0.466 | 0.453 | -2.8 | 0.65 | 0.67 | 0.63 | 225 | 223 | 0.046 |
| L.4 | 0.461 | 0.449 | 0.443 | -1.4 | 0.75 | 0.75 | 0.79 | 220 | 228 | 0.039 |
| P.1 | 0.466 | 0.482 | 0.485 | +0.6 | 0.20 | 0.50 | 0.51 | 186 | 179 | 0.066 |
| P.3 | 0.378 | 0.443 | 0.439 | -1.0 | 0.75 | 0.79 | 0.74 | 184 | 180 | 0.038 |
| P.4 | 0.406 | 0.463 | 0.466 | +0.6 | 0.78 | 0.71 | 0.69 | 184 | 178 | 0.047 |

**Table 9.** The fifteen real specimens on the surrogate, against the
finite-element route it was trained to reproduce and against the experiment.

![Real specimens on the surrogate](../../figures/ni/final/surrogate_real.png)

**Figure 14.** The fifteen real specimens: surrogate against the
finite-element failure strain (left), against the experimental failure strain
(middle) and against the experimental fracture location (right).

Against the finite-element route the surrogate's failure strain has an RMS
difference of 3.1%, 0.015 in strain, somewhat above the held-out error of
2.2%. Fourteen specimens are within 5% and ten within 2%. The largest
differences are 10.4 at $+8.2$%, 9.3 at $-4.9$% and 10.3 at $+4.5$%. The peak
stress differs by 3.6% RMS and the fracture location by 0.06 of the gage, with
9.4 the largest at 0.15.

All fifteen queries lie inside the designed ranges, and the ensemble spread is
0.02 to 0.07, below the flag level of 0.1. The spread does not single out 10.4:
its value of 0.043 is exceeded by three specimens that the surrogate predicts
within 3%.

Against the experiment the surrogate inherits the finite-element result and
adds little of its own: the RMSE in failure strain is 0.042 over the fifteen
specimens against 0.039 for the finite-element route. P.3 and P.4 are
predicted late and 9.3 and 9.4 early by the surrogate, as by the model it
learned.

The surrogate does not write a curve for a queried specimen, only the curve's
peak and the three scalars, so no stress-strain curve of a real specimen is
shown for this route.

### 8.5 The routes compared

| Route | Scope | Failure strain against FEM | Time |
|---|---|---|---|
| FEM | 15 specimens | reference | 6.7 h summed, 19 to 33 min per specimen |
| PINN+FEM | 15 specimens | equal to six digits | 8.4 h summed |
| PINN | L.1 | $-0.7$% | 336 Newton iterations against 704; run times not comparable |
| Surrogate | 15 specimens, none in training | 3.1% RMS, 8.2% at most | 0.9 s per specimen, after 15.0 h of training runs and 92 s of training |

**Table 10.** The four routes with the single parameter set. Run times are
summed over the runs, each on four threads of the same host; the surrogate was
trained and queried on eight threads.

The thirty training runs take 2 h 22 min of wall time with seven at a time,
and 15.0 hours summed, 20 to 44 minutes each. The surrogate deck runs in 92 s:
28 s for the baseline curves, 48 s for the training of the five networks and
the rest, 0.9 s per specimen, for the fifteen queries, most of it the
material-point curve of the queried orientation. The surrogate therefore
costs about as much as 33 finite-element runs before its first answer, and
about 1,700 times less than a run for each answer after that.

On this problem the two solver routes do not improve on the finite-element
route. The load path is smooth at 0.1% bar strain per step, the linear
extrapolation meets the residual target on most steps, and the network has
little left to correct; the network settings were kept as they are. The PINN
route halves the Newton iterations on L.1 at a small cost in accuracy, and its
saving is that of a looser equilibrium tolerance. The surrogate is the one
route that changes the cost of a prediction, by three orders of magnitude, and
its accuracy of 2 to 3% in failure strain is set by the thirty training
specimens.

Memory use was not recorded on any route.
