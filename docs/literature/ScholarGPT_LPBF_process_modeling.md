# **Research report: Process modeling of laser powder bed fusion (LPBF)**

## **Multiphysics, multiscale FEM–FFT frameworks using JAX, JAX-FEM, JAX-CP, MOOSE, and SIERRA**

## **1. Executive summary**

Laser powder bed fusion (LPBF) is a complex additive manufacturing (AM) process involving highly transient heat transfer, melting and solidification, fluid dynamics, phase transformations, microstructure evolution, thermomechanical deformation, crystallographic plasticity, and residual stress development.

A comprehensive computational framework for LPBF should ideally establish a **process–structure–property–performance (PSPP) relationship**, connecting manufacturing parameters to the final mechanical properties of the manufactured component.

Five computational platforms are particularly relevant to developing such a framework:

| Framework         | Primary capability                                                  | Potential application in LPBF                                |
| ----------------- | ------------------------------------------------------------------- | ------------------------------------------------------------ |
| JAX               | Automatic differentiation, GPU computing, numerical transformations | Differentiable multiphysics and inverse modeling             |
| JAX-FEM           | Differentiable finite element analysis                              | Thermal and thermomechanical process simulations             |
| JAX-CP            | Crystal plasticity modeling using JAX                               | Grain-scale deformation and mechanical properties            |
| MOOSE             | Parallel multiphysics finite element analysis                       | Thermal, phase-field, microstructure and mechanical modeling |
| SIERRA            | Large-scale computational multiphysics                              | High-fidelity thermal and structural process simulation      |
| FFT-based solvers | Spectral micromechanical homogenization                             | Polycrystalline deformation and residual stress analysis     |

The most promising research direction is a **differentiable, multiscale FEM–FFT computational framework** in which finite element simulations predict the manufacturing thermal history, microstructure models predict grain morphology and crystallographic texture, and crystal plasticity FEM or FFT solvers predict local and effective mechanical properties.

Such a framework could also support inverse identification of material parameters, process optimization, uncertainty quantification, and physics-based machine learning.

### **Proposed multiscale LPBF modeling architecture**

Proposed research architecture, rather than an existing unified software package.

The following report separates governing physical equations, numerical methods, software architectures, published research, and proposed simulation configurations. It also distinguishes capabilities that are established in the literature from those that would require original research and implementation.

## **2. Fundamental physics of laser powder bed fusion**

LPBF involves physical phenomena operating over substantially different spatial and temporal scales. This makes a single, fully coupled, high-resolution model computationally expensive, particularly for component-scale simulations.

[Agnew Group | University of Virginia School of Engineering and Applied Science](https://images.openai.com/static-rsc-4/Q_11LRiGb1n32l60OW32bauGcwBP9Jg5xZY1v_Kw0iP2QRqV22lLDJiTpRMbQYPLySYPVkxGDqUaK-YYMGP9OGbwC6s1Njk_5fbRGU1Ujf59mzE9APYfHlIAj-Ej5mpLftmfpgYdrgZ5R3TX4_tztbD8u3jKDu5WgFXBI_lBh7s?purpose=inline)

[Computational Assessment of Thermokinetics and Associated Microstructural Evolution in Laser Powder Bed Fusion Manufacturing of Ti6Al4V Alloy | Scientific Reports](https://images.openai.com/static-rsc-4/XYFNbvhfxE6UavpzKdxgHuLqdOCsbsaCcmrncJ7PkuGA8gyPHJVxrUZajl3ycpyxXXNG_EPunop6678YLhuda27hFKAN2CJo15cUeXYwWCym5vZlQx4vPP4GzUEk_eG2HJnUWzG5l8-QW1I9UH-Le1XYnF_VgmNT0LpU9xcTyp4?purpose=inline)

[Anisotropy of Additively Manufactured Metallic Materials](https://images.openai.com/static-rsc-4/172oMnkqUYxH-XDV3svWIxdLT5gOgcJlN9Y1vvkVtZiTYBw2iVxu5ZzKWIrmaoBzapl-0V-9PPYIjAeTtmcP9r5sN6HuG58eCacx8HFd0rjn6vTtGynBJvqUbMzdHjcQOaAG0UFnLlUaH18ijk2dRmZkCtlLa0Jstnt4KYgx6BA?purpose=inline)

14

### **2.1 Multiscale modeling hierarchy**

| Modeling scale       | Characteristic length | Physical phenomena                  | Numerical methods   |
| -------------------- | --------------------- | ----------------------------------- | ------------------- |
| Powder scale         | 10–100 µm             | Particle packing, laser absorption  | DEM, ray tracing    |
| Melt-pool scale      | 50–1,000 µm           | Melting, fluid flow, evaporation    | FEM, FVM, CFD       |
| Microstructure scale | 0.1–500 µm            | Nucleation, grain growth, texture   | Phase field, CA     |
| Grain scale          | 1–500 µm              | Slip, hardening, grain interactions | CP-FEM, CP-FFT      |
| Component scale      | 1–1,000 mm            | Residual stress, distortion         | FEM                 |
| Structural scale     | mm–m                  | Mechanical performance              | FEM, homogenization |

These length scales are indicative and vary with alloy, processing conditions, and model resolution.

### **2.2 Governing thermal equation**

The transient heat conduction equation is the foundation of most LPBF process models:

\\[ \rho(T)c_p(T)\frac{\partial T}{\partial t} = \nabla\cdot\left(k(T)\nabla T\right)+Q \\]

where:

| Symbol     | Description            |
| ---------- | ---------------------- |
| \\(T\\)    | Temperature            |
| \\(\rho\\) | Material density       |
| \\(c_p\\)  | Specific heat capacity |
| \\(k\\)    | Thermal conductivity   |
| \\(Q\\)    | Volumetric heat source |
| \\(t\\)    | Time                   |

To account for melting and solidification, an enthalpy formulation is preferable:

\\[ \rho\frac{\partial H}{\partial t} = \nabla\cdot(k\nabla T)+Q \\]

\\[ H(T)=\int\_{T_0}^{T}c_p(\theta)d\theta+L_f f_l(T) \\]

Here, \\(L_f\\) is the latent heat of fusion and \\(f_l\\) is the liquid fraction.

A simple solid–liquid transition model is:

\\[ f_l(T)= \begin{cases} 0 & T\leq T_s\\\ \dfrac{T-T_s}{T_l-T_s} & T_s\<T\<T_l\\\ 1 & T\geq T_l \end{cases} \\]

where \\(T_s\\) and \\(T_l\\) are the solidus and liquidus temperatures.

The equivalent apparent heat capacity becomes:

\\[ c_p^{\mathrm{eff}}=c_p+L_f\frac{df_l}{dT} \\]

This formulation is suitable for implementation in JAX-FEM, MOOSE, and SIERRA thermal solvers.

### **2.3 Moving Gaussian laser heat source**

A commonly used surface heat source is:

\\[ q\_{\mathrm{laser}}(x,y,t)= \frac{2\eta P}{\pi r_b^2} \exp\left[ -\frac{2((x-x_l)^2+(y-y_l)^2)}{r_b^2} \right] \\]

where:

- \\(P\\): laser power.
- \\(\eta\\): effective laser absorptivity.
- \\(r_b\\): laser beam radius under the \\(1/e^2\\) convention.
- \\(x_l,y_l\\): instantaneous laser coordinates.
- \\(v\\): laser scanning velocity.

For a straight scan in the \\(x\\)-direction:

\\[ x_l(t)=x_0+vt,\qquad y_l(t)=y_0 \\]

### **Interactive Gaussian laser heat-source model**

Laser power

280 W

Beam radius

35 µm

Absorptivity

0.35

Scan velocity

1050 mm/s

Peak heat flux

# **50.93**

GW/m²

Absorbed power

# **98.0 W**

Absorbed line energy

## **0.0933 J/mm**

Beam interaction time

## **66.7 µs**

Analytical incident heat-flux distribution. This is not a solved temperature or melt-pool field.

### **2.4 Thermal boundary conditions**

The exposed surfaces lose heat through convection and radiation:

\\[ -k\nabla T\cdot\mathbf n = h_c(T-T\_\infty) + \epsilon\sigma\_{\mathrm{SB}} (T^4-T\_\infty^4) -q\_{\mathrm{laser}} \\]

where \\(h_c\\) is the convection coefficient, \\(\epsilon\\) is emissivity, and \\(\sigma\_{\mathrm{SB}}\\) is the Stefan–Boltzmann constant.

The initial condition is:

\\[ T(\mathbf x,0)=T\_{\mathrm{preheat}} \\]

At the baseplate, a prescribed temperature, thermal contact conductance, or an explicitly meshed substrate can be employed.

For research-grade simulations, explicitly modeling the substrate is recommended when cooling rates and residual stresses are important.

### **2.5 Volumetric energy density**

A frequently used processing parameter is:

\\[ E_v=\frac{P}{vht_l} \\]

where \\(h\\) is hatch spacing and \\(t_l\\) is layer thickness.

However, volumetric energy density alone cannot uniquely characterize melt-pool behavior. Different combinations of power and velocity can produce the same \\(E_v\\) while resulting in different thermal gradients, peak temperatures, melt-pool geometries, and porosity.

## **3. High-fidelity melt-pool modeling**

A conduction-only FEM model predicts temperature distributions efficiently but does not explicitly capture melt-pool circulation, recoil pressure, keyhole formation, spattering, or free-surface deformation.

These phenomena require coupled fluid mechanics and heat transfer.

### **3.1 Conservation of mass**

For an incompressible liquid:

\\[ \nabla\cdot\mathbf u=0 \\]

### **3.2 Conservation of momentum**

\\[ \rho\left( \frac{\partial\mathbf u}{\partial t} +\mathbf u\cdot\nabla\mathbf u \right) = -\nabla p +\nabla\cdot\boldsymbol{\tau} +\rho\mathbf g +\mathbf f\_\sigma +\mathbf f_d \\]

The terms represent pressure, viscous stresses, gravity, surface tension, and solidification-induced damping.

An enthalpy–porosity approximation for the mushy region is:

\\[ \mathbf f_d= -C_m\frac{(1-f_l)^2}{f_l^3+\varepsilon} \mathbf u \\]

where \\(C_m\\) is a mushy-zone coefficient.

### **3.3 Energy conservation including advection**

\\[ \rho\left( \frac{\partial H}{\partial t} +\mathbf u\cdot\nabla H \right) = \nabla\cdot(k\nabla T)+Q \\]

### **3.4 Marangoni convection**

Surface-tension gradients generate shear stresses along the melt-pool surface:

\\[ \boldsymbol{\tau}\_M= \frac{d\gamma}{dT}\nabla_s T \\]

where \\(\gamma\\) is surface tension.

[Systematic approach to process parameter optimization for laser powder bed fusion of low-alloy steel based on melting modes | The International Journal of Advanced Manufacturing Technology | Springer Nature Link](https://images.openai.com/static-rsc-4/NzmkYrbwhBu3-FcvrWaY4AOAqmVPuqY3MeReS5V-Q6ljNlEHnMP5ho_K08kvfHlEjUymlkLas-i4irgtj2M1COLeLZ3pAJZE05ltIf3ECi4fNPcHfCIhj409BtVt5MqusckO2dkNAsUlAA06vfBB7N8PL36hi2hE2_qNaKg0NaY?purpose=inline)

[Absorbance measurement for in situ process regime identification in laser processing | The International Journal of Advanced Manufacturing Technology | Springer Nature Link](https://images.openai.com/static-rsc-4/XWLYSAOgd5upkrwSg3MpFH7KnT_3fv-wQ_9jyGhmAwNOvYpV0gRgZAo5V5F_ok5paHEzTQqkqxlcao20AKQ7o0nQi5dJo_n9OQGZPe2n_cEKA2nAbtETGHl730Q7MatiELC4JHF9IkdHnnr2W6dhLpPZT10lERz4LRC-AfywIP0?purpose=inline)

[Mechanism of keyhole pore formation in metal additive manufacturing | npj Computational Materials](https://images.openai.com/static-rsc-4/AID5zy4uThndVb1KcJBfLLgRa9IGMnnRsuqjDiI30kNSC1IcbcQ-sqZ062ZSLapAyJjlIAQBFxeoNTgx73zskXSruw40tgECswEiDARnGxGwnT_dX5_y_-ggCD8H_zkdH9C1tHs0doP8HtQm8vGAuOZEm6L2FOBiKpnq9seKSKY?purpose=inline)

17

These effects can substantially alter melt-pool geometry and solidification conditions.

### **3.5 Numerical methods for melt-pool modeling**

| Method                   | Advantages                        | Limitations                             |
| ------------------------ | --------------------------------- | --------------------------------------- |
| Conduction FEM           | Efficient, scalable               | No explicit fluid flow                  |
| CFD with volume of fluid | Free-surface and keyhole tracking | High computational cost                 |
| Level-set method         | Accurate interface representation | Mass conservation challenges            |
| Enthalpy–porosity        | Convenient phase-change treatment | Mushy-zone parameter sensitivity        |
| Lattice Boltzmann        | Local operations, GPU suitability | Calibration and multiphysics complexity |
| DEM–CFD                  | Explicit powder particles         | Expensive for multilayer builds         |

For the proposed FEM–FFT research framework, a calibrated thermal conduction model is a practical starting point. High-fidelity melt-pool CFD can subsequently provide reference solutions for improved heat-source calibration.

## **4. Solidification and microstructure evolution**

The predicted temperature field provides the primary input for grain-scale microstructure modeling.

### **4.1 Solidification parameters**

The thermal gradient is:

\\[ G=\\|\nabla T\\| \\]

The local normal velocity of a solidification isotherm can be estimated by:

\\[ R=\frac{|\partial T/\partial t|}{\\|\nabla T\\|} \\]

evaluated at the appropriate solidification interface.

The cooling rate is then:

\\[ |\dot T|\approx GR \\]

The ratio \\(G/R\\) is associated with morphological stability and solidification mode, while \\(GR\\) influences characteristic microstructural length scales.

### **Thermal conditions and grain morphology**

[Checking your browser - reCAPTCHA](https://images.openai.com/static-rsc-4/v41lssQ8QqCDOAxiosj3lHydC8l3gCpOCyVuAkcvIkoBrg-UR6hBBgYX8CMy9qNwFH3IHX1pbwJOPkw5DR2GBVSOrsfyb5rI8w1f5HAtcWc2sQdPXbZ0q4G6ROW1vAb7zzJAyg9iW_CGUQxfLmTpQ6KzgwwMwqbRBy-M0LpqZk0?purpose=inline)

Cellular

[Helical milling hole finishing of Inconel 718 parts manufactured by laser powder bed fusion | Progress in Additive Manufacturing | Springer Nature Link](https://images.openai.com/static-rsc-4/i0g4vgipbF22qdm5EoCjZ1z_29ZW78EpQ7ot2np2vpoU0mOSYUx46Eh4zn7rXBLdwnBpoPrtl0B5oIgwDW2OPlN90SaZPk84jPtutTqkJBIBTBIMxis9JEMbIlGCi4lh5wd0fcvsfscT3b_aFRiy1gMiSyY60q3sQ0trSzu6Xx0?purpose=inline)

Columnar

[Influence of deformation banding instabilities on small scale yielding of a Mg–Nd alloy revealed by in-situ digital image correlation | Scientific Reports](https://images.openai.com/static-rsc-4/INpDgzqJgCyiPGwdslKzXPQLnH5OtrkPqPsqv5wQutWSM9-85IQU784pq8M_2aatRIkVSY01tY9ZawQ7CuA_p6JK9l3rd-yUiEoXIfNUqfmHrH7ECTkjrcpQPNeQs7_AH34ebV2IiyjRO7pwGlFQTfWAkB1r6HJIeI52RriaC4Q?purpose=inline)

Equiaxed

Qualitative solidification map. Actual morphology transitions depend on composition, nucleation undercooling, anisotropy, and alloy-specific kinetics.

### **4.2 Phase-field formulation**

A general free-energy functional is:

\\[ \mathcal F= \int\_\Omega \left[ f(\phi,c,T) +\frac{\kappa\_\phi}{2}|\nabla\phi|^2 +\frac{\kappa_c}{2}|\nabla c|^2 \right]dV \\]

where \\(\phi\\) is a phase-field order parameter and \\(c\\) is solute concentration.

The nonconserved order parameter evolves according to the Allen–Cahn equation:

\\[ \frac{\partial\phi}{\partial t} = -L\_\phi\frac{\delta\mathcal F}{\delta\phi} \\]

The conserved concentration field evolves according to the Cahn–Hilliard equation:

\\[ \frac{\partial c}{\partial t} = \nabla\cdot \left( M_c\nabla\frac{\delta\mathcal F}{\delta c} \right) \\]

Expanded:

\\[ \frac{\partial\phi}{\partial t} = -L\_\phi \left( \frac{\partial f}{\partial\phi} -\kappa\_\phi\nabla^2\phi \right) \\]

\\[ \frac{\partial c}{\partial t} = \nabla\cdot \left[ M_c\nabla \left( \frac{\partial f}{\partial c} -\kappa_c\nabla^2c \right) \right] \\]

The phase-field model can incorporate crystallographic anisotropy, grain-boundary energy, heterogeneous nucleation, and thermal coupling.

Quantitative alloy solidification generally requires additional thermodynamic and kinetic treatment, such as CALPHAD-based free energies and appropriate interface corrections.

### **4.3 Cellular automata method**

Cellular automata (CA) is often computationally more practical than phase field for large LPBF grain-structure simulations.

The computational domain is discretized into cells, each storing:

- Grain identification and crystallographic orientation.
- Solid, liquid, or interface state.
- Local temperature and undercooling.
- Grain-growth front position.

A general growth law is:

\\[ v_g=a(\Delta T)^n \\]

\\[ \Delta T=T_l-T \\]

where \\(a\\) and \\(n\\) are alloy-dependent kinetic parameters.

The nucleation distribution can be represented by:

\\[ \frac{dN}{d\Delta T} = \frac{N\_{\max}}{\sqrt{2\pi}\sigma\_{\Delta T}} \exp\left[ -\frac{(\Delta T-\mu\_{\Delta T})^2} {2\sigma\_{\Delta T}^2} \right] \\]

The grain-growth front is advanced based on the local undercooling and crystallographic preferred growth directions.

For large-scale AM simulations, **FEM–CA–FFT** is an especially attractive computational architecture:

\\[ \boxed{ \text{Thermal FEM} \rightarrow \text{Cellular automata} \rightarrow \text{Crystal plasticity FFT} } \\]

It offers a reasonable compromise between computational efficiency, physical fidelity, and microstructure-sensitive mechanical predictions.

## **5. Thermomechanical process modeling**

After obtaining the temperature history, thermally induced deformation can be predicted using the mechanical equilibrium equation:

\\[ \nabla\cdot\boldsymbol\sigma+\mathbf b=0 \\]

The small-strain tensor is:

\\[ \boldsymbol\varepsilon= \frac12 \left( \nabla\mathbf u+\nabla\mathbf u^T \right) \\]

The total strain is decomposed into:

\\[ \boldsymbol\varepsilon= \boldsymbol\varepsilon^e+ \boldsymbol\varepsilon^p+ \boldsymbol\varepsilon^{th}+ \boldsymbol\varepsilon^{tr} \\]

where the last term represents transformation strain, when relevant.

The thermal strain is:

\\[ \boldsymbol\varepsilon^{th} = \left[ \int\_{T\_{\mathrm{ref}}}^{T}\alpha(\theta)d\theta \right]\mathbf I \\]

The stress–strain relationship is:

\\[ \boldsymbol\sigma= \mathbb C(T): \left( \boldsymbol\varepsilon- \boldsymbol\varepsilon^p- \boldsymbol\varepsilon^{th}- \boldsymbol\varepsilon^{tr} \right) \\]

### **5.1 J2 plasticity**

A conventional von Mises yield criterion is:

\\[ f= \sqrt{\frac32\mathbf s:\mathbf s} -\sigma_y(T,\bar\varepsilon^p) \\]

with the deviatoric stress:

\\[ \mathbf s= \boldsymbol\sigma- \frac13\operatorname{tr}(\boldsymbol\sigma)\mathbf I \\]

The associative plastic flow rule is:

\\[ \dot{\boldsymbol\varepsilon}^p = \dot\lambda \frac{\partial f}{\partial\boldsymbol\sigma} \\]

A simple isotropic hardening model is:

\\[ \sigma_y= \sigma\_{y0}(T)+H(T)\bar\varepsilon^p \\]

In an LPBF simulation, elastic modulus, yield stress, hardening, and thermal expansion should generally be temperature dependent.

### **5.2 Material activation and residual stress**

Additive manufacturing models must account for the progressive addition of material.

| Approach                | Description                                                                |
| ----------------------- | -------------------------------------------------------------------------- |
| Element birth and death | Elements become mechanically and thermally active during deposition        |
| Quiet element method    | Initially inactive elements have very low stiffness or conductivity        |
| Layer activation        | Entire layers are activated sequentially                                   |
| Track activation        | Individual scan tracks are activated progressively                         |
| Inherent strain         | Prescribed equivalent irreversible strains replace detailed thermal cycles |

For high-fidelity research, melting-induced stress relaxation must also be considered. Applying a conventional solid-state plasticity model without treating the liquid-state loss of shear stiffness can produce unrealistic residual stresses.

# **Part II. Computational frameworks and software development**

## **6. JAX: Differentiable scientific computing**

[GitHub - jax-ml/jax: Composable transformations of Python+NumPy programs: differentiate, vectorize, JIT to GPU/TPU, and more · GitHub](https://images.openai.com/static-rsc-4/LXPJ9qVWQK-YC7V4a4xlpPBoLR8AlFXu4ezE996ZLixoHpPYLCyytRmhS1UB9KmacUWqG884d4fR8rP9n4BuEh8oNE-L_NGjLys9RtXxZgLpKjoCpuRvTOYA7pWbv6p3zwa0GR9erd4JgcIbtxaCc3-MCWXNvvbdnNba_boKc5o?purpose=inline)

**JAX**

Automatic differentiation and accelerated numerical computing.

[Official documentation](https://docs.jax.dev/?utm_source=chatgpt.com) · [GitHub repository](https://github.com/jax-ml/jax?utm_source=chatgpt.com)

### **6.1 Development and motivation**

JAX was developed as a numerical computing framework combining a NumPy-like interface with automatic differentiation and accelerated execution through XLA compilation.

Unlike conventional FEM software, JAX is not itself a finite element solver. Instead, it provides the computational infrastructure required to develop differentiable finite element, finite difference, spectral, and constitutive modeling frameworks.

This makes JAX particularly valuable for computational materials science, where forward simulations are increasingly combined with parameter identification and inverse design.

### **6.2 Core computational methods**

| JAX operation    | Function                       | LPBF application               |
| ---------------- | ------------------------------ | ------------------------------ |
| `jax.grad`       | Scalar-function gradient       | Process optimization           |
| `jax.jacfwd`     | Forward-mode Jacobian          | Constitutive tangent           |
| `jax.jacrev`     | Reverse-mode Jacobian          | Residual differentiation       |
| `jax.jit`        | Just-in-time compilation       | Accelerated FEM kernels        |
| `jax.vmap`       | Automatic vectorization        | Element and voxel calculations |
| `jax.lax.scan`   | Compiled sequential operations | Time integration               |
| `jax.numpy.fft`  | Fast Fourier transforms        | Spectral crystal plasticity    |
| `jax.custom_vjp` | Custom reverse differentiation | Implicit nonlinear solvers     |

### **6.3 Automatic differentiation of governing equations**

Consider a nonlinear finite element system:

\\[ \mathbf R(\mathbf u,\boldsymbol\theta)=0 \\]

where \\(\mathbf u\\) is the solution and \\(\boldsymbol\theta\\) is a set of material or process parameters.

The Newton iteration is:

\\[ \mathbf K^{(k)}\Delta\mathbf u^{(k)} =-\mathbf R^{(k)} \\]

\\[ \mathbf u^{(k+1)} = \mathbf u^{(k)}+\Delta\mathbf u^{(k)} \\]

The tangent matrix is:

\\[ \mathbf K= \frac{\partial\mathbf R}{\partial\mathbf u} \\]

Conventional FEM development often requires analytical derivation and manual implementation of this matrix. JAX can compute element-level Jacobians using automatic differentiation.

### **6.4 Differentiable inverse problems**

Suppose the optimization objective is:

\\[ J(\boldsymbol\theta)= \frac12 \\|\mathbf T\_{\mathrm{sim}}- \mathbf T\_{\mathrm{exp}}\\|^2 +\lambda\mathcal R(\boldsymbol\theta) \\]

The parameters could include:

\\[ \boldsymbol\theta= [P,v,\eta,k,c_p,r_b] \\]

For an implicitly defined numerical solution:

\\[ \frac{d\mathbf u}{d\boldsymbol\theta} = -\left( \frac{\partial\mathbf R}{\partial\mathbf u} \right)^{-1} \frac{\partial\mathbf R}{\partial\boldsymbol\theta} \\]

Introducing the adjoint variable \\(\boldsymbol\lambda\\):

\\[ \left( \frac{\partial\mathbf R}{\partial\mathbf u} \right)^T\boldsymbol\lambda = \left( \frac{\partial J}{\partial\mathbf u} \right)^T \\]

yields:

\\[ \frac{dJ}{d\boldsymbol\theta} = \frac{\partial J}{\partial\boldsymbol\theta} - \boldsymbol\lambda^T \frac{\partial\mathbf R}{\partial\boldsymbol\theta} \\]

Adjoint differentiation is particularly useful when the number of design variables is large compared with the number of objectives.

For transient LPBF models, the adjoint must also account for the dependence of every time increment on previous thermal states and constitutive history variables.

### **6.5 Example JAX implementation of the LPBF heat source**

```
import jaximport jax.numpy as jnpjax.config.update("jax_enable_x64", True)@jax.jitdef laser_flux(x, y, t, power, speed,               radius, absorptivity):    xc = speed * t    r2 = (x - xc)**2 + y**2    return (        2.0 * absorptivity * power        / (jnp.pi * radius**2)        * jnp.exp(-2.0 * r2 / radius**2)    )dq_dpower = jax.grad(laser_flux, argnums=3)dq_dspeed = jax.grad(laser_flux, argnums=4)dq_dradius = jax.grad(laser_flux, argnums=5)
```

This code illustrates differentiability with respect to the laser parameters. It does not constitute a complete thermal solver.

### **6.6 Research advantages and limitations**

JAX provides an excellent environment for developing new constitutive models, differentiable numerical solvers, GPU-compatible material calculations, and inverse optimization methods.

Its principal limitations for LPBF are the need to develop substantial multiphysics infrastructure, manage GPU memory, implement material activation, and carefully treat nondifferentiable events such as nucleation, phase switching, and changing mesh topology.

## **7. JAX-FEM: Differentiable finite element modeling**

[Thermal mechanical control — JAX-FEM 0.0.10 documentation](https://images.openai.com/static-rsc-4/es3Mc7bZPYhUM8iE_jQTOO0geYbzzalcxbDafqsY0T3MrHxJqIi1mDrlSYMkRlu_kOD65n8OwHUD3ZPaE8xUWWv2gZ1RBVz60CerWQuijLGoULdTBzRN0YCcSc3b1-asVodutxHvnrkyvCSy22Eix1jS02FdcavGs6ba1CbRo8U?purpose=inline)

**JAX-FEM**

GPU-accelerated differentiable three-dimensional finite element analysis.

[GitHub](https://github.com/deepmodeling/jax-fem?utm_source=chatgpt.com) · [Documentation](https://deepmodeling.github.io/jax-fem/?utm_source=chatgpt.com)

### **7.1 Development**

JAX-FEM was developed by researchers including Tianju Xue and collaborators associated with computational mechanics and advanced manufacturing research.

A foundational publication is:

*Xue et al. (2023). JAX-FEM: A differentiable GPU-accelerated 3D finite element solver for automatic inverse design and mechanistic data science. Computer Physics Communications, 291, 108802.*

[DOI: 10.1016/j.cpc.2023.108802](https://doi.org/10.1016/j.cpc.2023.108802?utm_source=chatgpt.com)

The framework integrates finite element discretization, automatic differentiation, vectorized quadrature-point operations, and nonlinear equation solving.

### **7.2 Finite element formulation for LPBF**

Starting from the enthalpy heat equation:

\\[ \rho\dot H-\nabla\cdot(k\nabla T)-Q=0 \\]

Multiplying by a test function \\(w\\) and integrating over the domain gives:

\\[ \begin{aligned} R_T={}& \int\_\Omega w\rho\dot H\\,d\Omega\\\ &+\int\_\Omega\nabla w\cdot k\nabla T\\,d\Omega\\\ &-\int\_\Omega wQ\\,d\Omega\\\ &+\int\_{\Gamma_q}wq\_{\mathrm{out}}\\,d\Gamma =0 \end{aligned} \\]

The temperature approximation is:

\\[ T^h(\mathbf x,t)= \sum\_{a=1}^{n_e}N_a(\mathbf x)T_a(t) \\]

Using backward Euler:

\\[ \dot H\approx \frac{H(T^{n+1})-H(T^n)}{\Delta t} \\]

The element residual becomes:

\\[ \begin{aligned} R_a^e={}& \int\_{\Omega_e} N_a\rho \frac{H^{n+1}-H^n}{\Delta t}d\Omega\\\ &+\int\_{\Omega_e} \nabla N_a\cdot k\nabla T^{n+1}d\Omega\\\ &-\int\_{\Omega_e}N_aQ\\,d\Omega\\\ &+\int\_{\Gamma_e}N_aq\_{\mathrm{out}}d\Gamma \end{aligned} \\]

JAX-FEM can use automatic differentiation to construct:

\\[ K\_{ab}^e= \frac{\partial R_a^e} {\partial T_b^{n+1}} \\]

This eliminates much of the manual tangent derivation for nonlinear temperature-dependent properties.

### **7.3 Solution algorithm**

### **JAX-FEM thermal solution procedure**

The main computational operations are element mapping, quadrature integration, residual assembly, Jacobian evaluation, and sparse linear solving.

GPU acceleration of element calculations does not necessarily imply that every global sparse solve is executed entirely on the GPU. Solver backend selection and host–device communication must be evaluated separately.

### **7.4 Proposed LPBF implementation structure**

```
lpbf_jax_fem/
├── mesh.py
├── scan_path.py
├── material_properties.py
├── heat_source.py
├── thermal_problem.py
├── activation.py
├── mechanical_problem.py
├── time_integration.py
├── calibration.py
├── microstructure_export.py
└── postprocessing.py
```

Conceptually, the simulation loop is:

```
# Framework-independent pseudocodetemperature = initial_temperaturefor step in time_steps:    laser_position = scan_path(step.time)    active_material = activation(step.time)    temperature = solve_thermal_fem(        previous_temperature=temperature,        laser_position=laser_position,        material_mask=active_material,        dt=step.dt    )    save_temperature_history(temperature)    displacement, state = solve_mechanics(        temperature=temperature,        material_mask=active_material    )
```

The primary development requirement is connecting the thermal and mechanical constitutive formulations with the actual JAX-FEM API and handling state evolution and material activation.

### **7.5 Published results and research significance**

The foundational JAX-FEM publication demonstrates the usefulness of differentiable FEM for nonlinear mechanics, accelerated numerical simulation, and inverse design.

For LPBF, the important methodological contribution is the ability to differentiate simulation outputs with respect to material properties and design parameters.

Potential research objectives include:

- Calibrating absorptivity and effective conductivity against experimental melt pools.
- Optimizing laser power along the scan path.
- Minimizing peak thermal gradients and residual stress.
- Identifying temperature-dependent constitutive parameters.
- Training physics-based surrogate models using differentiable simulations.

These are proposed LPBF extensions; they should not all be interpreted as demonstrated built-in capabilities of the original JAX-FEM publication.

## **8. JAX-CP: Differentiable crystal plasticity**

### **8.1 Development and research context**

Crystal plasticity extends conventional plasticity by explicitly modeling the deformation mechanisms associated with individual crystallographic slip systems.

JAX is especially well suited to crystal plasticity because constitutive updates must be performed at large numbers of quadrature points or voxels, and those calculations are highly amenable to vectorization.

For this report, **JAX-CP denotes a JAX-based crystal plasticity computational implementation**. The exact independently maintained software package, associated repository, publication, and current version referred to by that name could not be verified. It is therefore important to distinguish the established crystal plasticity equations from proposed JAX implementations.

[Low-cycle fatigue life prediction of a polycrystalline nickel-base superalloy using crystal plasticity modelling approach](https://images.openai.com/static-rsc-4/RbRRd-jzRBYTQXMknzpDqtzvMVDnWWDcV7Y_tVsyC5pIAEksfudmAGrxTLOnZvBlfuXrJjhkzuYyhOhzekNVE3rVK0EcZpcvXbCFPExv2wzwAtM8zBHHZKAFO4efvrbg2JeY97vHzZj8WRMtA-JttXmFrR1wFfNij-LAY2FsKhY?purpose=inline)

[A Critical Review of von Mises Criterion for Compatible Deformation of Polycrystalline Materials | MDPI](https://images.openai.com/static-rsc-4/i9BGp8Daqk-xUQxC92_dvBl5cyW-x7JcoSIQBqdEYWCBHtv5D6WXuuV5Ij-iuAWUMRN7WNrfLvtom1z7A-0qE2z8gYahkDPtkZQEjlskIkpbFBhAX1rujrguhf_1C21beQof3XC15iDfShd-L9RzsIjVUInkSgKPi6gN6xf4zGA?purpose=inline)

[Gallery — AMITEX general  documentation](https://images.openai.com/static-rsc-4/iabv44bttOQ2XB7upsFqzjnPvqgDxsAZDrYAH-hzEYTw6cvZMWjPFZ7Kl1iTKeYn041ThzeFBleGL24lPjOWHx3_5kSLcOIUSdx5YvuUCNiSvhG9hAl5o94IIchnet3O1RouehlyjaKn986Z0bc_PZAH247Iakl9lJwe7-v8FcA?purpose=inline)

12

### **8.2 Finite-strain crystal plasticity formulation**

The deformation gradient is multiplicatively decomposed:

\\[ \mathbf F=\mathbf F^e\mathbf F^p \\]

For nonisothermal applications, a thermal contribution may be included:

\\[ \mathbf F=\mathbf F^e\mathbf F^p\mathbf F^{th} \\]

with the decomposition order chosen consistently with the constitutive formulation.

The plastic velocity gradient is:

\\[ \mathbf L^p= \dot{\mathbf F}^p(\mathbf F^p)^{-1} \\]

Plastic deformation occurs by crystallographic slip:

\\[ \mathbf L^p= \sum\_{\alpha=1}^{N_s} \dot\gamma^\alpha \mathbf s^\alpha\otimes\mathbf m^\alpha \\]

where \\(\mathbf s^\alpha\\) and \\(\mathbf m^\alpha\\) are the slip direction and slip-plane normal.

For FCC metals such as austenitic 316L stainless steel, twelve primary slip systems are commonly employed:

\\[ \\{111\\}\langle110\rangle \\]

### **8.3 Resolved shear stress**

A finite-strain formulation can define the resolved shear stress through the Mandel stress:

\\[ \tau^\alpha= \mathbf M: (\mathbf s^\alpha\otimes\mathbf m^\alpha) \\]

where:

\\[ \mathbf M=\mathbf F^{eT}\boldsymbol\tau \mathbf F^{e-T} \\]

and \\(\boldsymbol\tau\\) is Kirchhoff stress.

A widely used viscoplastic slip law is:

\\[ \dot\gamma^\alpha= \dot\gamma_0 \left| \frac{\tau^\alpha}{g^\alpha} \right|^{1/m} \operatorname{sign}(\tau^\alpha) \\]

Here \\(g^\alpha\\) is the slip resistance and \\(m\\) is the rate-sensitivity parameter.

### **8.4 Hardening evolution**

The slip resistance evolves according to:

\\[ \dot g^\alpha= \sum\_\beta h\_{\alpha\beta} |\dot\gamma^\beta| \\]

A Voce-type hardening law is:

\\[ h\_{\alpha\beta} = q\_{\alpha\beta}h_0 \left( 1-\frac{g^\beta}{g_s} \right)^a \\]

within its applicable range.

The hardening matrix \\(q\_{\alpha\beta}\\) accounts for self and latent hardening.

### **8.5 Plastic deformation gradient update**

An exponential integration scheme is:

\\[ \mathbf F^p\_{n+1} = \exp(\Delta t\mathbf L^p\_{n+1}) \mathbf F^p_n \\]

The elastic deformation gradient follows from:

\\[ \mathbf F^e\_{n+1} = \mathbf F\_{n+1} (\mathbf F^p\_{n+1})^{-1} \\]

The local constitutive equations are nonlinear and generally require an implicit iterative solution.

### **8.6 Automatic differentiation of crystal plasticity**

JAX can compute the constitutive tangent:

\\[ \mathbb A= \frac{d\mathbf P}{d\mathbf F} \\]

where \\(\mathbf P\\) is the first Piola–Kirchhoff stress.

For a local constitutive problem:

\\[ \mathbf r(\mathbf z,\mathbf F)=0 \\]

with \\(\mathbf z\\) representing slip increments or other internal variables, the consistent tangent is:

\\[ \frac{d\mathbf P}{d\mathbf F} = \frac{\partial\mathbf P}{\partial\mathbf F} - \frac{\partial\mathbf P}{\partial\mathbf z} \left( \frac{\partial\mathbf r}{\partial\mathbf z} \right)^{-1} \frac{\partial\mathbf r}{\partial\mathbf F} \\]

This formulation is important for stable Newton convergence in either CP-FEM or Newton-based CP-FFT.

### **8.7 Illustrative constitutive implementation**

```
import jaximport jax.numpy as jnp@jax.jitdef slip_rates(mandel_stress, schmid_tensors,               slip_resistance, gamma0, m):    tau = jnp.einsum(        "ij,aij->a",        mandel_stress,        schmid_tensors    )    return (        gamma0        * jnp.sign(tau)        * (jnp.abs(tau) / slip_resistance)**(1.0 / m)    )@jax.jitdef plastic_velocity_gradient(rates, schmid):    return jnp.einsum("a,aij->ij", rates, schmid)batched_slip_rates = jax.vmap(    slip_rates,    in_axes=(0, 0, 0, None, None))
```

The example demonstrates vectorized slip-rate evaluation, rather than a complete verified finite-strain crystal plasticity solver.

### **8.8 LPBF-specific crystal plasticity inputs**

| Input                   | Source                                     |
| ----------------------- | ------------------------------------------ |
| Grain orientation       | EBSD, CA, or phase field                   |
| Grain morphology        | Process-informed microstructure simulation |
| Initial residual stress | Thermomechanical FEM                       |
| Temperature history     | Thermal FEM                                |
| Elastic constants       | Experimental measurements or literature    |
| Initial slip resistance | Mechanical calibration                     |
| Hardening parameters    | Tensile tests and inverse identification   |
| Phase fractions         | Phase-field model or characterization      |

The most valuable contribution of crystal plasticity is predicting mechanical anisotropy, local strain concentrations, intergranular stresses, and microstructure-dependent deformation that conventional isotropic plasticity cannot resolve.

## **9. MOOSE: Multiphysics Object-Oriented Simulation Environment**

[Framework System Design Description | MOOSE](https://images.openai.com/static-rsc-4/UY2u-TRSqkfka1UcYSScubqHu__V4lBwQtFo8EXrZk23lvNVNR7Sc60Tc9uDx1PqZvpee8tkXVll6yqlZWwzig9f12UzJMYj1-GpjCaQA_b0u5kykQQ18hGfA-InAvY4cXYt-_8PibTtycUgjPLRAAPMAhg_abaS3XEmwKwo1Cw?purpose=inline)

**MOOSE**

Parallel, nonlinear, coupled multiphysics simulation.

[Official website](https://mooseframework.inl.gov/?utm_source=chatgpt.com) · [GitHub](https://github.com/idaholab/moose?utm_source=chatgpt.com)

### **9.1 Development**

MOOSE was developed at Idaho National Laboratory (INL) as a general-purpose, object-oriented multiphysics simulation framework.

Its software architecture is built around C++, finite element discretization through libMesh, and scalable numerical solution infrastructure through PETSc.

Important development publications include Gaston et al. (2009), Gaston et al. (2015), and Permann et al. (2020).

Unlike JAX-FEM, MOOSE offers extensive established multiphysics infrastructure, including heat transfer, solid mechanics, phase field, automatic differentiation, adaptive mesh refinement, and coupling between independently configured simulations.

### **9.2 MOOSE computational architecture**

### **MOOSE multiphysics software architecture**

### **9.3 Kernel-based mathematical formulation**

MOOSE represents differential equations using modular weak-form contributions known as kernels.

For thermal modeling:

\\[ R_i^T= R_i^{\mathrm{storage}} + R_i^{\mathrm{conduction}} + R_i^{\mathrm{source}} + R_i^{\mathrm{boundary}} \\]

For mechanics:

\\[ R_i^u= \int\_\Omega \nabla N_i:\boldsymbol\sigma\\,d\Omega - \int\_\Omega N_i\mathbf b\\,d\Omega - \int\_{\Gamma_t}N_i\bar{\mathbf t}\\,d\Gamma \\]

For a coupled thermal, mechanical, and phase-field problem:

\\[ \mathbf R= \begin{bmatrix} \mathbf R_T\\\ \mathbf R_u\\\ \mathbf R\_\phi\\\ \mathbf R_c \end{bmatrix} =\mathbf 0 \\]

The corresponding Jacobian has block structure:

\\[ \mathbf J= \begin{bmatrix} J\_{TT}&J\_{Tu}&J\_{T\phi}&J\_{Tc}\\\ J\_{uT}&J\_{uu}&J\_{u\phi}&J\_{uc}\\\ J\_{\phi T}&J\_{\phi u}&J\_{\phi\phi}&J\_{\phi c}\\\ J\_{cT}&J\_{cu}&J\_{c\phi}&J\_{cc} \end{bmatrix} \\]

MOOSE supports automatic differentiation of many residual contributions, including off-diagonal multiphysics couplings.

### **9.4 Jacobian-Free Newton–Krylov method**

A major computational method associated with MOOSE is Jacobian-Free Newton–Krylov (JFNK).

The nonlinear Newton equation is:

\\[ \mathbf J(\mathbf x_k)\Delta\mathbf x = -\mathbf R(\mathbf x_k) \\]

A Krylov solver can approximate Jacobian–vector products through:

\\[ \mathbf J\mathbf v\approx \frac{ \mathbf R(\mathbf x+\varepsilon\mathbf v) -\mathbf R(\mathbf x) }{\varepsilon} \\]

This avoids explicitly constructing the full exact Jacobian for the matrix-free operator, although an assembled approximate Jacobian may still be required for preconditioning.

The preconditioner is essential for efficient multiphysics solution.

A general coupled system can employ:

\\[ \mathbf M^{-1}\mathbf J\Delta\mathbf x = -\mathbf M^{-1}\mathbf R \\]

where \\(\mathbf M\\) is an approximation to the coupled Jacobian.

### **9.5 MOOSE modules relevant to LPBF**

| Module or capability      | LPBF application                         |
| ------------------------- | ---------------------------------------- |
| Heat Transfer             | Transient thermal modeling               |
| Solid Mechanics           | Residual stress and deformation          |
| Phase Field               | Grain evolution and phase transformation |
| Crystal Plasticity        | Slip-based constitutive modeling         |
| Navier–Stokes             | Melt-pool flow development               |
| MultiApps                 | Multiscale or coupled simulations        |
| Transfers                 | Field exchange between applications      |
| Mesh adaptivity           | Moving thermal-gradient refinement       |
| Automatic differentiation | Nonlinear coupled Jacobians              |
| PETSc solvers             | Scalable parallel solution               |

### **9.6 Illustrative MOOSE thermal input structure**

```
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 80
  ny = 24
  nz = 20

  xmin = 0
  xmax = 0.002
  ymin = 0
  ymax = 0.0006
  zmin = 0
  zmax = 0.0005
[]

[Variables]
  [T]
    initial_condition = 300
  []
[]

[Kernels]
  [heat_storage]
    type = HeatConductionTimeDerivative
    variable = T
  []

  [heat_conduction]
    type = HeatConduction
    variable = T
  []
[]

[Materials]
  [thermal_properties]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat density'
    prop_values = '25 600 7900'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 1e-5
  end_time = 0.003
[]

[Outputs]
  exodus = true
[]
```

This is an illustrative starting configuration. A complete LPBF application must add the moving heat-flux boundary condition, convection, radiation, latent heat, temperature-dependent properties, and material activation. Object names and required parameters also need validation against the selected MOOSE version.

### **9.7 MOOSE-based crystal plasticity**

MOOSE solid mechanics contains finite-strain crystal plasticity infrastructure, including stress-update and slip-rate constitutive objects.

A MOOSE crystal plasticity model would solve:

\\[ \nabla_0\cdot\mathbf P=0 \\]

with the constitutive update at each quadrature point determined by the slip evolution equations described in Section 8.

The typical nested solution process includes global mechanical equilibrium, local constitutive iteration, slip-resistance evolution, stress updating, and consistent or approximate tangent evaluation.

MOOSE therefore offers an important baseline for validating a newly developed JAX-based crystal plasticity implementation.

### **9.8 Published results**

The MOOSE foundational publications demonstrate a mature architecture for coupled nonlinear multiphysics, scientific software reuse, and large-scale parallel computation.

A particularly relevant related development is the MOOSE-based phase-field framework described by Tonks et al. (2012), which demonstrates the application of finite element multiphysics infrastructure to microstructure evolution.

For an LPBF research study, MOOSE is well suited to establishing a conventional numerical reference solution before developing a differentiable JAX implementation.

## **10. SIERRA: High-performance computational multiphysics**

[Sierra – Sandia National Laboratories](https://images.openai.com/static-rsc-4/yPC8MKoghxuCmdEqiHaZU93yO6IdpJV4fo2u1g9gwukawuwZB_bJTYEbddteS4fPAmjv7f2Owx76tX9jpMs-R4DwHCf444kxKgy-EeIQFrSrafGBU2DJbqmPOb6A-r9R9NM1KM6cvR0m7cSQ7EW3XJo3StkuaZou_8_yEvx6_5o?purpose=inline)

**SIERRA**

Large-scale computational mechanics and multiphysics software developed at Sandia National Laboratories.

[Sandia National Laboratories](https://www.sandia.gov/?utm_source=chatgpt.com)

### **10.1 Development**

SIERRA is a family of high-performance computational simulation applications developed at Sandia National Laboratories.

Its development has emphasized massively parallel computation, nonlinear material behavior, structural mechanics, thermal analysis, and complex engineering multiphysics.

It differs from JAX and MOOSE in its historical focus on large-scale engineering simulation and national-laboratory computational requirements.

### **10.2 Relevant SIERRA applications**

| Application                         | Primary capability               | LPBF relevance                             |
| ----------------------------------- | -------------------------------- | ------------------------------------------ |
| SIERRA/Aria                         | Thermal and coupled multiphysics | Temperature and heat transfer              |
| SIERRA/Adagio                       | Implicit solid mechanics         | Quasistatic residual stress and distortion |
| SIERRA/Presto                       | Explicit solid mechanics         | Dynamic structural response                |
| SIERRA/Fuego                        | Fluid and thermal transport      | Potential melt-pool flow modeling          |
| SIERRA/Thermal–Fluid infrastructure | Coupled transport analysis       | Multiphysics thermal modeling              |

Actual available capabilities and licensing depend on the application, release, and institutional access.

### **10.3 Thermal and mechanical coupling**

A sequential thermomechanical simulation can be expressed as:

\\[ \mathbf M_T(\mathbf T)\dot{\mathbf T} + \mathbf K_T(\mathbf T)\mathbf T = \mathbf f_T(t) \\]

followed by:

\\[ \mathbf R_u(\mathbf u,\mathbf T,\mathbf z)=0 \\]

where \\(\mathbf z\\) denotes material history variables.

For mechanical simulation:

\\[ \mathbf M_u\ddot{\mathbf u} + \mathbf C_u\dot{\mathbf u} + \mathbf f\_{\mathrm{int}} = \mathbf f\_{\mathrm{ext}} \\]

For quasistatic LPBF residual stress analysis, inertia is usually neglected:

\\[ \mathbf f\_{\mathrm{int}} = \mathbf f\_{\mathrm{ext}} \\]

An implicit mechanical solve then takes the form:

\\[ \mathbf K\_{\mathrm{tan}}\Delta\mathbf u = -\mathbf R_u \\]

### **10.4 Component-scale material addition**

A SIERRA-based LPBF simulation can conceptually employ progressive activation of material regions, temperature-dependent plasticity, thermal-to-mechanical field mapping, and sequential thermal cycles.

A possible workflow is:

Conceptual sequential coupling. The actual field-transfer mechanism must be verified for the available SIERRA release.

### **10.5 Research strengths**

SIERRA is particularly relevant when the investigation emphasizes component-scale residual stresses, distortion, large meshes, realistic engineering geometries, distributed-memory HPC, and established computational mechanics verification practices.

The main challenges for the proposed research are access to the software, available documentation, and the development of external differentiable or FFT-based coupling.

SIERRA should therefore be considered an industrial-scale or laboratory-scale reference framework rather than the primary development environment for an open differentiable solver, unless institutional access is already available.

### **10.6 SIERRA versus MOOSE**

| Criterion                     | MOOSE                                      | SIERRA                                |
| ----------------------------- | ------------------------------------------ | ------------------------------------- |
| Primary focus                 | General scientific multiphysics            | Large-scale engineering multiphysics  |
| Core implementation           | C++                                        | Compiled HPC applications             |
| Parallel computation          | MPI/PETSc                                  | Distributed-memory HPC                |
| Microstructure modeling       | Established phase-field module             | Application-dependent                 |
| Crystal plasticity            | Established solid mechanics infrastructure | Material-model availability dependent |
| Differentiable inverse design | Requires additional development            | Requires additional development       |
| Research accessibility        | Public open-source framework               | Application-dependent access          |
| Recommended role              | Reference multiphysics solver              | Component-scale engineering reference |

# **Part III. FEM–FFT multiscale modeling**

## **11. FFT-based computational micromechanics**

Fast Fourier Transform (FFT) methods offer an efficient alternative to finite element methods for computing the mechanical response of heterogeneous periodic materials.

The foundational spectral homogenization method was introduced by Moulinec and Suquet in the 1990s and subsequently extended to nonlinear materials, viscoplastic polycrystals, finite-strain mechanics, and Galerkin formulations.

FFT methods are especially attractive for LPBF research because experimentally reconstructed and computationally generated microstructures are often represented on regular voxel grids.

[TC2 Meeting](https://images.openai.com/static-rsc-4/uxFqithHGP4Tdtuyn74UOvIGBrB4AQ6fnNMfGzLoCKwPKfWfennE5xOZWO1Sjh6SfX6Q0VrQr4-IOO_GKC6m_OOlAZdz_dJyTZBgmuFulcsGhhYPk3hE4fIcmFdb1vM2iLbPxDifqggWNH8smzlLGXRUUjm84gW5QWfp8jsEYNk?purpose=inline)

**CP-FEM**

Grain-resolved finite element discretization.

[The Düsseldorf Advanced Material Simulation Kit:  DAMASK](https://images.openai.com/static-rsc-4/Qbni6Yfc0cPDsSUHUs_FgDkOY4033qDmUdEvRvze0OZofVwTW34w4LmnMsLqPt89qWaNnpsersueIUdZE_F_Y5JffWErvLneg1JDIahM1E_oHC-8p4s2WgTRuCQfC9oSFi4vzWijpLZrXH9PcKyOzjSibwBvGF6BoXFcbyVteLI?purpose=inline)

**CP-FFT**

Regular voxel grid for spectral homogenization.

### **11.1 Governing equations**

For a periodic representative volume element (RVE):

\\[ \nabla\cdot\boldsymbol\sigma(\mathbf x)=0 \\]

The strain field is decomposed into a macroscopic strain and a periodic fluctuation:

\\[ \boldsymbol\varepsilon(\mathbf x) = \bar{\boldsymbol\varepsilon} + \tilde{\boldsymbol\varepsilon}(\mathbf x) \\]

with:

\\[ \langle\tilde{\boldsymbol\varepsilon}\rangle=0 \\]

The macroscopic stress is:

\\[ \bar{\boldsymbol\sigma} = \frac1{|\Omega|} \int\_\Omega \boldsymbol\sigma(\mathbf x)d\Omega \\]

### **11.2 Reference-medium formulation**

The heterogeneous constitutive equation is:

\\[ \boldsymbol\sigma(\mathbf x) = \mathbb C(\mathbf x): \boldsymbol\varepsilon(\mathbf x) \\]

Introducing a homogeneous reference stiffness \\(\mathbb C^0\\), the polarization field is:

\\[ \boldsymbol\tau(\mathbf x) = \left[ \mathbb C(\mathbf x)-\mathbb C^0 \right]: \boldsymbol\varepsilon(\mathbf x) \\]

Therefore:

\\[ \boldsymbol\sigma = \mathbb C^0:\boldsymbol\varepsilon +\boldsymbol\tau \\]

The corresponding Lippmann–Schwinger equation is:

\\[ \boldsymbol\varepsilon(\mathbf x) = \bar{\boldsymbol\varepsilon} - \boldsymbol\Gamma^0\* \boldsymbol\tau(\mathbf x) \\]

where \\(\boldsymbol\Gamma^0\\) is the Green operator.

Transforming to Fourier space:

\\[ \hat{\boldsymbol\varepsilon}(\boldsymbol\xi) = -\hat{\boldsymbol\Gamma}^0(\boldsymbol\xi): \hat{\boldsymbol\tau}(\boldsymbol\xi) \\]

for nonzero Fourier frequencies, while:

\\[ \hat{\boldsymbol\varepsilon}(\mathbf0) = \bar{\boldsymbol\varepsilon} \\]

### **11.3 Green operator for an isotropic reference medium**

For Lamé constants \\(\lambda_0\\) and \\(\mu_0\\):

\\[ \begin{aligned} \hat\Gamma^0\_{ijkl} ={}& \frac{ \delta\_{ik}\xi_j\xi_l+ \delta\_{il}\xi_j\xi_k+ \delta\_{jk}\xi_i\xi_l+ \delta\_{jl}\xi_i\xi_k }{ 4\mu_0|\boldsymbol\xi|^2 }\\\ &- \frac{\lambda_0+\mu_0} {\mu_0(\lambda_0+2\mu_0)} \frac{\xi_i\xi_j\xi_k\xi_l} {|\boldsymbol\xi|^4} \end{aligned} \\]

The zero-frequency mode is treated separately to enforce the macroscopic loading condition.

### **11.4 Classical FFT iterative solution**

The basic fixed-point algorithm proceeds through local stress calculations, polarization evaluation, Fourier transformation, Green-operator application, inverse transformation, and equilibrium checking.

### **FFT micromechanical solution cycle**

The convergence criterion can be based on the equilibrium residual:

\\[ e\_{\mathrm{eq}}= \frac{\\|\nabla\cdot\boldsymbol\sigma\\|\_{L^2}} {\\|\bar{\boldsymbol\sigma}\\|/L} \\]

where \\(L\\) is a characteristic RVE length. Near zero macroscopic stress, an alternative nonzero normalization is needed.

### **11.5 Modern Galerkin FFT**

For nonlinear finite-strain crystal plasticity, a Galerkin FFT formulation is particularly attractive.

The weak equilibrium equation is:

\\[ \int\_{\Omega_0} \delta\mathbf F:\mathbf P(\mathbf F)\\,dV=0 \\]

Using a spectral projection operator \\(\mathcal G\\):

\\[ \mathcal G:\mathbf P(\mathbf F)=0 \\]

Newton linearization gives:

\\[ \mathcal G: \left[ \mathbb A(\mathbf F^k): \Delta\mathbf F \right] = -\mathcal G:\mathbf P(\mathbf F^k) \\]

where:

\\[ \mathbb A=\frac{d\mathbf P}{d\mathbf F} \\]

The linearized system can be solved through a matrix-free Krylov method.

This formulation fits naturally with JAX because the stress and constitutive tangent can be evaluated through automatic differentiation, while compatibility projection is implemented through FFT operations.

### **11.6 Proposed JAX-based FFT implementation**

```
import jaximport jax.numpy as jnp@jax.jitdef spectral_convolution(polarization, green):    polarization_hat = jnp.fft.fftn(        polarization,        axes=(0, 1, 2)    )    strain_hat = -jnp.einsum(        "...ijkl,...kl->...ij",        green,        polarization_hat    )    strain_hat = strain_hat.at[0, 0, 0].set(0)    fluctuation = jnp.fft.ifftn(        strain_hat,        axes=(0, 1, 2)    ).real    return fluctuation
```

A complete solver would include the Green operator or compatibility projection, constitutive update, nonlinear iteration, macroscopic boundary conditions, time integration, and convergence monitoring.

### **11.7 Existing reference software**

The proposed JAX-based FFT framework should be compared against established open computational micromechanics tools, particularly:

| Software          | Main purpose                                          | Research role                |
| ----------------- | ----------------------------------------------------- | ---------------------------- |
| DAMASK            | Crystal plasticity spectral and mesh-based simulation | Primary CP-FFT reference     |
| ExaCA             | Exascale cellular automata grain growth               | AM microstructure generation |
| ExaConstit        | High-performance crystal plasticity FEM               | CP-FEM reference             |
| PRISMS-Plasticity | Crystal plasticity FEM                                | Constitutive verification    |
| PRISMS-PF         | Parallel phase-field modeling                         | Microstructure evolution     |

Useful repositories and documentation:

[DAMASK](https://damask-multiphysics.org/?utm_source=chatgpt.com) · [ExaCA](https://github.com/LLNL/ExaCA?utm_source=chatgpt.com) · [ExaConstit](https://github.com/LLNL/ExaConstit?utm_source=chatgpt.com) · [PRISMS-Plasticity](https://github.com/prisms-center/plasticity?utm_source=chatgpt.com)

### **11.8 Computational scalability**

A regular-grid FFT solver generally requires approximately:

\\[ O(N\log N) \\]

operations per FFT, where \\(N\\) is the total number of voxels.

The overall nonlinear solution cost also depends on constitutive iterations, Krylov convergence, memory bandwidth, and communication overhead.

### **FFT microstructure resolution and memory requirements**

Voxel grid

32³64³128³256³

Scalar-equivalent fields per voxel

50

Numerical precision

FP32

FP64

Total voxels

# **2,097,152**

Field storage

# **0.78 GiB**

32³

0.01 GiB

64³

0.10 GiB

128³

0.78 GiB

256³

6.25 GiB

Calculated field-storage estimate only. FFT workspaces, complex arrays, solver vectors, constitutive Jacobians, and automatic differentiation can substantially increase peak memory consumption.

# **Part IV. Proposed research simulation**

## **12. Reference material: LPBF 316L stainless steel**

For an initial research investigation, 316L stainless steel is a suitable material because of its widespread use in LPBF, availability of experimental literature, FCC crystal structure, and well-established crystal plasticity modeling approaches.

IN718 and Ti-6Al-4V would also be important extensions, but introduce additional challenges associated with precipitation or solid-state phase transformations.

### **12.1 Proposed baseline processing parameters**

| Parameter              | Symbol         | Initial value |
| ---------------------- | -------------- | ------------- |
| Laser power            | \\(P\\)        | 200 W         |
| Scan velocity          | \\(v\\)        | 800 mm/s      |
| Laser beam radius      | \\(r_b\\)      | 40 µm         |
| Laser absorptivity     | \\(\eta\\)     | 0.35          |
| Hatch spacing          | \\(h\\)        | 90 µm         |
| Powder layer thickness | \\(t_l\\)      | 40 µm         |
| Initial temperature    | \\(T_0\\)      | 300 K         |
| Solidus temperature    | \\(T_s\\)      | 1,650 K       |
| Liquidus temperature   | \\(T_l\\)      | 1,720 K       |
| Latent heat            | \\(L_f\\)      | 260 kJ/kg     |
| Density                | \\(\rho\\)     | 7,900 kg/m³   |
| Specific heat          | \\(c_p\\)      | 600 J/kg·K    |
| Thermal conductivity   | \\(k\\)        | 25 W/m·K      |
| Emissivity             | \\(\epsilon\\) | 0.35          |
| Convection coefficient | \\(h_c\\)      | 20 W/m²·K     |

These are illustrative starting parameters, not a validated 316L material database. Temperature-dependent experimental properties should replace the constant values in a quantitative study.

### **12.2 Calculated characteristic parameters**

Using the proposed baseline conditions, the calculated absorbed laser power is:

\\[ P\_{\mathrm{abs}}=\eta P=70\ \mathrm W \\]

The volumetric energy density is:

\\[ E_v= \frac{200}{800(0.09)(0.04)} \\]

\\[ \boxed{E_v=69.44\ \mathrm{J/mm^3}} \\]

The characteristic thermal diffusivity is:

\\[ \alpha\_{\mathrm{th}} = \frac{k}{\rho c_p} = 5.27\times10^{-6}\ \mathrm{m^2/s} \\]

A beam-radius-based Peclet number is:

\\[ Pe=\frac{vr_b}{2\alpha\_{\mathrm{th}}} \approx3.03 \\]

The beam interaction time is:

\\[ t_b=\frac{2r_b}{v}=100\ \mu s \\]

Volumetric energy density

# **69.44**

J/mm³

Peak absorbed heat flux

# **27.85**

GW/m²

Peclet number

# **3.03**

Radius-based definition

Laser interaction time

# **100**

µs

These values were calculated analytically from the assumed inputs; they are not published or simulated melt-pool results.

## **13. Simulation geometry and discretization**

### **Proposed single-track LPBF model**

Illustrative geometry; layer thickness and laser footprint are not drawn to scale.

### **13.1 Initial mesh and solver settings**

| Setting                     | Recommended starting configuration     |
| --------------------------- | -------------------------------------- |
| Computational domain        | 2 × 0.6 × 0.5 mm                       |
| Scan length                 | 1.5 mm                                 |
| Fine mesh size              | 5–10 µm near melt pool                 |
| Coarse mesh size            | 25–50 µm away from melt pool           |
| Element type                | Hexahedral or tetrahedral FEM          |
| Time integration            | Backward Euler                         |
| Heating time step           | 2–5 µs                                 |
| Cooling time step           | Adaptive, progressively increased      |
| Nonlinear solver            | Newton or JFNK                         |
| Initial nonlinear tolerance | \\(10^{-6}\\), appropriately scaled    |
| Mechanical formulation      | Quasistatic thermoplasticity           |
| Microstructure voxel size   | Initially 1–5 µm                       |
| CP-FFT grid                 | \\(64^3\\), progressing to \\(128^3\\) |

The laser displacement within a time increment should be monitored:

\\[ \Delta x_l=v\Delta t \\]

At a time step of 5 µs:

\\[ \Delta x_l=4\ \mu m \\]

This provides a reasonable initial temporal resolution for a refined mesh, subject to numerical convergence testing.

## **14. FEM-to-FFT data transfer**

The interface between the FEM and FFT models is arguably the most important original development opportunity in this research.

### **14.1 Temperature-field mapping**

A FEM temperature field is represented by:

\\[ T^h(\mathbf x,t) = \sum_aN_a(\mathbf x)T_a(t) \\]

For each FFT voxel center \\(\mathbf x_v\\):

\\[ T_v(t)= \sum_aN_a(\mathbf x_v)T_a(t) \\]

This requires locating the corresponding finite element and evaluating its shape functions.

### **14.2 Mechanical-field mapping**

The following state variables may be transferred:

\\[ \mathcal D_v= \\{ T,\dot T,G,R, \mathbf q, g,\boldsymbol\sigma, \mathbf F^p, \boldsymbol\varepsilon^p \\} \\]

where \\(\mathbf q\\) is a crystal orientation quaternion and \\(g\\) is a grain identifier.

A critical issue is that interpolation of stresses from FEM to a voxel grid does not guarantee that the resulting field is in equilibrium:

\\[ \nabla\cdot\boldsymbol\sigma_v \neq0 \\]

Consequently, the FFT solver must perform a mechanically consistent equilibration step.

Similarly, transferring macroscopic plastic strain directly into crystallographic slip variables is not unique. A constitutively consistent initialization or thermal-history-driven crystal plasticity simulation is required.

### **14.3 One-way versus two-way coupling**

**One-way FEM–FFT**

Efficient sequential prediction of process-informed properties.

**Two-way FE²–FFT**

Concurrent homogenization and macroscale constitutive feedback.

In a concurrent multiscale formulation, macroscopic deformation is prescribed to the microscale RVE:

\\[ \bar{\mathbf F} = \langle\mathbf F\rangle \\]

The homogenized stress is returned to the macroscale model:

\\[ \bar{\mathbf P} = \langle\mathbf P\rangle \\]

The homogenized tangent is:

\\[ \bar{\mathbb A} = \frac{\partial\bar{\mathbf P}} {\partial\bar{\mathbf F}} \\]

Energetic consistency is governed by the Hill–Mandel condition:

\\[ \bar{\mathbf P}:\dot{\bar{\mathbf F}} = \langle\mathbf P:\dot{\mathbf F}\rangle \\]

A one-way FEM–FFT formulation is recommended for the initial study, with two-way FE²–FFT developed subsequently.

## **15. Numerical verification and validation**

An important aspect of the research is demonstrating that differences between JAX-FEM, MOOSE, SIERRA, and FFT solvers arise from physical and numerical assumptions rather than implementation errors.

### **15.1 Verification hierarchy**

| Test                       | Reference                  | Quantity of interest           |
| -------------------------- | -------------------------- | ------------------------------ |
| One-dimensional conduction | Analytical solution        | Temperature error              |
| Moving heat source         | Rosenthal approximation    | Temperature distribution       |
| Phase-change conduction    | Stefan problem             | Interface position             |
| Linear elasticity          | Analytical benchmark       | Stress and displacement        |
| Single-crystal loading     | Constitutive reference     | Slip activation                |
| Two-phase composite        | Homogenization reference   | Effective stiffness            |
| Polycrystal                | DAMASK or CP-FEM           | Stress–strain response         |
| Single LPBF track          | Experimental cross-section | Melt-pool width and depth      |
| Multilayer build           | AM-Bench experiments       | Thermal history and distortion |

### **15.2 Experimental validation quantities**

The thermal solver should be validated against melt-pool width, depth, length, temperature measurements, and cooling rates where available.

The microstructure model should be validated against grain size, aspect ratio, orientation distribution functions, pole figures, and EBSD maps.

The mechanical solver should be validated against tensile stress–strain curves, yield stress, elastic anisotropy, residual stress distributions, and local strain measurements.

Useful error metrics include:

\\[ e_W= \frac{|W\_{\mathrm{sim}}-W\_{\mathrm{exp}}|} {W\_{\mathrm{exp}}}\times100\\% \\]

\\[ \mathrm{RMSE}= \sqrt{ \frac1N \sum\_{i=1}^{N} (y_i^{\mathrm{sim}}-y_i^{\mathrm{exp}})^2 } \\]

### **15.3 Numerical convergence study**

### **Recommended convergence test matrix**

| Resolution     | Coarse | Medium | Fine   |
| -------------- | ------ | ------ | ------ |
| FEM mesh       | 20 µm  | 10 µm  | 5 µm   |
| Time step      | 10 µs  | 5 µs   | 2.5 µs |
| FFT grid       | 32³    | 64³    | 128³   |
| RVE replicates | 3      | 5      | 10     |

Proposed test conditions. Spatial resolution, time step, RVE physical dimensions, and statistical sampling should be varied independently where possible.

### **15.4 Experimental benchmark source**

The NIST Additive Manufacturing Benchmark Test Series (AM-Bench) is particularly valuable for rigorous validation.

Potential datasets include single-track melting, thermal measurements, residual stresses, part distortion, and material microstructure characterization.

[NIST AM-Bench official website](https://www.nist.gov/ambench?utm_source=chatgpt.com)

A study using AM-Bench can provide greater reproducibility and scientific value than one validated solely against unspecified literature images.

# **Part V. Published literature and state of the art**

## **16. Important published research**

The following literature establishes the scientific basis for the proposed framework.

Live literature retrieval is unavailable in this session. The publications below are established references that I can identify from existing knowledge, rather than the results of a new systematic database search. Consequently, the most recent JAX-CP developments, 2025–2026 publications, software versions, and exact quantitative performance claims require further verification.

### **16.1 LPBF thermal, fluid, and process modeling**

| Publication              | Method                                  | Principal contribution                                                                       |
| ------------------------ | --------------------------------------- | -------------------------------------------------------------------------------------------- |
| King et al. (2015)       | Review of laser–material interactions   | Established the importance of melt-pool physics, process regimes, and computational modeling |
| Khairallah et al. (2016) | High-fidelity powder-scale multiphysics | Demonstrated the importance of recoil pressure and Marangoni flow in defect formation        |
| DebRoy et al. (2018)     | Comprehensive AM review                 | Connected process parameters, transport phenomena, microstructure, and properties            |
| Markl and Körner (2016)  | Multiscale simulation review            | Examined computational methods across powder, melt-pool, and component scales                |
| Rodgers et al. (2017)    | Grain-growth simulation                 | Demonstrated computational modeling of AM microstructure morphology                          |
| Lindgren et al. (2016)   | Thermomechanical simulation review      | Established constitutive modeling requirements for manufacturing process simulation          |

### **16.2 Differentiable numerical methods and multiphysics**

| Publication            | Framework                 | Principal contribution                                                 |
| ---------------------- | ------------------------- | ---------------------------------------------------------------------- |
| Bradbury et al. (2018) | JAX                       | Composable automatic differentiation and numerical transformations     |
| Xue et al. (2023)      | JAX-FEM                   | Differentiable GPU-oriented finite element analysis and inverse design |
| Gaston et al. (2009)   | MOOSE                     | General parallel multiphysics numerical infrastructure                 |
| Gaston et al. (2015)   | MOOSE                     | Multiphysics coupling architecture                                     |
| Permann et al. (2020)  | MOOSE                     | Mature open-source multiphysics computational framework                |
| Tonks et al. (2012)    | MOOSE-related phase field | Coupled phase-field finite element modeling                            |

### **16.3 Crystal plasticity and FFT modeling**

| Publication                | Numerical method              | Principal contribution                                                |
| -------------------------- | ----------------------------- | --------------------------------------------------------------------- |
| Roters et al. (2010)       | CP-FEM                        | Comprehensive finite-strain crystal plasticity constitutive framework |
| Moulinec and Suquet (1998) | FFT homogenization            | Nonlinear heterogeneous material simulation through spectral methods  |
| Lebensohn (2001)           | Viscoplastic FFT              | Spatially resolved deformation of polycrystals                        |
| Eisenlohr et al. (2013)    | Finite-strain spectral method | Crystal plasticity through spectral equilibrium solving               |
| de Geus et al. (2017)      | Galerkin FFT                  | Finite-strain nonlinear FFT formulation                               |
| Roters et al. (2019)       | DAMASK                        | Unified crystal plasticity and multiphysics materials simulation      |

### **16.4 Collective interpretation of published results**

Four important conclusions follow from this body of literature.

First, high-fidelity LPBF research establishes that conduction-only models cannot reproduce every melt-pool phenomenon, particularly keyholing and defect formation. However, appropriately calibrated conduction models remain valuable for computing large-scale thermal histories.

Second, microstructure prediction requires more information than the maximum temperature alone. Thermal gradients, interface velocity, remelting, and crystallographic competition are essential.

Third, crystal plasticity spectral methods offer an efficient means of resolving heterogeneous mechanical fields in regular-grid microstructures.

Finally, differentiable FEM introduces a significant opportunity to integrate parameter identification and inverse design into conventional process modeling.

What remains scientifically interesting is the integration of these capabilities into a verified, process-informed, differentiable multiscale framework.

## **17. Comparative assessment of the five frameworks**

The following assessment reflects their architectural suitability for the proposed study, not measured performance benchmarks.

| Capability                    | JAX                    | JAX-FEM               | JAX-CP approach          | MOOSE                      | SIERRA                |
| ----------------------------- | ---------------------- | --------------------- | ------------------------ | -------------------------- | --------------------- |
| Automatic differentiation     | Native                 | Native through JAX    | Native through JAX       | Supported for many objects | Application-dependent |
| GPU-oriented development      | Strong                 | Strong                | Strong potential         | Configuration-dependent    | Application-dependent |
| Thermal FEM                   | Develop                | Suitable              | Not primary              | Established                | Established           |
| Thermomechanics               | Develop                | Suitable              | Grain-scale focus        | Established                | Established           |
| Phase field                   | Develop                | Develop or extend     | Not primary              | Established                | Verify                |
| Crystal plasticity            | Develop                | Extend                | Primary purpose          | Established                | Verify                |
| FFT solver                    | Native FFT operations  | External coupling     | Suitable for development | External coupling          | External coupling     |
| Distributed-memory HPC        | Supported abstractions | Verify implementation | Develop                  | Established                | Established           |
| Differentiable inverse design | Strong                 | Strong                | Strong potential         | Additional work            | Additional work       |
| LPBF-specific activation      | Develop                | Extend                | Not primary              | Configure or develop       | Verify                |
| Open research development     | Strong                 | Strong                | Package-dependent        | Strong                     | Access-dependent      |

### **Recommended division of computational responsibilities**

**JAX-FEM**

Primary thermal and thermomechanical research solver

**JAX-based CP-FFT**

Primary original micromechanical solver development

**MOOSE**

Reference multiphysics and constitutive verification

**DAMASK**

Independent crystal plasticity FFT validation

**SIERRA**

Optional large-scale engineering validation

**ExaCA / phase field**

Process-informed microstructure generation

# **Part VI. Proposed research contribution**

## **18. Research gap and novelty**

A generic FEM-to-FFT coupling is not sufficient by itself to establish scientific novelty. Sequential thermal modeling, microstructure generation, and crystal plasticity homogenization are already established research directions.

A stronger contribution would be to develop one or more of the following innovations.

| Research direction                   | Potential contribution                                             | Difficulty  |
| ------------------------------------ | ------------------------------------------------------------------ | ----------- |
| Differentiable FEM–FFT               | End-to-end process sensitivity analysis                            | High        |
| GPU-native CP-FFT                    | Accelerated micromechanics in JAX                                  | Medium–high |
| Process-informed initial states      | Physically consistent residual stress and plastic history transfer | High        |
| Inverse thermal calibration          | Identification of effective LPBF material and laser parameters     | Medium      |
| Thermally coupled crystal plasticity | Prediction of grain-scale stress during manufacturing              | High        |
| Concurrent FE²–FFT                   | Microstructure-sensitive component-scale response                  | Very high   |
| Multifidelity process modeling       | Coupling calibrated conduction models with selective CFD           | High        |

The most promising direction is a differentiable FEM–FFT framework that connects process parameters to microstructure-sensitive mechanical behavior.

## **19. Proposed research title**

**Development of a Differentiable Multiscale FEM–FFT Framework for Process–Structure–Property Modeling of Laser Powder Bed Fusion**

An alternative emphasizing computational mechanics:

**GPU-Accelerated Differentiable Thermal Finite Element and Crystal Plasticity Spectral Modeling of Additively Manufactured Metals**

### **19.1 Research hypothesis**

A process-informed, differentiable FEM–FFT framework can predict microstructure-dependent mechanical responses in LPBF materials while enabling efficient calibration and sensitivity analysis of process and constitutive parameters.

This is a hypothesis to be tested rather than an established result.

### **19.2 Research objectives**

The proposed objectives are to develop a JAX-based transient LPBF thermal solver, verify its numerical predictions against MOOSE, generate three-dimensional microstructures using process-derived thermal histories, implement a JAX-based crystal plasticity FFT solver, validate its constitutive behavior against DAMASK or CP-FEM, and demonstrate the influence of LPBF process parameters on mechanical response.

A subsequent objective is to integrate automatic differentiation across the computational workflow for inverse process design.

## **20. Differentiable process optimization**

The process parameter vector can be written as:

\\[ \boldsymbol\theta= [P,v,h,t_l,\eta,T\_{\mathrm{preheat}}] \\]

The simulation chain is:

\\[ \mathbf T=\mathcal M_T(\boldsymbol\theta) \\]

\\[ \mathcal S=\mathcal M_S(\mathbf T) \\]

\\[ \boldsymbol\sigma= \mathcal M\_{\mathrm{CP}}(\mathcal S,\mathbf T) \\]

where \\(\mathcal S\\) represents the microstructure.

A possible multiobjective optimization is:

\\[ \begin{aligned} \min\_{\boldsymbol\theta}\ J={}& w_1J\_{\mathrm{meltpool}} +w_2J\_{\mathrm{stress}}\\\ &+w_3J\_{\mathrm{anisotropy}} +w_4J\_{\mathrm{cost}} \end{aligned} \\]

with constraints on acceptable melt-pool geometry, processing windows, material properties, and productivity.

For example, the residual stress objective could be:

\\[ J\_{\mathrm{stress}}= \frac1{|\Omega|} \int\_\Omega \sigma\_{\mathrm{vm}}^2dV \\]

Mechanical anisotropy can be quantified by:

\\[ A_y= \frac{ |\sigma_y^{BD}-\sigma_y^{TD}| }{ (\sigma_y^{BD}+\sigma_y^{TD})/2 } \\]

where BD and TD denote build and transverse directions.

### **20.1 The differentiability challenge**

A key scientific challenge is that conventional microstructure simulation involves discrete events such as grain nucleation, grain competition, remelting, and changes in grain labels.

Although JAX can differentiate smooth numerical operations, a conventional cellular automata model is not automatically differentiable in a physically meaningful way.

Possible research approaches include continuous phase-field models, smooth microstructure descriptors, differentiable surrogate models, stochastic sensitivity estimation, or hybrid gradient-based and derivative-free optimization.

This challenge could itself form a substantial original contribution.

## **21. Suggested computational experiments**

### **21.1 Process parameter study**

A first factorial study could use three laser powers and three scanning velocities.

### **Proposed nine-case LPBF parameter matrix**

Volumetric energy density, J/mm³; hatch = 90 µm, layer = 40 µm.Calculated energy densities, not predictions of process quality or porosity.

For each case, the research should extract melt-pool width and depth, peak temperature, cooling rate, thermal gradient, solidification velocity, grain morphology, texture, residual stress, yield stress, and mechanical anisotropy.

### **21.2 Solver benchmarking study**

| Experiment                | JAX-FEM | MOOSE                   | SIERRA   | JAX-FFT  | DAMASK    |
| ------------------------- | ------- | ----------------------- | -------- | -------- | --------- |
| Single-track thermal      | Yes     | Yes                     | Optional | —        | —         |
| Multilayer thermal        | Yes     | Yes                     | Optional | —        | —         |
| Residual stress           | Yes     | Yes                     | Optional | Optional | Optional  |
| Single-crystal plasticity | Yes     | Yes                     | Verify   | Yes      | Yes       |
| Polycrystal tension       | Yes     | Yes                     | Verify   | Yes      | Yes       |
| Parameter identification  | Yes     | Reference               | —        | Yes      | Reference |
| GPU scalability           | Yes     | Configuration-dependent | Verify   | Yes      | Reference |

Benchmarking should record wall-clock time, peak memory, mesh or voxel resolution, solver iterations, numerical error, and hardware specifications.

JAX compilation time should be reported separately from subsequent execution time so that performance comparisons are meaningful.

## **22. Research development roadmap**

### **Illustrative 18-month research program**

Proposed project schedule, assuming access to computational resources and appropriate experimental reference data.

### **22.1 Potential publication strategy**

| Paper   | Proposed subject                         | Main contribution                                                    |
| ------- | ---------------------------------------- | -------------------------------------------------------------------- |
| Paper 1 | Differentiable thermal FEM for LPBF      | Numerical implementation, verification, thermal calibration          |
| Paper 2 | JAX-based CP-FFT solver                  | Constitutive differentiation, spectral solution, GPU benchmarking    |
| Paper 3 | Process-informed FEM–FFT modeling        | Thermal history, microstructure, and mechanical property linkage     |
| Paper 4 | Differentiable LPBF process optimization | Process sensitivity, inverse design, microstructure-aware objectives |

Potential journals include *Additive Manufacturing*, *Computational Materials Science*, *Computer Methods in Applied Mechanics and Engineering*, *International Journal of Plasticity*, and *Computer Physics Communications*.

# **Part VII. References**

## **23. Selected academic bibliography**

### **Additive manufacturing and process physics**

1\. King, W. E., et al. (2015). Laser powder bed fusion additive manufacturing of metals; physics, computational, and materials challenges. *Applied Physics Reviews, 2*, 041304. [DOI](https://doi.org/10.1063/1.4937809?utm_source=chatgpt.com)

2\. Khairallah, S. A., Anderson, A. T., Rubenchik, A., & King, W. E. (2016). Laser powder-bed fusion additive manufacturing: Physics of complex melt flow and formation mechanisms of pores, spatter, and denudation zones. *Acta Materialia, 108*, 36–45. [DOI](https://doi.org/10.1016/j.actamat.2016.02.014?utm_source=chatgpt.com)

3\. DebRoy, T., et al. (2018). Additive manufacturing of metallic components – Process, structure and properties. *Progress in Materials Science, 92*, 112–224. [DOI](https://doi.org/10.1016/j.pmatsci.2017.10.001?utm_source=chatgpt.com)

4\. Markl, M., & Körner, C. (2016). Multiscale modeling of powder bed–based additive manufacturing. *Annual Review of Materials Research, 46*, 93–123. [DOI](https://doi.org/10.1146/annurev-matsci-070115-032158?utm_source=chatgpt.com)

5\. Rodgers, T. M., Madison, J. D., & Tikare, V. (2017). Simulation of metal additive manufacturing microstructures using kinetic Monte Carlo. *Computational Materials Science, 135*, 78–89.

6\. Lindgren, L.-E., et al. (2016). Simulation of additive manufacturing using coupled constitutive and microstructure models. *Additive Manufacturing, 12*, 144–158.

### **JAX and differentiable finite element modeling**

7\. Bradbury, J., et al. (2018). JAX: Composable transformations of Python+NumPy programs. Software. [Repository](https://github.com/jax-ml/jax?utm_source=chatgpt.com)

8\. Xue, T., et al. (2023). JAX-FEM: A differentiable GPU-accelerated 3D finite element solver for automatic inverse design and mechanistic data science. *Computer Physics Communications, 291*, 108802. [DOI](https://doi.org/10.1016/j.cpc.2023.108802?utm_source=chatgpt.com)

### **MOOSE and phase-field modeling**

9\. Gaston, D., Newman, C., Hansen, G., & Lebrun-Grandié, D. (2009). MOOSE: A parallel computational framework for coupled systems of nonlinear equations. *Nuclear Engineering and Design, 239*, 1768–1778. [DOI](https://doi.org/10.1016/j.nucengdes.2009.05.021?utm_source=chatgpt.com)

10\. Gaston, D. R., et al. (2015). Physics-based multiscale coupling for full core nuclear reactor simulation. *Annals of Nuclear Energy, 84*, 45–54.

11\. Permann, C. J., et al. (2020). MOOSE: Enabling massively parallel multiphysics simulation. *SoftwareX, 11*, 100430. [DOI](https://doi.org/10.1016/j.softx.2020.100430?utm_source=chatgpt.com)

12\. Tonks, M. R., Gaston, D., Millett, P. C., Andrs, D., & Talbot, P. (2012). An object-oriented finite element framework for multiphysics phase field simulations. *Computational Materials Science, 51*, 20–29. [DOI](https://doi.org/10.1016/j.commatsci.2011.07.028?utm_source=chatgpt.com)

13\. Chen, L.-Q. (2002). Phase-field models for microstructure evolution. *Annual Review of Materials Research, 32*, 113–140. [DOI](https://doi.org/10.1146/annurev.matsci.32.112001.132041?utm_source=chatgpt.com)

### **Crystal plasticity and spectral homogenization**

14\. Roters, F., et al. (2010). Overview of constitutive laws, kinematics, homogenization and multiscale methods in crystal plasticity finite-element modeling: Theory, experiments, applications. *Acta Materialia, 58*, 1152–1211. [DOI](https://doi.org/10.1016/j.actamat.2009.10.058?utm_source=chatgpt.com)

15\. Moulinec, H., & Suquet, P. (1998). A numerical method for computing the overall response of nonlinear composites with complex microstructure. *Computer Methods in Applied Mechanics and Engineering, 157*, 69–94. [DOI](https://doi.org/10.1016/S0045-7825\(97\)00218-1?utm_source=chatgpt.com)

16\. Lebensohn, R. A. (2001). N-site modeling of a 3D viscoplastic polycrystal using fast Fourier transform. *Acta Materialia, 49*, 2723–2737. [DOI](https://doi.org/10.1016/S1359-6454\(01\)00172-0?utm_source=chatgpt.com)

17\. Eisenlohr, P., Diehl, M., Lebensohn, R. A., & Roters, F. (2013). A spectral method solution to crystal elasto-viscoplasticity at finite strains. *International Journal of Plasticity, 46*, 37–53. [DOI](https://doi.org/10.1016/j.ijplas.2012.09.012?utm_source=chatgpt.com)

18\. de Geus, T. W. J., Vondřejc, J., Zeman, J., Peerlings, R. H. J., & Geers, M. G. D. (2017). Finite strain FFT-based non-linear solvers made simple. *Computer Methods in Applied Mechanics and Engineering, 318*, 412–430. [DOI](https://doi.org/10.1016/j.cma.2016.12.032?utm_source=chatgpt.com)

19\. Roters, F., et al. (2019). DAMASK – The Düsseldorf Advanced Material Simulation Kit for modeling multi-physics crystal plasticity, thermal, and damage phenomena. *Computational Materials Science, 158*, 420–478. [DOI](https://doi.org/10.1016/j.commatsci.2018.04.030?utm_source=chatgpt.com)

20\. Vondřejc, J., Zeman, J., & Marek, I. (2014). An FFT-based Galerkin method for homogenization of periodic media. *Computers & Mathematics with Applications, 68*, 156–173.

### **23.1 Recommended reading priority**

| Priority | Reference                | Reason                                    |
| -------- | ------------------------ | ----------------------------------------- |
| 1        | Xue et al. (2023)        | Primary differentiable FEM architecture   |
| 2        | DebRoy et al. (2018)     | Comprehensive LPBF process physics        |
| 3        | de Geus et al. (2017)    | Modern nonlinear FFT formulation          |
| 4        | Roters et al. (2019)     | CP-FFT reference implementation           |
| 5        | Roters et al. (2010)     | Crystal plasticity fundamentals           |
| 6        | Permann et al. (2020)    | MOOSE development and methods             |
| 7        | Khairallah et al. (2016) | Melt-pool physics and limitations         |
| 8        | Tonks et al. (2012)      | Multiphysics microstructure modeling      |
| 9        | Eisenlohr et al. (2013)  | Finite-strain spectral crystal plasticity |
| 10       | Markl and Körner (2016)  | Multiscale AM research context            |

## **24. Research limitations and outstanding literature verification**

This report provides the theoretical formulations, computational architecture, proposed simulation settings, and foundational literature necessary to initiate the research. However, several limitations should be explicitly recognized:

- **JAX-CP:** The exact named project and its associated development publication have not been independently identified. The report therefore covers established crystal plasticity theory and a proposed JAX-based implementation.
- **SIERRA:** Specific LPBF application publications, user manuals, material-activation capabilities, and release-dependent solver features require verification.
- **Published numerical results:** Exact reported speedups, mesh sizes, melt-pool dimensions, and quantitative validation errors have not been quoted without access to the original publications.
- **Current literature:** The bibliography is not an exhaustive systematic review through November 2026, and the DOI links have not been live-tested.
- **Simulation setup:** The baseline 316L parameters and proposed code examples are research starting points, not validated executable models.
- **Research novelty:** A formal literature gap analysis is required before claiming that the proposed differentiable FEM–FFT framework is the first of its kind.

These limitations are particularly important if this report is used to develop a doctoral proposal or journal publication.

## **25. Final conclusion and recommendation**

For conducting original LPBF process modeling research, I recommend the following computational framework:

**JAX-FEM should be the principal development environment, MOOSE the reference multiphysics framework, and JAX-based CP-FFT the principal original solver contribution.** DAMASK should be used for crystal plasticity validation, while SIERRA can provide additional engineering-scale verification if access is available.

The strongest scientific contribution would not simply be reproducing conventional LPBF thermal modeling. It would be establishing a computationally efficient and differentiable relationship between manufacturing conditions, thermal history, grain structure, crystallographic deformation, and final mechanical properties.

This direction combines computational mechanics, materials science, additive manufacturing, high-performance computing, and inverse design within a coherent research program capable of supporting several substantive journal publications.