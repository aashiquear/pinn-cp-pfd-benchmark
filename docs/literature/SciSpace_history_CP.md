# Mechanics-Informed Neural Networks for Plasticity, Elastoplasticity, and Viscoplasticity: A Critical Literature Review (2020–2026)

Research focus: Constitutive modeling · Thermodynamics-informed learning · History-dependent material behavior · Crystal plasticity · FEM acceleration

## 1. Executive summary

Mechanics-informed neural networks have developed substantially from 2020 to October 2026. The most important advances are not simply improvements in stress–strain prediction accuracy, but new ways to represent internal material states, plastic flow, dissipation, and loading-history dependence while maintaining compatibility with finite element methods (FEM).

The literature reveals five major directions:

1. Thermodynamics-informed constitutive networks: Learn free-energy and dissipation potentials while enforcing mechanical consistency.
2. Physics-informed stress-integration networks: Replace or accelerate traditional elastoplastic return-mapping algorithms.
3. Recurrent and state-space neural networks: Learn loading-history dependence through hidden internal states or physically interpretable plastic variables.
4. Continuous-time neural evolution models: Use neural ordinary or controlled differential equations to reduce dependence on the numerical time discretization.
5. Physics-embedded and multiscale surrogates: Retain selected analytical constitutive models inside neural architectures or replace expensive crystal-plasticity homogenization calculations.

### The most important finding for crystal plasticity

History dependence is not fundamentally the problem that prevents conventional physics-informed neural networks (PINNs) from solving crystal plasticity. The difficulty arises from how the history is represented, how irreversible constitutive evolution is integrated, and how those updates interact with the nonlinear FEM equilibrium solver.

A conventional PINN approximating stress as

$$ \boldsymbol{\sigma}\_\theta = f\_\theta(\boldsymbol{\varepsilon},t,\mathbf{x}) $$

does not generally have enough information to distinguish two material points with identical current strain but different previous plastic loading.

A more appropriate representation is

$$ \boxed{ \begin{aligned} \dot{\mathbf z}&= g\_\theta(\mathbf F,\dot{\mathbf F},\mathbf z,T,\boldsymbol\mu),\\\ \boldsymbol\sigma&= f\_\theta(\mathbf F,\mathbf z,T,\boldsymbol\mu), \end{aligned} } $$

where $\mathbf z$ contains the internal variables encoding material memory and $\boldsymbol\mu$ contains material and microstructural parameters.

This is the central idea behind continuous-time thermodynamics-based networks, recurrent neural operators, and several recent crystal-plasticity surrogates.

Three especially relevant publications are:

- [Evolution TANN — Masi and Stefanou (2023)](https://scispace.com/papers/evolution-tann-and-the-identification-of-internal-variables-2uuqc013?utm_source=chatgpt.com), which learns continuous-time evolution equations for internal variables.
- [Learning macroscopic internal variables and history dependence — Liu et al. (2023)](https://scispace.com/papers/learning-macroscopic-internal-variables-and-history-3bss1ord?utm_source=chatgpt.com), which develops a recurrent neural operator for crystal-plasticity-based homogenization.
- [Incremental Neural Controlled Differential Equations — He and Semnani (2024)](https://scispace.com/papers/incremental-neural-controlled-differential-equations-for-1z3p77kxw4?utm_source=chatgpt.com), which directly addresses stability, convergence, and time-discretization sensitivity of history-dependent neural constitutive models embedded in FEM.

However, eliminating the local plastic return-mapping iteration does not automatically eliminate the global FEM Newton iterations. The latter enforce structural equilibrium and may still be required.

### Scope and evidence limitations

This review covers foundational work from 2020–2022 and developments through October 2026, emphasizing elastoplasticity, viscoplasticity, finite-strain constitutive laws, and crystal plasticity.

The review draws on searches of the SciSpace scholarly index and available research-paper excerpts. It is a broad, targeted literature synthesis, not a certified exhaustive systematic review. Some 2025–2026 studies are preprints or repository publications rather than independently replicated results.

The equations below provide a substantial common mathematical framework and important model-specific formulations. They are not a verbatim reproduction of every equation or parameter table in every publication; verifying that level of completeness would require full-text and supplementary-material extraction for each paper.

Likewise, I distinguish numerical results explicitly reported in the available sources from qualitative claims, unreported hardware specifications, and analytical comparisons developed in this review.

## 2. Literature landscape: major research developments

The most useful way to classify this literature is by the part of classical computational mechanics that the neural network replaces.

Where neural learning enters the mechanics workflow

Loading history and boundary conditions

$\mathbf F(t),\boldsymbol\varepsilon(t),T(t),\boldsymbol\mu$

Classical constitutive solver

Elastic trial state

Plastic yield check

Slip / return mapping and hardening updates

Stress and consistent tangent

Mechanics-informed surrogate

Learned state evolution

Energy / yield / flow constraints

Direct update or neural ODE integration

Stress and neural tangent

Global FEM equilibrium or neural field solver

May still require Newton iterations, equilibrium checks and history-state commitment

Conceptual comparison. Hybrid models may retain some or all of the classical constitutive integration within the neural architecture.

### 2.1 Foundational and established studies, 2020–2024

| Study                                                                                                                                                                 | Main innovation                                                                         | Material models and benchmarks                                                                      |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| [Masi et al. (2021), TANN](https://scispace.com/papers/thermodynamics-based-artificial-neural-networks-for-3mnm9o4bfv?utm_source=chatgpt.com)                         | Free energy, dissipation and thermodynamic relationships encoded in neural architecture | Hyperplasticity, hypoplasticity, strain hardening and softening; comparison with conventional ANN   |
| [Benabou (2021)](https://scispace.com/papers/development-of-lstm-networks-for-predicting-viscoplasticity-2ja77c95mg?utm_source=chatgpt.com)                           | LSTM memory for rate- and temperature-history effects                                   | Anand-law viscoplastic solder, strain-rate jumps, temperature changes, loading/unloading            |
| [Borkowski et al. (2022)](https://scispace.com/papers/recurrent-neural-network-based-multiaxial-plasticity-model-17nevmthbh?utm_source=chatgpt.com)                   | Multiaxial RNN with nonnegative plastic-power regularization                            | Radial-return-generated plasticity histories                                                        |
| [Masi and Stefanou (2022), multiscale TANN](https://scispace.com/papers/multiscale-modeling-of-inelastic-materials-with-22m6sipk?utm_source=chatgpt.com)              | Learned thermodynamically meaningful microstructural state variables                    | 2D/3D inelastic lattices, cyclic histories, FEM×TANN                                                |
| [Bonatti, Berisha and Mohr (2022)](https://scispace.com/papers/from-cp-fft-to-cp-rnn-recurrent-neural-network-surrogate-6syx05fw?utm_source=chatgpt.com)              | Discretization-robust RNN surrogate for crystal plasticity                              | 2D/3D CP-FFT aluminum-alloy simulations and explicit FE applications                                |
| [Masi and Stefanou (2023), eTANN](https://scispace.com/papers/evolution-tann-and-the-identification-of-internal-variables-2uuqc013?utm_source=chatgpt.com)            | Continuous-time evolution equations for internal states                                 | Plasticity, viscosity, damage, coupled inelastic behaviors                                          |
| [Liu et al. (2023), RNO](https://scispace.com/papers/learning-macroscopic-internal-variables-and-history-3bss1ord?utm_source=chatgpt.com)                             | Learns low-dimensional macroscopic memory and a continuous-time response operator       | Laminates and 2D/3D crystal-plasticity polycrystals                                                 |
| [Eghtesad, Fuhg and Bouklas (2023), NN-EVP](https://scispace.com/papers/nn-evp-a-physics-informed-neural-network-based-elasto-1lnc8ejxad?utm_source=chatgpt.com)      | Convex neural potentials for grain-size-dependent elasto-viscoplasticity                | Synthetic power-law hardening, experimental metal flow curves, Hall–Petch behavior                  |
| [He and Semnani (2024), INCDE](https://scispace.com/papers/incremental-neural-controlled-differential-equations-for-1z3p77kxw4?utm_source=chatgpt.com)                | Continuous/incremental neural differential equations designed for FE stability          | $J_2$ and Drucker–Prager plasticity; cyclic, monotonic and boundary-value tests                 |
| [Hildebrand and Klinge (2024)](https://scispace.com/papers/hybrid-data-driven-and-physics-informed-regularized-learning-2altaqiqfl?utm_source=chatgpt.com)            | Physics-regularized state-update network with explicit backstress memory                | 3D Armstrong–Frederick cyclic plasticity                                                            |
| [Jadoon, Meyer and Fuhg (2024)](https://scispace.com/papers/automated-model-discovery-of-finite-strain-elastoplasticity-jdzglxgnw28q?utm_source=chatgpt.com)          | Finite-strain thermodynamic potentials for learned isotropic/kinematic hardening        | Armstrong–Frederick, Ohno–Wang and cyclic mild-steel experiments                                    |
| [Zhang, Yin and Sheil (2024), si-PiNet](https://scispace.com/papers/si-pinet-a-novel-stress-integration-method-for-elastoplastic-7h7dx0gfn4mm?utm_source=chatgpt.com) | Physics-constrained elastoplastic stress integration                                    | von Mises, Mohr–Coulomb and Modified Cam-clay                                                       |
| [Maia et al. (2024), PRNN](https://scispace.com/papers/physically-recurrent-neural-network-for-rate-and-path-2lpixeb7is?utm_source=chatgpt.com)                       | Embedded constitutive material points provide physical memory                           | Finite-strain composite, elastic fibers, elasto-viscoplastic matrix, cyclic and rate-changing paths |
| [Benady, Baranger and Chamoin (2024)](https://scispace.com/papers/unsupervised-learning-of-history-dependent-constitutive-2d3syu3h7y?utm_source=chatgpt.com)          | Identifies history-dependent constitutive potentials from partial measurements          | Elastoplastic and elastoviscoplastic tests, displacement/strain observations, noisy data            |

### 2.2 Recent developments, 2025–2026

| Study                                                                                                                                                                                        | Main contribution                                                                           | Practical significance                                                                                    |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| [Flaschel et al. (2025), Generalized Standard Material Networks](https://scispace.com/papers/convex-neural-networks-learn-generalized-standard-material-lb0jwojc6usp?utm_source=chatgpt.com) | Convex free-energy and dissipation networks across elastic, plastic and viscoplastic models | Thermodynamic guarantees, but retains implicit integration and Newton–Raphson                             |
| [Gavris and Sun (2025)](https://scispace.com/papers/discovering-neural-elastoplasticity-from-kinematic-mnxk34wi2a83?utm_source=chatgpt.com)                                                  | Differentiable inverse identification using kinematic observations                          | Reduces dependence on experimentally measured stress labels                                               |
| [Heidenreich (2026), extended minimal state cell](https://scispace.com/papers/recurrent-neural-networks-for-data-driven-constitutive-6i56uixcwt7q?utm_source=chatgpt.com)                    | Compact state-aware recurrent models and transfer learning                                  | Reduced training-data requirements for complex 3D inelastic behavior                                      |
| [Pingulkar, Gupta and Singh (2026), No-Prop PINN](https://scispace.com/papers/no-prop-learning-based-physics-informed-neural-networks-for-1xv5akle7q3n?utm_source=chatgpt.com)               | Replaces local nonlinear return-mapping solution with neural prediction                     | Up to two orders of magnitude faster training than conventional gradient-based training in reported tests |
| [Jones and Fuhg (2026)](https://scispace.com/papers/a-hierarchy-of-thermodynamics-learning-frameworks-for-e64fizia6qqx?utm_source=chatgpt.com)                                               | Compares thermodynamic frameworks within a common neural-ODE architecture                   | Tests include a rate-dependent crystal-plasticity polycrystal                                             |
| [Moon et al. (2026)](https://scispace.com/papers/physics-informed-discovery-of-yield-functions-in-plasticity-wto5ros79axy?utm_source=chatgpt.com)                                            | Learns convex anisotropic yield surfaces from displacement/force observations               | von Mises, Hill 1948 and Yld2000-2d benchmark cases                                                       |
| [Wang et al. (2026), Plasolver](https://scispace.com/papers/plasolver-physics-informed-neural-operators-for-qnb49jyjdok9?utm_source=chatgpt.com)                                             | Physics-pretrained neural operator and optional solver warm start                           | Approximately 1% relative prediction errors and around 100× standalone acceleration claimed               |
| [Hooda (2026), CP-PINN](https://scispace.com/papers/physics-informed-neural-network-with-decoupled-pathways-for-tcz6v93ld8l1?utm_source=chatgpt.com)                                         | Multimaterial crystal-plasticity LSTM pathways                                              | FCC Cu, BCC Fe and HCP Ti, with stress and CRSS prediction                                                |
| [Kocbay (2026)](https://scispace.com/papers/learning-1d-stress-resultant-plasticity-from-2d-continuum-96n8c80uhift?utm_source=chatgpt.com)                                                   | Physics-informed beam-model reduction retaining exact return mapping                        | Microsecond-scale constitutive evaluation reported for reduced 1D models                                  |

The literature is not moving toward one universal neural architecture. It is moving toward specialized models that preserve the minimum physics needed for the intended application.

## 3. Governing equations: common mathematical foundation

The following equations establish the notation needed to compare the papers. They represent standard constitutive formulations and synthesis-level neural approximations; they should not be interpreted as a claim that every cited paper uses precisely the same equations.

### 3.1 Small-strain elastoplasticity

The additive decomposition of strain is

$$ \boxed{\boldsymbol\varepsilon =\boldsymbol\varepsilon^e+\boldsymbol\varepsilon^p} $$

with elastic stress response

$$ \boldsymbol\sigma =\mathbb C: (\boldsymbol\varepsilon-\boldsymbol\varepsilon^p). $$

For isotropic elasticity,

$$ \mathbb C:\boldsymbol\varepsilon^e =\lambda\\,\operatorname{tr}(\boldsymbol\varepsilon^e)\mathbf I +2G\boldsymbol\varepsilon^e, $$

where

$$ G=\frac{E}{2(1+\nu)},\qquad \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}. $$

For $J_2$ plasticity with mixed hardening, let

$$ \mathbf s=\operatorname{dev}\boldsymbol\sigma, \qquad \boldsymbol\xi=\mathbf s-\mathbf X, $$

where $\mathbf X$ is the deviatoric backstress. Then

$$ \sigma\_{\mathrm{eq}} =\sqrt{\frac32\boldsymbol\xi:\boldsymbol\xi}, $$

and

$$ \boxed{ f=\sigma\_{\mathrm{eq}}- [\sigma\_{y0}+R(p)]\leq 0. } $$

The associative plastic-flow rule is

$$ \dot{\boldsymbol\varepsilon}^p =\dot\gamma\frac{\partial f}{\partial\boldsymbol\sigma} =\dot\gamma\frac32 \frac{\boldsymbol\xi}{\sigma\_{\mathrm{eq}}}, $$

with accumulated plastic strain rate

$$ \dot p= \sqrt{\frac23 \dot{\boldsymbol\varepsilon}^{p}: \dot{\boldsymbol\varepsilon}^{p}} =\dot\gamma. $$

Rate-independent plasticity additionally requires the Karush–Kuhn–Tucker conditions:

$$ \boxed{ f\leq0,\qquad \dot\gamma\geq0,\qquad \dot\gamma f=0. } $$

These inequalities create a switching system: the material behaves elastically until yielding becomes active.

That switching behavior is one reason unrestricted smooth neural approximations can struggle near elastic–plastic transitions.

### 3.2 Isotropic and kinematic hardening

A common nonlinear isotropic hardening law is

$$ R(p)=Q\_\infty(1-e^{-bp}), $$

where $Q\_\infty$ controls saturation strength and $b$ controls the hardening rate.

For Armstrong–Frederick kinematic hardening,

$$ \boxed{ \dot{\mathbf X} =\frac23 C_k\dot{\boldsymbol\varepsilon}^p -\gamma_k\mathbf X\dot p. } $$

Here $C_k$ determines the initial kinematic hardening stiffness and $\gamma_k$ controls dynamic recovery.

The backstress is essential for representing the Bauschinger effect: yield behavior changes when the loading direction is reversed.

A neural model that stores only total strain and a scalar accumulated plastic strain generally cannot reproduce arbitrary multiaxial kinematic hardening. It needs directional memory, such as backstress tensors or an equivalent learned state.

This is central to the methods of Hildebrand and Klinge, Jadoon et al., and the recurrent state-space literature.

### 3.3 Viscoplastic constitutive equations

For an overstress-type viscoplastic model, a representative Perzyna law is

$$ \boxed{ \dot\gamma =\dot\gamma_0 \left\langle \frac{f}{K} \right\rangle^n, } $$

where $\langle x\rangle=\max(x,0)$, $\dot\gamma_0$ is a reference rate, $K$ is an overstress scale, and $n$ describes rate sensitivity.

The associative flow direction may remain

$$ \dot{\boldsymbol\varepsilon}^{vp} =\dot\gamma\frac{\partial f} {\partial\boldsymbol\sigma}. $$

Unlike ideal rate-independent plasticity, the material may sustain a positive overstress while deforming viscoplastically.

A representative Anand-type viscoplastic model has the structure

$$ \dot p =A\exp\left(-\frac{Q}{R_gT}\right) \left[ \sinh\left(\frac{\xi\sigma\_{\mathrm{eq}}}{s}\right) \right]^{1/m}, $$

$$ \dot s =h_0 \left|1-\frac{s}{s^\*}\right|^a \operatorname{sign}\left(1-\frac{s}{s^\*}\right) \dot p. $$

A commonly used saturation-strength relation is

$$ s^\*= \hat s \left[ \frac{\dot p}{A} \exp\left(\frac{Q}{R_gT}\right) \right]^{n_A}. $$

Here $s$ is a deformation resistance, $Q$ is activation energy, $R_g$ the gas constant, $T$ absolute temperature, $m$ and $n_A$ rate-related exponents, and $h_0,a,\hat s,\xi,A$ additional model parameters. Variants of the Anand model use different conventions.

Benabou's 2021 study is particularly relevant because it learns the stress responses generated by an Anand-law material under changing temperature and strain rate.

For a neural viscoplasticity model, $\Delta t$ is not merely an arbitrary input: physical time affects the amount of creep or flow.

### 3.4 Finite-strain crystal plasticity

At finite strain, crystal plasticity usually begins with the multiplicative decomposition

$$ \boxed{ \mathbf F=\mathbf F^e\mathbf F^p. } $$

The plastic velocity gradient is

$$ \mathbf L^p =\dot{\mathbf F}^p(\mathbf F^p)^{-1} =\sum\_{\alpha=1}^{N_s} \dot\gamma^\alpha \mathbf s^\alpha\otimes\mathbf m^\alpha, $$

where $\mathbf s^\alpha$ and $\mathbf m^\alpha$ are the slip direction and slip-plane normal of slip system $\alpha$.

The resolved shear stress is obtained by projecting an appropriate stress measure, commonly the Mandel stress $\mathbf M$, onto the Schmid tensor:

$$ \tau^\alpha =\mathbf M: (\mathbf s^\alpha\otimes\mathbf m^\alpha). $$

A common rate-sensitive slip law is

$$ \boxed{ \dot\gamma^\alpha =\dot\gamma_0 \left| \frac{\tau^\alpha}{g^\alpha} \right|^{1/m} \operatorname{sign}(\tau^\alpha), } $$

with slip-system strength evolution

$$ \dot g^\alpha =\sum\_{\beta=1}^{N_s} h\_{\alpha\beta} |\dot\gamma^\beta|. $$

The key crystal-plasticity states therefore include $\mathbf F^p$, slip resistance $g^\alpha$, accumulated slip, and sometimes lattice orientation, dislocation densities, and backstress on individual slip systems.

A model for HCP alloys may additionally require twinning variables.

The important distinction is that polycrystal stress depends on the evolution of all relevant microstructural states, not only the current macroscopic deformation gradient.

### 3.5 Parameter catalogue

| Parameter                                    | Physical meaning                                             | Relevant model families                        |
| -------------------------------------------- | ------------------------------------------------------------ | ---------------------------------------------- |
| $E,\nu,\lambda,G$                        | Elastic stiffness                                            | Elastoplasticity, CP, energy networks          |
| $\sigma\_{y0}$                           | Initial yield strength                                       | $J_2$, cyclic plasticity                   |
| $H,Q\_\infty,b$                          | Hardening slope/saturation                                   | Isotropic-hardening models                     |
| $C_k,\gamma_k$                           | Kinematic hardening and dynamic recovery                     | Armstrong–Frederick-type models                |
| $K,n,\dot\gamma_0$                       | Viscoplastic overstress and rate parameters                  | Perzyna and related laws                       |
| $A,Q,R_g,\xi,m,s,h_0,a,\hat s,n_A$       | Anand-type thermal/viscoplastic parameters                   | Temperature- and rate-history models           |
| $g^\alpha,h\_{\alpha\beta},N_s$          | Slip resistance, hardening interaction and slip-system count | Crystal plasticity                             |
| $\mathbf Q$ or orientation descriptors   | Crystal orientation / texture                                | CP and microstructure-informed models          |
| $d,k\_{\mathrm{HP}}$                     | Grain size and Hall–Petch strengthening coefficient          | NN-EVP                                         |
| $\mathbf z\in\mathbb R^k$                | Learned internal-state vector                                | RNO, eTANN, INCDE, RNN                         |
| $\theta,L,w$                             | NN weights, depth and width                                  | All neural models                              |
| $\Delta t,\Delta\boldsymbol\varepsilon$  | Time and deformation increments                              | Incremental networks and numerical integration |
| $w\_\sigma,w\_\psi,w_D,w\_{\mathrm{eq}}$ | Training loss weights                                        | Soft-constrained physics-informed learning     |

These are model variables and parameters, not a compiled set of calibrated numerical values. Different materials require different parameters, and no universal parameter values are justified by these papers.

## 4. How neural networks encode the mechanics

### 4.1 Thermodynamics-informed neural networks

TANN and related models introduce a learned Helmholtz free-energy density,

$$ \psi\_\theta =\psi\_\theta(\boldsymbol\varepsilon,\mathbf z,T). $$

For an isothermal strain-based material with internal variables, the stress is derived from the potential:

$$ \boxed{ \boldsymbol\sigma =\frac{\partial\psi\_\theta} {\partial\boldsymbol\varepsilon}. } $$

Define the thermodynamic forces conjugate to internal variables by

$$ \mathbf A =-\frac{\partial\psi\_\theta}{\partial\mathbf z}. $$

The dissipation inequality is

$$ \boxed{ \mathcal D=\mathbf A\cdot\dot{\mathbf z}\geq0. } $$

A dual dissipation potential $\phi\_\theta^\*(\mathbf A)$ can generate evolution through

$$ \dot{\mathbf z} =\frac{\partial\phi\_\theta^\*}{\partial\mathbf A}. $$

If the potential has the necessary convexity, minimum and normalization properties, thermodynamic dissipation can be enforced by construction.

A representative network structure is

$$ \psi\_\theta(\boldsymbol\varepsilon,\mathbf z), \qquad \phi\_\theta^\*(\mathbf A), \qquad \boldsymbol\sigma=\partial\_\varepsilon\psi\_\theta, \qquad \dot{\mathbf z}=\partial\_{\mathbf A}\phi\_\theta^\*. $$

This is especially useful for viscoplasticity because it connects stress prediction, dissipation and state evolution.

Input-convex neural networks are often used to enforce required convexity. In schematic form, an input-convex layer may be represented as

$$ \mathbf h\_{\ell+1} =a\_\ell( \mathbf W\_\ell^+\mathbf h\_\ell+ \mathbf U\_\ell\mathbf x+\mathbf b\_\ell), $$

with suitable convex, nondecreasing activations and nonnegative weights $\mathbf W\_\ell^+$.

Convexity is not automatically necessary in every input of the free-energy function; it is enforced in the arguments required by the chosen thermodynamic framework.

Main limitation: Thermodynamic consistency does not guarantee accurate identification of unobserved state variables or a computationally inexpensive local update.

For example, the generalized standard material networks of Flaschel et al. explicitly retain an implicit time integration and Newton–Raphson solution process. They improve constitutive expressiveness and admissibility but do not eliminate iterative integration.

### 4.2 Continuous-time evolution TANN

The eTANN approach of Masi and Stefanou addresses a weakness of increment-based neural models: dependence on the chosen numerical step size.

An incremental network learns

$$ \mathbf z\_{n+1} =\mathbf z_n+ \mathcal N\_\theta( \mathbf z_n,\Delta\boldsymbol\varepsilon_n,\Delta t_n). $$

If the learned update has not been constrained appropriately, changing $\Delta t$ or subdividing the same loading path changes the predicted response.

eTANN instead learns

$$ \boxed{ \dot{\mathbf z} =g\_\theta(\boldsymbol\varepsilon,\mathbf z,T) } $$

together with

$$ \boldsymbol\sigma =\partial\_\varepsilon\psi\_\theta. $$

The physical response is then obtained by numerical integration:

$$ \mathbf z\_{n+1} =\mathbf z_n+ \int\_{t_n}^{t\_{n+1}} g\_\theta(\boldsymbol\varepsilon(t),\mathbf z(t),T(t)) \\,dt. $$

The paper also uses an autoencoder to identify suitable latent internal variables from detailed microscopic fields.

This has an important practical consequence: the autoencoder may be needed during training, but the resulting constitutive model can operate using only the learned internal variables, evolution network and free-energy network.

What changes: The neural representation is no longer intrinsically tied to a particular training increment.

What remains: An ODE must still be integrated during inference. Extremely stiff material dynamics may require implicit integration or adaptive substepping.

### 4.3 Recurrent neural operators for crystal plasticity

Liu et al. introduce a particularly relevant history-dependent representation:

$$ \boxed{ \begin{aligned} \boldsymbol\sigma(t) &=f\_\theta(\mathbf F(t),\mathbf z(t)),\\\ \dot z_i(t) &=g\_{\theta,i}(\mathbf F(t),\mathbf z(t)). \end{aligned} } $$

The learned state $\mathbf z(t)$ summarizes past loading.

The paper reports that low-dimensional latent states can effectively represent much more complex microscopic histories: two learned internal variables for its 2D crystal-plasticity setting and five for its 3D setting.

This result is important because a conventional crystal-plasticity RVE may possess hundreds or thousands of underlying slip-related state variables.

It suggests that the macroscopic observable response can sometimes be represented on a much lower-dimensional state manifold.

However, those learned variables need not have a unique interpretation. Multiple internal-state parameterizations may produce the same stress response, and a dimension sufficient for one material or loading family may not generalize to another.

### 4.4 Incremental neural controlled differential equations

The work of He and Semnani directly targets a problem that ordinary recurrent networks encounter when used inside FEM: small changes in a strain increment can produce disproportionately large changes in the predicted stress.

A conceptual neural controlled differential equation has the form

$$ d\mathbf z =\mathbf V\_\theta(\mathbf z,\boldsymbol\varepsilon) \\,d\boldsymbol\varepsilon. $$

For rate-dependent problems, physical time evolution can be added explicitly:

$$ d\mathbf z =\mathbf V\_\theta(\mathbf z,\boldsymbol\varepsilon) \\,d\boldsymbol\varepsilon +\mathbf r\_\theta(\mathbf z,\boldsymbol\varepsilon,T)\\,dt. $$

These are generic forms illustrating the idea, rather than a verbatim transcription of the paper's INCDE algorithm.

For a sufficiently regular learned vector field, a small incremental change should induce a controlled change in the state:

$$ \\|\delta\mathbf z\_{n+1}\\| \lesssim C\\|\delta\mathbf z_n\\|+ L\\|\delta\boldsymbol\varepsilon_n\\|. $$

This matters because a nonlinear FEM solver repeatedly evaluates its constitutive law at trial strain states.

A neural model that produces an inconsistent stress Jacobian or a nonsmooth response may prevent global equilibrium convergence, even if its average supervised stress error is small.

INCDE explicitly investigates these stability and convergence issues in $J_2$ and Drucker–Prager benchmark problems.

### 4.5 Physically recurrent neural networks

The PRNN model of Maia et al. takes a different approach. Instead of asking a generic recurrent network to discover all hidden history dynamics, it embeds conventional constitutive models within a learned encoder–material–decoder architecture.

A schematic representation is

$$ \mathbf F_j^{\mathrm{local}} =\mathcal E\_{\theta,j}(\overline{\mathbf F}), $$

$$ (\boldsymbol\sigma_j,\mathbf q\_{j,n+1}) = \mathcal M_j( \mathbf F_j^{\mathrm{local}}, \mathbf q\_{j,n},\Delta t), $$

$$ \boxed{ \overline{\boldsymbol\sigma} =\mathcal D\_\theta( \boldsymbol\sigma_1,\ldots,\boldsymbol\sigma_N). } $$

Here $\mathcal M_j$ is a physical constitutive model and $\mathbf q_j$ its conventional internal variables.

The network learns a compressed representation of how local material points contribute to the macroscopic response. History is provided by the internal variables already present in those physical models.

This approach often generalizes better than an entirely unconstrained RNN, but it is critical to understand its limitation:

PRNN does not eliminate all plasticity predictor/corrector calculations. The embedded viscoplastic constitutive models still update their internal variables. The computational savings arise primarily from replacing the much larger micromechanical boundary-value problem with a small set of fictitious material points.

## 5. Numerical benchmarks and published performance

Reported computational speedups must be interpreted carefully. A model that accelerates local constitutive integration by 100× has not necessarily accelerated the full FE simulation by 100×.

The following comparison separates accuracy, model size, hardware, and actual reported computational improvements.

### 5.1 Quantitative benchmark comparison

| Paper / model                                                                                                                                                   | Benchmark and quantitative findings                                                                                                                                            | Hardware and timing evidence                                                          |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------- |
| [Bonatti et al., CP-RNN (2022)](https://scispace.com/papers/from-cp-fft-to-cp-rnn-recurrent-neural-network-surrogate-6syx05fw?utm_source=chatgpt.com)           | 2D/3D aluminum-alloy CP surrogate; approximately 5,000–100,000 NN parameters; FE instability mitigation investigated                                                           | No comparable end-to-end wall time verified                                           |
| [Liu et al., RNO (2023)](https://scispace.com/papers/learning-macroscopic-internal-variables-and-history-3bss1ord?utm_source=chatgpt.com)                       | 2 latent variables for 2D CP, 5 for 3D CP; up to 2,400 training realizations; RNO cost assessed against conventional material models                                           | Single-CPU FE evaluation reported; comparable numerical wall times not extracted      |
| [Single-LSTM heterogeneity model (2022)](https://scispace.com/papers/a-single-long-short-term-memory-network-for-enhancing-the-3j33v1ap?utm_source=chatgpt.com) | 2D $J_2$ heterogeneous and anisotropic plasticity; approximately 2% error reported for a tested path-dependent setting                                                     | Hardware/wall time not verified                                                       |
| [Maia et al., PRNN (2024)](https://scispace.com/papers/physically-recurrent-neural-network-for-rate-and-path-2lpixeb7is?utm_source=chatgpt.com)                 | Finite-strain composite micromodel; roughly $10^3$× faster homogenized stress evaluation; tests include cyclic loading, relaxation and rate changes                        | Single Xeon E5-2630 v4 CPU core; cluster node with 128 GB RAM                         |
| [Hildebrand and Klinge (2024)](https://scispace.com/papers/hybrid-data-driven-and-physics-informed-regularized-learning-2altaqiqfl?utm_source=chatgpt.com)      | Armstrong–Frederick cyclic plasticity; physics-regularized networks remained stable in reported tests of up to 50 cycles, whereas an unconstrained model drifted               | No verified runtime                                                                   |
| [Heidenreich, EMSC (2026)](https://scispace.com/papers/recurrent-neural-networks-for-data-driven-constitutive-6i56uixcwt7q?utm_source=chatgpt.com)              | Models with fewer than 25,000 parameters; transfer learning reduces training requirements from thousands to hundreds or, in shared-core settings, fewer than 100 loading paths | Training-data reduction demonstrated; no hardware-normalized wall times verified      |
| [Hooda, CP-PINN (2026)](https://scispace.com/papers/physics-informed-neural-network-with-decoupled-pathways-for-tcz6v93ld8l1?utm_source=chatgpt.com)            | Stress $R^2=0.970$; CRSS $R^2=0.988$; 96.2% uncertainty coverage at a nominal two-sigma level                                                                          | No independently established end-to-end FE speedup                                    |
| [Pingulkar et al., No-Prop (2026)](https://scispace.com/papers/no-prop-learning-based-physics-informed-neural-networks-for-1xv5akle7q3n?utm_source=chatgpt.com) | Up to two orders of magnitude faster training than conventional gradient-based learning                                                                                        | Training speedup, not a demonstrated 100× full-FEM speedup                            |
| [Wang et al., Plasolver (2026)](https://scispace.com/papers/plasolver-physics-informed-neural-operators-for-qnb49jyjdok9?utm_source=chatgpt.com)                | Relative prediction errors on the order of 1%, approximately 100× faster standalone predictions; optional warm start reduces classical solver iterations by approximately 50%  | Reported FNO vs Transolver training time: 6 vs 21 minutes per epoch in one experiment |
| [Kocbay (2026)](https://scispace.com/papers/learning-1d-stress-resultant-plasticity-from-2d-continuum-96n8c80uhift?utm_source=chatgpt.com)                      | 1D beam moment–curvature surrogate with microsecond-scale constitutive evaluations; validated against 2D plane-stress references                                               | Major gains partly arise from dimensional model reduction                             |

NR or “not verified” means the accessible search results and research excerpts did not establish a comparable numerical value. It does not prove the full paper contains no such measurement. $R^2=0.970$ is a goodness-of-fit statistic, not a 3% relative stress-error measurement. Speedup and accuracy figures come from different workloads and should not be directly ranked.

### 5.2 A particularly useful hardware-controlled benchmark

Maia et al.'s physically recurrent finite-strain model provides unusually concrete computing details.

The micromechanical reference has approximately 4,992 integration points and 7,860 degrees of freedom. Its network replaces repeated homogenization analysis with encoded fictitious material-point evaluations.

The study reports execution on a single core of an Intel Xeon E5-2630 v4 processor, on a node with 128 GB RAM running CentOS 7.

Its principal result is an approximately three-orders-of-magnitude reduction in homogenized stress-evaluation runtime relative to the reference micromodel.

That makes the result important for practical engineering: it does not depend on an unspecified large GPU. Nevertheless, it measures the surrogate against a comparatively expensive micromechanical calculation, not against a simple analytically integrated $J_2$ material law.

### 5.3 Evidence synthesis

| No. | Source                                                                                                                                                      | Key insight and qualification                                                                                                                                 |
| --- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | [Masi and Stefanou, eTANN](https://scispace.com/papers/evolution-tann-and-the-identification-of-internal-variables-2uuqc013?utm_source=chatgpt.com)         | Continuous-time state evolution reduces dependence on arbitrary neural training increments, but time integration remains necessary                            |
| 2   | [Liu et al., RNO](https://scispace.com/papers/learning-macroscopic-internal-variables-and-history-3bss1ord?utm_source=chatgpt.com)                          | Low-dimensional learned states can reproduce history-dependent polycrystal response at computational cost approaching conventional constitutive modeling      |
| 3   | [He and Semnani, INCDE](https://scispace.com/papers/incremental-neural-controlled-differential-equations-for-1z3p77kxw4?utm_source=chatgpt.com)             | FEM robustness requires stable stress responses under strain perturbations, not merely accurate standalone stress trajectories                                |
| 4   | [Maia et al., PRNN](https://scispace.com/papers/physically-recurrent-neural-network-for-rate-and-path-2lpixeb7is?utm_source=chatgpt.com)                    | Strong acceleration is possible when neural homogenization replaces a large microscopic FE solve while preserving physical history variables                  |
| 5   | [Flaschel et al., GSM networks](https://scispace.com/papers/convex-neural-networks-learn-generalized-standard-material-lb0jwojc6usp?utm_source=chatgpt.com) | Thermodynamically constrained potentials provide physical consistency, but their learned constitutive response can still require iterative Newton integration |
| 6   | [Pingulkar et al., No-Prop](https://scispace.com/papers/no-prop-learning-based-physics-informed-neural-networks-for-1xv5akle7q3n?utm_source=chatgpt.com)    | Physics-residual learning can replace local return mapping; efficient output-layer fitting can substantially reduce neural training expense                   |
| 7   | [Wang et al., Plasolver](https://scispace.com/papers/plasolver-physics-informed-neural-operators-for-qnb49jyjdok9?utm_source=chatgpt.com)                   | Neural prediction plus a conventional numerical correction provides a route to speed and controllable final solver accuracy                                   |

### 5.4 What accuracy is actually demonstrated?

The strongest claims in this literature usually concern one of four targets:

- Material-point stress-history reconstruction against a reference constitutive law.
- Homogenized stress prediction against an RVE or CP-FFT reference.
- FE force–displacement and field agreement.
- Stability across cyclic loading or changes in increment size.

These are different evaluation tasks.

For example, a network achieving 2% mean stress error on an isolated material-point history may still fail inside implicit FEM if its constitutive tangent is irregular. Conversely, a model with slightly larger pointwise stress error may yield a more robust FE solution because its tangent and internal-state updates are smoother.

For engineering applications, accuracy should be measured at both the material and structural levels.

## 6. Practical computational cost: hardware, time, accuracy and deployment

### 6.1 Hardware requirements

The following is an engineering assessment based on algorithm complexity and reported use cases. The categories are not experimentally standardized hardware specifications.

| Model architecture                        | Training hardware                                | Inference hardware                                | Main practical bottleneck                        |
| ----------------------------------------- | ------------------------------------------------ | ------------------------------------------------- | ------------------------------------------------ |
| Small feedforward PINN / stress-update NN | CPU feasible; GPU useful for extensive training  | CPU practical                                     | Differentiation and tangent computation          |
| TANN / convex potential NN                | CPU or GPU; derivative-intensive                 | CPU possible                                      | Automatic differentiation and implicit evolution |
| LSTM / GRU constitutive network           | GPU advantageous for training long histories     | CPU or GPU                                        | Sequential recurrence and training memory        |
| eTANN / neural ODE                        | GPU commonly beneficial for fitting              | CPU possible for small state spaces               | ODE integration, stiffness and tolerances        |
| INCDE                                     | GPU useful during sequence training              | CPU potentially practical                         | Controlled integration and tangent regularity    |
| PRNN                                      | CPU demonstrated                                 | CPU demonstrated                                  | Embedded constitutive material-point updates     |
| Crystal-plasticity RNO                    | GPU useful for fitting large trajectory datasets | CPU possible                                      | Stress Jacobians, integration and FE coupling    |
| Spatial neural operator                   | GPU advantageous                                 | GPU often advantageous for high-throughput fields | Field size, memory traffic and model evaluation  |

The PRNN single-core result is direct evidence of CPU feasibility. Other classifications describe computational tendencies rather than measured thresholds.

### 6.2 Why a neural material law may be slower than classical plasticity

Consider a standard isotropic $J_2$ plasticity model with linear hardening.

The elastic predictor is

$$ \boldsymbol\sigma^{tr}\_{n+1} =\mathbb C: (\boldsymbol\varepsilon\_{n+1} -\boldsymbol\varepsilon^p_n). $$

If its yield function is nonpositive, the step is elastic.

For the plastic case, isotropic hardening with modulus $H$ gives the classical radial-return relation

$$ \boxed{ \Delta\gamma= \frac{\sigma^{tr}\_{eq}-\sigma\_{y,n}} {3G+H}. } $$

This update can be evaluated using relatively inexpensive tensor algebra. There is little advantage in replacing such a closed-form expression with a large neural network.

For nonlinear hardening, the plastic multiplier instead satisfies an equation such as

$$ r(\Delta\gamma) = \sigma^{tr}\_{eq} -3G\Delta\gamma -\sigma\_{y0} -R(p_n+\Delta\gamma) =0. $$

This may require Newton iteration, but it is only a scalar nonlinear solve.

An extensive neural architecture may still be slower than this numerical procedure, particularly when the network must also produce an accurate algorithmic stiffness.

The greatest opportunity for neural acceleration comes from expensive constitutive updates, many interacting slip systems, repeated microscopic equilibrium solves, or complex multiscale material models.

### 6.3 Total computation time

Let

$$ T\_{\mathrm{offline}} = T\_{\mathrm{data}} +T\_{\mathrm{training}} +T\_{\mathrm{validation}} +T\_{\mathrm{implementation}}. $$

The classical and neural deployment costs can be decomposed approximately as

$$ T\_{\mathrm{classical}} = T\_{\mathrm{global}} +T\_{\mathrm{constitutive}}, $$

$$ T\_{\mathrm{neural}} = T\_{\mathrm{global}}^{NN} +T\_{\mathrm{NN\\,evaluation}} +T\_{\mathrm{NN\\,tangent}} +T\_{\mathrm{state\\,integration}}. $$

The last three terms matter: neural inference time is not necessarily just the forward-pass time.

The total cost over $N$ repeated engineering simulations is approximately

$$ T\_{\mathrm{total}}^{NN} = T\_{\mathrm{offline}}+ N T\_{\mathrm{neural}}. $$

A surrogate becomes economically advantageous when

$$ \boxed{ N> \frac{T\_{\mathrm{offline}}} {T\_{\mathrm{classical}}-T\_{\mathrm{neural}}} } $$

assuming the surrogate is accurate for all the required simulations and the denominator is positive.

This equation explains why surrogate models are often more attractive for optimization, sensitivity analysis and repeated multiscale simulation than for a single material test.

### 6.4 How much of the local speedup survives at the FEM level?

A useful estimate comes from Amdahl's law.

Let $p$ be the fraction of baseline computation time spent in the part being accelerated, and $S\_{\mathrm{local}}$ be its speedup.

Then the best-case total acceleration, ignoring added overhead, is

$$ \boxed{ S\_{\mathrm{total}} = \frac{1}{(1-p)+p/S\_{\mathrm{local}}}. } $$

Interactive estimate: local surrogate speedup vs whole-solver speedup

Baseline runtime spent in constitutive / micromechanical evaluation

90%

Speedup of that local computation

1000×

10×100×1,000×

Local computation

# 1,000×

Estimated full-solver speedup

# 9.91×

Accelerated portion: 90%Unaccelerated: 10%

This is an analytical scenario, not a measured speedup from any individual paper. It assumes unchanged work outside the accelerated region and omits surrogate overhead and convergence changes.

Even a 1,000× local speedup produces only about a 10× full-solver speedup if 90% of baseline runtime was spent in the accelerated computation.

For FEM×RVE crystal plasticity, that may still be an excellent result. For an ordinary elastic–plastic FE model where material integration is cheap, the gains can be much smaller.

### 6.5 Other computational overheads that the papers do not always compare

A fair hardware and time comparison should also account for the cost of generating CP-FEM training trajectories, storing the internal state at every integration point, constructing consistent tangents, transfer between CPUs and GPUs, and repeated global Newton iterations.

Neural training itself can become expensive because plasticity is sequential. Backpropagation through a long loading history may require retaining many intermediate activations.

For a network containing $P$ parameters and a history of $N_t$ time steps, simply counting $P$ does not determine training memory; intermediate states, computational graphs and constitutive differentiation contribute substantially.

An important conclusion from Heidenreich's 2026 work is that transfer learning and compact state-space representations may reduce data requirements. However, fewer training paths do not necessarily mean proportionally less wall-clock time, especially if the generating crystal-plasticity simulations have different costs.

## 7. Why crystal-plasticity PINNs struggle with loading history and FEM predictor steps

This is arguably the most important part of your question.

The difficulties are related, but there are four separate numerical and physical issues that should not be conflated.

### 7.1 Issue 1: Crystal plasticity is not a memoryless stress–strain problem

A conventional neural mapping might use

$$ \boldsymbol\sigma\_{n+1} =f\_\theta(\mathbf F\_{n+1},T\_{n+1}). $$

However, there is generally no unique stress corresponding to a particular $\mathbf F\_{n+1}$.

Consider two loading trajectories:

Same final strain, different loading histories

Path A: monotonic loading to 0.5%

Path B: load to 1%, then unload to 0.5%

Illustrative loading paths, not a reproduced experimental result. Both histories finish at the same strain, but can have different plastic strain, backstress, slip resistance and lattice state.

For crystal plasticity, the stress response depends on stored history:

$$ \boldsymbol\sigma\_{n+1} = \mathcal F( \mathbf F\_{n+1}, \mathbf F^p_n, \\{g_n^\alpha\\}, \\{\gamma_n^\alpha\\}, \mathbf Q_n,\ldots). $$

A neural model must either retain those physical variables or learn a sufficiently informative equivalent state.

How the literature addresses this:

- LSTM and GRU models use recurrent hidden states.
- RNO and eTANN use learned internal-variable evolution equations.
- PRNN retains physical constitutive states inside embedded material models.
- Physics-regularized cyclic plasticity networks propagate explicit backstress and plastic-strain variables.

The important innovation is stateful constitutive learning, rather than a conventional PINN applied independently to every loading step.

### 7.2 Issue 2: The FEM elastic predictor and plastic corrector

A conventional local constitutive update in CP-FEM can be represented schematically as follows.

At the beginning of the increment, the previous plastic state is known:

$$ \mathbf F^p_n,\qquad g^\alpha_n,\qquad \gamma^\alpha_n. $$

Given a candidate deformation gradient, the elastic trial state is

$$ \boxed{ \mathbf F^{e,tr}\_{n+1} = \mathbf F\_{n+1} (\mathbf F^p_n)^{-1}. } $$

The trial stresses and resolved shear stresses are computed from the elastic state.

The plastic corrector then solves nonlinear evolution equations for the slip increments:

$$ \Delta\gamma^\alpha -\Delta t\\, \dot\gamma^\alpha( \boldsymbol\tau\_{n+1}, \mathbf g\_{n+1}) =0. $$

Together with the hardening equations,

$$ \mathbf r( \Delta\boldsymbol\gamma, \mathbf g\_{n+1}, \mathbf F^p\_{n+1})=\mathbf0. $$

A local Newton method may be used:

$$ \mathbf J_r^{(k)} \Delta\mathbf y^{(k)} = -\mathbf r^{(k)}, $$

where $\mathbf y$ collects the unknown slip and hardening updates.

The updated plastic deformation gradient can be obtained with an exponential update of the plastic velocity gradient:

$$ \mathbf F^p\_{n+1} = \exp(\Delta t\\,\mathbf L^p\_{n+1}) \mathbf F^p_n, $$

subject to the assumptions and integration scheme of the chosen CP model.

This process is repeated across many integration points and load steps. The cost increases substantially in multislip, multiphase or multiscale settings.

### 7.3 How neural methods replace or retain the predictor/corrector

| Approach                               | Treatment of local predictor/corrector                                          | Treatment of history                                     |
| -------------------------------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------- |
| Direct stateful stress-update NN       | Can replace the complete local stress update                                    | Learned or prescribed internal states                    |
| No-Prop PINN                           | Explicitly targets replacement of the nonlinear return-mapping solve            | Previous constitutive states remain required             |
| eTANN                                  | Replaces the constitutive evolution law with a learned ODE                      | Integrated internal variables                            |
| RNO / CP-RNN                           | Learns a surrogate for the constitutive history operator                        | Continuous or discretization-aware latent memory         |
| INCDE                                  | Uses incremental neural differential evolution                                  | Evolving latent state designed for incremental stability |
| PRNN                                   | Retains constitutive integration in embedded material cells                     | Physical internal variables                              |
| Generalized standard material networks | Retains implicit integration and Newton solution                                | Thermodynamic internal variables                         |
| Plasolver                              | Learns a spatial elastoplastic solution operator; optional classical correction | Loading-history and internal-variable information        |

Therefore, there is no single answer to whether neural models eliminate the FEM predictor step.

Some remove the local nonlinear return-mapping solve. Others retain that physics deliberately. A neural-ODE method may avoid the classical return-mapping formulation while still requiring numerical integration of a state evolution equation.

### 7.4 Issue 3: Global FEM equilibrium is a separate problem

In implicit nonlinear FEM, the global equilibrium residual is

$$ \mathbf R(\mathbf u) = \mathbf f\_{\mathrm{int}}(\mathbf u) -\mathbf f\_{\mathrm{ext}} =\mathbf0. $$

The internal force depends on the stresses at the material integration points:

$$ \mathbf f\_{\mathrm{int}} = \int\_\Omega \mathbf B^T\boldsymbol\sigma\\,dV. $$

Newton iteration requires a tangent stiffness,

$$ \mathbf K = \frac{\partial\mathbf R}{\partial\mathbf u}, $$

which involves the algorithmic constitutive tangent

$$ \boxed{ \mathbb C\_{\mathrm{alg}} = \frac{\partial\boldsymbol\sigma\_{n+1}} {\partial\boldsymbol\varepsilon\_{n+1}}. } $$

A neural surrogate might produce an accurate stress value but a poor tangent.

This can lead to oscillatory global iterations, excessive substepping, or failure of the equilibrium solve.

For a neural state-update model,

$$ \boldsymbol\sigma\_{n+1} = f\_\theta(\boldsymbol\varepsilon\_{n+1}, \mathbf z\_{n+1}), $$

with

$$ \mathbf z\_{n+1} =\mathcal U\_\theta( \mathbf z_n,\boldsymbol\varepsilon\_{n+1}), $$

the consistent derivative requires the full chain rule:

$$ \boxed{ \mathbb C\_{\mathrm{alg}} = \frac{\partial f\_\theta}{\partial\boldsymbol\varepsilon} + \frac{\partial f\_\theta}{\partial\mathbf z} \frac{\partial\mathcal U\_\theta} {\partial\boldsymbol\varepsilon}. } $$

Ignoring the state-update derivative may produce a tangent inconsistent with the predicted stress.

This issue is particularly important for implicit CP-FEM, where coupling between lattice elasticity, plastic slips and state evolution is strong.

In explicit FEM, there is generally no global Newton equilibrium iteration, but history updates, local stability and the global critical time step still matter.

### 7.5 Issue 4: Trial states must not corrupt committed history

During a global implicit Newton iteration, the solver evaluates multiple trial configurations for the same physical time increment.

A constitutive network should not permanently update its history every time it receives one of those trial strains.

The correct logic is:

Committed state at $t_n$

$\mathbf z_n,\mathbf F^p_n,\mathbf g_n$

Global Newton trial $k$

Compute candidate $\mathbf F^{(k)}\_{n+1}$

Evaluate neural/physical constitutive update from the same committed state

Return trial stress, tangent and provisional updated state

Global equilibrium convergence?

No: reject provisional state and evaluate another trial

Yes: commit $\mathbf z\_{n+1}$ and other history variables

This is a very practical implementation issue. A physically correct neural evolution law can still produce a completely incorrect FE calculation if internal states are inadvertently updated during unsuccessful trial evaluations.

It is one of the reasons standalone neural stress-history benchmarks do not establish compatibility with industrial implicit FE solvers.

## 8. Can the entire history be compressed into a small neural state?

This is where recurrent neural operators and minimal-state networks provide the most compelling evidence.

A general history-dependent material can be written formally as

$$ \boldsymbol\sigma(t) = \mathcal G[ \mathbf F(s),\\,0\leq s\leq t]. $$

The operator depends on a function: the entire previous deformation history.

A recurrent state representation attempts to replace that infinite-dimensional history with a finite-dimensional vector:

$$ \dot{\mathbf z} =g\_\theta(\mathbf F,\mathbf z), \qquad \boldsymbol\sigma=f\_\theta(\mathbf F,\mathbf z). $$

This effectively asks the network to identify a compact state that is sufficient to predict future response.

Liu et al. demonstrate this possibility for crystal-plasticity-based homogenization. Heidenreich's extended minimal state cell approach further investigates how compact recurrent states can reproduce complex rate- and temperature-dependent material responses.

However, a latent state learned from one family of histories may omit variables needed for a new family of histories.

For instance, an RNO trained predominantly on monotonic tension might find a compressed state adequate for work hardening, but insufficient for nonproportional cyclic loading, lattice rotation, twinning or strain-path changes.

This is why the state dimension must be treated as a constitutive model-design choice, not simply a hyperparameter chosen from one test-set error.

A useful formal objective is to minimize predictive error while penalizing model complexity:

$$ \min\_{\theta,k} \left[ \mathcal L\_{\mathrm{history}} +\lambda_k k \right], $$

where $k=\dim(\mathbf z)$.

The reduced model should then be tested for path changes, step-size refinement, microstructure changes and long-term recursive stability.

## 9. Recommended architecture for crystal-plasticity-informed modeling

Based on the reviewed studies, I would not prioritize a conventional PINN that directly learns

$$ (\mathbf x,t,\mathbf F)\mapsto\boldsymbol\sigma. $$

A stronger research direction is a thermodynamically constrained, continuous-time state-space surrogate with explicit elastic structure.

### 9.1 Proposed hybrid model

Proposed CP surrogate architecture

Inputs

$\mathbf F,\dot{\mathbf F},T,\mathbf Q,\mathbf z_n,\boldsymbol\mu$

Physics-preserving elastic branch

Learned free-energy potential with symmetry / objectivity constraints

$\psi\_\theta(\mathbf C^e,\mathbf z,T)$

Continuous internal-state branch

Neural ODE or controlled differential equation for slip-related memory

$\dot{\mathbf z}=g\_\theta(\cdot)$

Thermodynamic / crystallographic constraints

Plastic dissipation · slip-system orientation · objectivity · admissible hardening · history consistency

Constitutive output for FEM

Stress · updated internal state · consistent algorithmic tangent · error indicator

FEM integration with fallback

Use surrogate for admissible states; recover classical CP integration when an error/stability criterion is violated

The physical elastic branch can use

$$ \mathbf C^e=(\mathbf F^e)^T\mathbf F^e, $$

and a learned free-energy potential

$$ \psi\_\theta =\psi\_\theta(\mathbf C^e,\mathbf z,T). $$

An appropriate second Piola-type stress in the elastic intermediate configuration is then related to the potential derivative by

$$ \mathbf S^e =2\frac{\partial\psi\_\theta}{\partial\mathbf C^e}. $$

It must be transformed into the stress measure used by the FE solver.

The internal-state evolution could be represented by

$$ \dot{\mathbf z} = g\_\theta( \mathbf z,\mathbf F,T, \\{\tau^\alpha\\},\boldsymbol\mu). $$

A dissipation condition should be imposed in a form compatible with finite-strain thermodynamics and the chosen internal variables.

For example, using appropriate conjugate plastic forces,

$$ \mathcal D =\mathbf M:\mathbf L^p -\text{stored-hardening-energy rate} \geq0. $$

For a slip-based model, the plastic-power contribution is

$$ \mathbf M:\mathbf L^p =\sum\_\alpha\tau^\alpha\dot\gamma^\alpha. $$

The architecture should not assume that nonnegative slip power alone guarantees the full thermodynamic inequality when stored hardening energy and other internal mechanisms are present.

### 9.2 Suggested loss function

A model combining CP-FEM training data and physical constraints can use

$$ \boxed{ \mathcal L= w\_\sigma\mathcal L\_\sigma +w_z\mathcal L_z +w_D\mathcal L_D +w\_{\mathrm{flow}}\mathcal L\_{\mathrm{flow}} +w\_{\mathrm{ref}}\mathcal L\_{\mathrm{ref}} +w\_{\mathrm{sym}}\mathcal L\_{\mathrm{sym}}. } $$

For example,

$$ \mathcal L\_\sigma = \frac{ \sum_n \\|\boldsymbol\sigma_n^{NN} -\boldsymbol\sigma_n^{CP}\\|\_F^2 }{ \sum_n \\|\boldsymbol\sigma_n^{CP}\\|\_F^2+\epsilon }, $$

$$ \mathcal L_D = \frac1N \sum_n \left[\max(0,-\mathcal D_n)\right]^2, $$

$$ \mathcal L\_{\mathrm{flow}} = \frac1N \sum\_{n,\alpha} \left| \dot\gamma\_{n,\theta}^\alpha -\dot\gamma\_{\mathrm{CP}}^\alpha \right|^2. $$

A useful history-refinement penalty compares identical physical loading trajectories represented with different time subdivisions:

$$ \mathcal L\_{\mathrm{ref}} = \left\\| \mathcal S\_\theta[\mathbf F;\Delta t] - \mathcal S\_\theta[\mathbf F;\Delta t/2] \right\\|^2. $$

This is only meaningful when the compared trajectories represent the same physical history, including the same strain rate for genuinely rate-dependent materials.

Additional terms can enforce frame indifference or material symmetry. For example, a physically objective constitutive model must satisfy the appropriate stress-transformation relationship under superposed rigid rotations.

One should preferably enforce important physical identities through the architecture rather than relying exclusively on finite penalty weights.

### 9.3 Training strategy

I would prioritize the following sequence:

1. Generate stress, state, and tangent data from a validated CP-FEM or CP-FFT reference using multiple microstructures, loading directions, rates and unloading histories.
2. Fit a low-dimensional state-evolution model, starting with the minimum state dimension that accurately captures long-term cyclic and nonproportional loading.
3. Add hard thermodynamic and symmetry constraints, and train stress prediction jointly with history evolution.
4. Test state consistency under time-step refinement and complete recursive rollout, not merely one-step prediction.
5. Deploy in a material-point interface with provisional and committed internal states separated.
6. Construct a consistent tangent by differentiating the complete stress-update map.
7. Use an out-of-domain detector and classical CP integration fallback for uncertain or unstable states.

This hybrid strategy avoids requiring a neural model to extrapolate reliably into completely unknown plasticity regimes.

## 10. A rigorous benchmark suite for practical evaluation

A major weakness of the available literature is that the reported benchmarks are not standardized across architectures.

For a credible comparison between crystal plasticity, RNO, eTANN, INCDE, PRNN and conventional PINNs, I recommend a common suite.

| Benchmark                          | What to compare                                                  | Why it matters              |
| ---------------------------------- | ---------------------------------------------------------------- | --------------------------- |
| Monotonic tension/compression      | Stress, tangent, yield and hardening                             | Baseline accuracy           |
| Cyclic tension–compression         | Hysteresis, Bauschinger effect, energy dissipation               | Loading-memory validation   |
| Nonproportional multiaxial loading | All stress components, slip evolution, latent states             | Tensorial memory            |
| Strain-rate jumps                  | Stress overshoot, relaxation and state update                    | Rate-dependent physics      |
| Temperature changes                | Flow and history dependence                                      | Thermoviscoplastic coupling |
| Loading-path reversal              | Backstress, CRSS and reverse yielding                            | Kinematic hardening         |
| Time-step refinement               | Stress/state variation under halved increments                   | Discretization consistency  |
| Single-element FEM                 | Stress, tangent, convergence and iteration count                 | Initial FE compatibility    |
| Polycrystal RVE                    | Homogenized stress, texture and latent variables                 | Crystal-plasticity validity |
| Full FE boundary-value problem     | Force–displacement, local fields, Newton convergence and runtime | Engineering practicality    |

### 10.1 Recommended error metrics

Use normalized stress error

$$ e\_\sigma = \sqrt{ \frac{\sum_n \\|\boldsymbol\sigma_n^{NN} -\boldsymbol\sigma_n^{ref}\\|\_F^2} {\sum_n\\|\boldsymbol\sigma_n^{ref}\\|\_F^2} }. $$

For stress-history errors,

$$ e\_{\max} = \max_n \\|\boldsymbol\sigma_n^{NN} -\boldsymbol\sigma_n^{ref}\\|\_F. $$

For the cyclic energy response,

$$ e_W = \frac{|W^{NN}-W^{ref}|} {|W^{ref}|+\epsilon}, $$

where $W$ is the accumulated mechanical work or an appropriately defined dissipated-energy measure.

For numerical cost,

$$ S\_{\mathrm{runtime}} =\frac{T\_{\mathrm{reference}}} {T\_{\mathrm{neural}}}. $$

Additionally, report the number of global Newton iterations, local integration failures, rejected increments, numerical substeps and constitutive-tangent evaluations.

Most importantly, compare CPU-only inference, GPU inference, and the total cost of training-data generation separately.

### 10.2 Proposed acceptance criteria

There is no universally accepted stress-error threshold for all plasticity simulations. A model accurate enough for a smooth tensile response may be unsuitable for localization, fatigue, contact or failure prediction.

Therefore, validation tolerances should be set by the engineering quantity of interest. For example, an application may tolerate several percent stress error but almost no error in the strain at which instability begins.

For safety-critical or high-cycle applications, long-term error accumulation and robustness are often more important than one-step mean-square error.

## 11. Where the field still has unresolved problems

Several research gaps remain even in the most recent publications.

Reliable extrapolation beyond the training domain. A learned constitutive law may obey thermodynamics yet be inaccurate for unseen loading orientations, temperatures or slip mechanisms.

Loss of physical detail during homogenization. A two- or five-variable latent state can reproduce macroscopic stress while missing local stress concentrations, shear bands, texture evolution or grain-level damage.

Identifiability of internal variables. Different combinations of hidden variables and potentials may predict the same observed stress histories. Flaschel et al. explicitly discuss nonuniqueness of learned thermodynamic potentials.

Consistent FE tangents. Stress accuracy is insufficient unless the derivative of the complete state update is suitable for equilibrium iteration.

Cost of generating training data. CP-FEM and CP-FFT trajectories may remain the most expensive part of the full learning workflow, particularly for parameterized microstructures.

Stiff history dynamics. Neural ODEs reduce dependence on arbitrary time increments but do not eliminate numerical stiffness, substepping or integration error.

Rate independence versus physical time dependence. Plasticity driven by a loading path and viscoplasticity driven by time and rate have different mathematical requirements. The same discretization-independent architecture cannot automatically guarantee both without appropriate evolution structure.

Limited like-for-like hardware benchmarks. Many papers compare neural inference against an expensive multiscale reference, without separately reporting data generation, training, FE tangent evaluation and full structural solve costs.

The 2026 comparison by [Jones and Fuhg](https://scispace.com/papers/a-hierarchy-of-thermodynamics-learning-frameworks-for-e64fizia6qqx?utm_source=chatgpt.com) is especially relevant here: it examines how different thermodynamic assumptions affect expressiveness, stability and learnability for inelastic datasets, including a rate-dependent crystal-plasticity polycrystal.

## 12. Overall research assessment

The following is a qualitative synthesis of available validation evidence, not a head-to-head numerical ranking.

| Approach                                | History handling                         | Practical advantage                                        | Main limitation                                   |
| --------------------------------------- | ---------------------------------------- | ---------------------------------------------------------- | ------------------------------------------------- |
| Standard PINN without state evolution   | Weak for irreversible histories          | Useful for simpler field/inverse problems                  | Incomplete constitutive memory                    |
| LSTM / GRU                              | Learned recurrent memory                 | Simple stress-history surrogate                            | Discretization and rollout stability              |
| TANN / generalized standard material NN | Thermodynamic internal variables         | Strong mechanical constraints                              | Potentially expensive derivatives and integration |
| eTANN                                   | Continuous-time internal-state ODE       | Better time-step adaptability                              | ODE solver cost and state identification          |
| INCDE                                   | Controlled incremental history evolution | FE-oriented stability formulation                          | Needs broader CP and industrial validation        |
| PRNN                                    | Embedded physical constitutive states    | Strong extrapolation potential and multiscale acceleration | Retains local constitutive solvers                |
| RNO / CP-RNN                            | Compressed learned internal states       | Efficient history-dependent polycrystal surrogate          | Latent-state and tangent robustness               |
| Direct return-mapping surrogate         | Explicit previous-state input            | Low local inference cost                                   | Risk near yield transitions and novel loading     |
| Neural operator plus solver correction  | History-aware spatial prediction         | Speed with correction-based accuracy                       | Classical solver still participates               |

### Which methods are most promising?

For viscoplasticity with temperature and rate history, eTANN and compact continuous-time recurrent models provide an attractive mathematical formulation. Physics-informed viscoplastic potentials are valuable when thermodynamic admissibility is a priority.

For cyclic elastoplasticity, explicitly propagated plastic strains and backstresses with physics-informed constraints are preferable to a completely opaque stress-only recurrent network when the material's hardening structure is already understood.

For crystal plasticity and multiscale FE calculations, the strongest starting points are RNO, CP-RNN and PRNN. Their different treatment of internal states offers a useful trade-off between pure learned acceleration and preservation of physical constitutive laws.

For replacing expensive local predictor/corrector solves, methods such as No-Prop are conceptually direct. However, their success on a rate-independent nonlinear-hardening model should not be treated as validation for the much more complex multislip crystal-plasticity problem.

For whole-structure acceleration, physics-informed neural operators and solver warm starts are promising because they address more than constitutive evaluation. The 2026 Plasolver work offers an example, although its preprint results need further application-specific benchmarking.

## 13. Conclusion

Mechanics-informed neural networks are becoming practical for elastoplastic, viscoplastic and crystal-plasticity modeling, particularly when replacing expensive microscopic simulations or complex history-dependent constitutive updates. Reported gains range from substantial training-data reductions to approximately 1,000× faster local homogenized stress evaluations, but comparable full-FEM hardware and timing evidence remains limited.

The most promising resolution to crystal plasticity's history-dependence problem is continuous-time internal-state learning combined with mechanically consistent stress prediction, stable state integration, and a consistent FEM tangent. These methods can replace expensive local plasticity calculations or compress complex microscopic histories, but they do not generally eliminate the need to solve structural equilibrium. For a robust research implementation, a hybrid eTANN/RNO-type surrogate with explicit physical constraints, state-commit discipline and fallback to conventional CP integration offers a strong balance between computational efficiency and reliability.

### 🔍 Ready to go beyond a quick answer?

#### 🚀 Plan and execute research workflows via SciSpace AI Agent

1\. Compare mechanics-informed plasticity models\
▸ gather ▸ compare ▸ document • 📝 report • ⏩ [Run Task](https://scispace.com?utm_campaign=scispacegpt\&q=Compare+mechanics-informed+plasticity+models.+gather,+compare+and+document.+deliver+a+report.\&utm_source=chatgpt.com)

2\. Benchmark crystal plasticity surrogate runtimes\
▸ collect ▸ measure ▸ compare • 📊 chart • ⏩ [Run Task](https://scispace.com?utm_campaign=scispacegpt\&q=Benchmark+crystal+plasticity+surrogate+runtimes.+collect,+measure+and+compare.+deliver+a+chart.\&utm_source=chatgpt.com)

3\. Design a history-aware crystal plasticity surrogate\
▸ specify ▸ prototype ▸ validate • 📄 draft • ⏩ [Run Task](https://scispace.com?utm_campaign=scispacegpt\&q=Design+a+history-aware+crystal+plasticity+surrogate.+specify,+prototype+and+validate.+deliver+a+draft.\&utm_source=chatgpt.com)

Try [SciSpace](https://scispace.com?utm_campaign=scispacegpt\&utm_source=chatgpt.com) – the leading AI platform for academic research. It's your all-in-one workspace to discover, read, and analyze scientific literature.