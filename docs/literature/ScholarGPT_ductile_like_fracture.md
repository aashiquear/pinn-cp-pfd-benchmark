# **Literature review and modeling assessment: Orientation-dependent ductile fracture in single-crystal nickel**

**Scope:** Finite-strain FCC crystal plasticity coupled to an AT2 phase-field fracture model, with particular emphasis on orientation-dependent tensile failure, plastic-work driving forces, dislocation storage, and regularization effects.

Literature verification incomplete

**Evidence limitation.** Live scholarly retrieval is unavailable in this session. Consequently, I cannot responsibly claim to have verified the cited papers, their exact governing equations, their DOI records, or their experimental validation datasets. I will distinguish established mathematical consequences of your supplied model from candidate constitutive modifications requiring literature verification. I will not fabricate a bibliography or convert unverified recollections into citations.

## **1. Summary of findings**

Your observations expose two different potential problems: **the fracture model may use the wrong internal variable to determine damage initiation**, or **the regularization length may be sufficiently large to alter the structural failure mechanism**. Both deserve investigation before assigning orientation-dependent fracture parameters.

The central findings from analyzing your equations are:

1\. **A universal plastic-work criterion cannot distinguish orientations on the basis of their total work alone.** If the ductile-like orientation accumulates more crack-driving work at every corresponding strain, a common monotonically increasing damage law with identical resistance cannot selectively delay its damage under otherwise equivalent mechanical conditions.

2\. **Total plastic work is not the same as stored energy available to create a crack.** A model using all accumulated plastic work as a fracture driver effectively assigns the entire plastic-work history to that driver. A dislocation-storage-based variable or a state-dependent stored-energy fraction provides a physically motivated route to distinguishing orientations with different slip activity.

3\. **A plastic-strain-controlled initiation threshold can decouple fracture initiation from flow stress.** This is especially relevant if two orientations have comparable accumulated plastic strain but very different accumulated plastic work. However, a universal threshold does not automatically reproduce orientation-dependent fracture strains.

4\. **Crystallographic fracture resistance can distinguish orientations even when their flow stresses are similar.** An orientation-dependent crack-plane energy, or a slip-system-resolved damage formulation, can represent the influence of active slip systems and available crack planes. It requires evidence that material separation is crystallographically controlled.

5\. **A local axial strain near 1.0 does not establish an intrinsic fracture strain of 1.0.** A localized shear or necking band may form through plastic instability and subsequently undergo material separation. Damage initiation, localization, and final specimen separation must be treated as distinct events.

6\. **Your length scale is a major confounding variable.** With \\(l_0=0.6\\) mm and a gauge thickness of \\(0.35\\) mm,

\\[ \boxed{\frac{l_0}{t}=1.71} \\]

The phase-field regularization length exceeds the gauge thickness. This does not mathematically invalidate the simulation, but it makes the interpretation of \\(G_c\\), the predicted damage-zone geometry, and the sensitivity to thickness especially important.

7\. **A constant plastic-work multiplier is not, strictly speaking, identical to rescaling&#x20;**\\(G_c\\)**&#x20;in your full model.** They become approximately interchangeable in a plastic-work-dominated regime with fixed geometry and \\(l_0\\). The tensile elastic contribution and coupled evolution prevent exact equivalence.

### **Preliminary ranking**

For your particular observations, I would prioritize a **slip-state-dependent stored-energy driving force**, followed by a **plastic-strain-controlled initiation gate**, with a **stress-state or crystallographic separation mechanism** chosen according to the experimental fracture morphology. Before calibrating any of these, a length-scale sensitivity study is essential.

## **2. Reference model and mathematical diagnosis**

Your crack-driving history variable is

\\[ H(\mathbf X,t)=\max\_{s\leq t} \left[ \psi_e^+(\mathbf X,s)+W_p(\mathbf X,s) \right], \\]

with

\\[ W_p(t)= \int_0^t\sum\_{\alpha=1}^{12} \tau^\alpha\dot\gamma^\alpha\\,\mathrm dt . \\]

The damage equilibrium equation is

\\[ \frac{G_c}{l_0}d-G_cl_0\nabla^2d =2(1-d)H. \tag{1} \\]

Here and below, the expressions use your supplied convention for the tensile elastic energy and plastic-work accumulation.

### **2.1 Why plastic work can penalize a high-strength orientation**

If plastic work dominates the driving energy, then approximately

\\[ W_p\sim\int \sigma\_{\mathrm{flow}}\\,\mathrm d\varepsilon_p. \\]

At comparable accumulated plastic strains, the orientation with higher flow stress can have greater \\(W_p\\), even if it is experimentally more ductile.

In the spatially homogeneous limit, \\(\nabla d=0\\), equation (1) gives

\\[ d=\frac{2H}{G_c/l_0+2H}. \tag{2} \\]

For fixed \\(G_c\\) and \\(l_0\\),

\\[ \frac{\partial d}{\partial H}>0. \\]

Thus a higher driving energy produces greater equilibrium damage in this simplified comparison.

This is a meaningful limitation of the constitutive choice, although it is not a proof that the higher-work specimen must always fail first in a fully coupled, nonuniform boundary-value problem. Localization, lattice rotation, damage redistribution, and the specimen's thickness profile can change the ordering.

### **2.2 An important distinction concerning&#x20;**\\(\zeta\\)**&#x20;and&#x20;**\\(G_c\\)

Divide equation (1) by \\(G_c/l_0\\):

\\[ d-l_0^2\nabla^2d = \frac{2l_0}{G_c}(1-d) \left(\psi_e^++\zeta W_p\right). \tag{3} \\]

In the limiting case \\(\psi_e^+\ll\zeta W_p\\),

\\[ \frac{2l_0\zeta}{G_c}W_p \\]

controls the dimensionless driving term. Therefore, changing a global \\(\zeta\\) is approximately equivalent to changing \\(G_c\\) inversely, as far as this term is concerned.

But when the elastic term is appreciable,

\\[ \frac{\psi_e^+}{G_c} +\frac{\zeta W_p}{G_c} \\]

contains two different ratios. A change in \\(\zeta\\) affects one and not the other.

**The more important identifiability problem is that neither a universal&#x20;**\\(\zeta\\)**&#x20;nor a universal&#x20;**\\(G_c\\)**&#x20;necessarily supplies the missing orientation selectivity.**

### **2.3 Why AT2 is particularly relevant**

Your model has no finite elastic damage-initiation stage.

From equation (2), for any positive \\(H\\),

\\[ d>0. \\]

That is intrinsic to the quadratic AT2 crack-density formulation under this type of driving force. The calculated damage can begin accumulating before any experimentally observable fracture event.

Consequently, one must distinguish:

- the appearance of small diffuse \\(d\\);
- acceleration of damage in a localization band;
- peak load and structural instability;
- loss of load-carrying capacity;
- final physical separation.

These are not necessarily coincident. Moving the last event by changing a damage parameter does not necessarily correct the physical onset mechanism.

## **3. Candidate modifications of crack driving force and fracture resistance**

The following table is a **constitutive comparison, not a verified literature attribution table**. Its equations are schematic implementations or direct adaptations of the model you supplied. They are intended to identify mechanisms and enable controlled numerical tests, not to represent exact formulations from particular papers.

Define

\\[ p(t)=\int_0^t\dot{\bar\varepsilon}\_p\\,\mathrm dt,\qquad P\_\alpha(t)=\int_0^t\tau^\alpha\dot\gamma^\alpha\\,\mathrm dt, \\]

and let \\(\rho^\alpha\\) denote slip-system dislocation density, \\(\mathbf n\\) a prospective crack-plane normal, and \\(\mathbf Q\\) the crystal orientation.

For compactness, \\(H=\max_t[\psi_e^++D_p]\\), unless a separate evolution law is indicated. The accumulated scalar plastic strain \\(p\\) must be defined consistently with finite-strain crystal plasticity; it is not automatically equivalent to the sum of absolute slips.

| Modification                                        | Representative mathematical form                                                                                                                                                              | New parameters and physical meaning                                                                                                    | Crystal plasticity application                                                   | Single-crystal tensile validation                                               |
| --------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| **1. Plastic-work prefactor**                       | \\(D_p=\zeta W_p\\)                                                                                                                                                                           | \\(\zeta\\): fraction or weighting of accumulated plastic work contributing to damage                                                  | Already present in your CP model                                                 | Your benchmark comparison is supplied, but independent validation is unverified |
| **2. Plastic-work threshold**                       | \\(D_p=\zeta\langle W_p-W_0\rangle\_+\\)                                                                                                                                                      | \\(W_0\\): accumulated work required before plastic work contributes                                                                   | Compatible; published application unverified                                     | Unverified                                                                      |
| **3. Plastic-strain initiation threshold**          | \\(D_p=\zeta\int_0^t \mathcal I[p(s)\geq p_0]\\,\mathrm dW_p(s)\\)                                                                                                                            | \\(p_0\\): accumulated plastic strain at activation                                                                                    | Compatible; published application unverified                                     | Unverified                                                                      |
| **4. Stored fraction of plastic work**              | \\(D_p=U\_{\rm store}=\int_0^t\chi(s)\\,\mathrm dW_p(s)\\)                                                                                                                                    | \\(\chi\\): stored fraction; \\(\beta=1-\chi\\) as a simplified Taylor–Quinney interpretation                                          | Compatible; specific application unverified                                      | Unverified                                                                      |
| **5. Dislocation-density stored energy**            | \\(D_p=\langle U\_\rho-U\_{\rho,c}\rangle\_+\\), with \\(U\_\rho\propto\mu b^2\sum\_\alpha\rho^\alpha\\) in a simplified dislocation-energy estimate                                          | Dislocation-energy coefficient, interaction assumptions, and \\(U\_{\rho,c}\\): energetic storage and activation                       | Directly compatible with dislocation-density CP; published validation unverified | Unverified                                                                      |
| **6. Plastic-strain-dependent degradation**         | \\(g(d,p)=(1-d)^2q(p)\\), with \\(q(p)>0\\)                                                                                                                                                   | Parameters defining \\(q\\): additional mechanical softening associated with plastic strain                                            | Mathematically compatible; thermodynamic treatment required                      | Unverified                                                                      |
| **7. Plastic-strain-dependent fracture resistance** | \\(G_c=G\_{c0}F(p)\\), often considering \\(F'(p)<0\\) for degradation                                                                                                                        | Softening rate and characteristic strain: progressive toughness reduction                                                              | Compatible; published CP application unverified                                  | Unverified                                                                      |
| **8. Triaxiality/Lode-dependent damage**            | \\(D_p=\int_0^t f(T,L)\\,\mathrm dW_p\\), \\(T=\sigma_m/\sigma\_{\rm eq}\\)                                                                                                                   | Parameters in \\(f\\): dependence on hydrostatic tension and deviatoric stress state                                                   | Compatible with CP stress histories; specific validation unverified              | Unverified                                                                      |
| **9. Anisotropic crack resistance**                 | \\(G_c=G_c(\mathbf n,\mathbf Q)\\)                                                                                                                                                            | Directional toughness or fracture-energy coefficients: preferred crystallographic fracture planes                                      | Compatible with CP orientation updates                                           | Unverified                                                                      |
| **10. Slip-plane-resolved phase fields**            | \\(\displaystyle \mathcal G\_\Gamma=\sum\_\alpha \frac{G_c^\alpha}{c\_\alpha}\int\_\Omega\left[\frac{w\_\alpha(d\_\alpha)}{l\_\alpha}+l\_\alpha\|\nabla d\_\alpha\|^2\right]\mathrm dV\\)     | \\(d\_\alpha,G_c^\alpha,l\_\alpha\\): separate damage mechanisms or candidate slip/crack planes; \\(c\_\alpha\\) denotes normalization | Structurally suited to CP; published applications unverified                     | Unverified                                                                      |
| **11. Void growth and porosity**                    | \\(\dot f_v=(1-f_v)\operatorname{tr}\mathbf D_p+\dot f\_{\rm nuc}\\), with suitable restrictions on the growth term                                                                           | Initial porosity, nucleation, coalescence and stress-state parameters: internal void evolution                                         | Can be coupled to CP; specific formulations unverified                           | Unverified                                                                      |
| **12. AT1 / cohesive phase-field formulation**      | \\(\displaystyle\mathcal G\_\Gamma=\frac{G_c}{c_w}\int\_\Omega\left[\frac{w(d)}{l_0}+l_0\|\nabla d\|^2\right]\mathrm dV\\), with \\(w(d)=d\\) for AT1; PF-CZM uses model-specific degradation | Fracture strength, critical initiation energy and cohesive softening parameters                                                        | Compatible in principle; CP applications unverified                              | Unverified                                                                      |

Notation and limitations: \\(\langle x\rangle\_+=\max(x,0)\\). The porosity growth expression represents a simplified plastically incompressible matrix kinematic setting, not a full void-nucleation/coalescence model. The multiphase-field equation is an energy-structure illustration, not a complete constitutive formulation. In particular, slip-plane fracture requires suitable directional gradient terms, coupling, and degradation rules. No row above is asserted to reproduce an individually verified publication.

### **3.1 Which mechanisms can actually reverse the orientation trend?**

This is the most useful way to classify the modifications.

| Modification                                    | Can it selectively delay a higher-work orientation? | Qualification                                                                                       |
| ----------------------------------------------- | --------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| One common \\(\zeta\\)                          | Generally no                                        | It applies the same work weighting across orientations                                              |
| One common \\(W_0\\)                            | Generally no                                        | A higher-work specimen crosses the same work threshold no later, if its local histories are ordered |
| One common plastic-strain threshold \\(p_0\\)   | **Potentially**                                     | Works only if the ordering of plastic strain differs meaningfully from the ordering of plastic work |
| Variable stored-work fraction \\(\chi\\)        | **Yes**                                             | Requires a physically justified difference in storage rather than arbitrary orientation tuning      |
| Dislocation-energy criterion                    | **Yes**                                             | High accumulated work need not imply equally high stored dislocation energy                         |
| \\(G_c(p)\\) decreasing with plastic strain     | Generally no                                        | If the ductile-like orientation also accumulates more \\(p\\), it is weakened earlier               |
| Triaxiality/Lode weighting                      | **Yes, conditionally**                              | Requires distinct local stress-state histories                                                      |
| Crystallographic \\(G_c(\mathbf n,\mathbf Q)\\) | **Yes**                                             | Requires differences in accessible fracture planes or directional resistance                        |
| Slip-resolved damage                            | **Yes**                                             | Can distinguish the distribution of plastic activity across systems                                 |
| Porosity coupling                               | **Yes, conditionally**                              | Requires distinct void-growth/nucleation conditions                                                 |
| AT1/PF-CZM threshold                            | **Potentially**                                     | Improves initiation behavior but does not by itself guarantee reversed orientation ordering         |

This classification separates a model that adds a physical distinguishing variable from one that simply adds another shared fitting parameter.

### **3.2 The issue with Taylor–Quinney fractions**

A useful energy balance is

\\[ W_p=Q+U\_{\rm store}, \\]

where \\(Q\\) is the plastic work converted to heat and \\(U\_{\rm store}\\) is energy retained in the evolving microstructure, under the simplified partition assumed here.

For a constant differential heat fraction \\(\beta\\),

\\[ \mathrm dU\_{\rm store}=(1-\beta)\\,\mathrm dW_p. \\]

An associated candidate crack-driving term is therefore

\\[ D_p=\int_0^t(1-\beta(s))\\,\mathrm dW_p(s). \tag{4} \\]

The distinction between differential and cumulative fractions matters.

If \\(\beta\\) depends on dislocation density, plastic strain, or slip-system activity, it generally is not correct to calculate stored energy simply as

\\[ U\_{\rm store}(t)=[1-\beta(t)]W_p(t). \\]

The correct history integral depends on the fraction at the time each increment of work occurs.

Further, the full accumulated plastic work cannot automatically be interpreted as energy available for fracture. Your current \\(W_p\\) is an accumulated dissipation-related quantity; promoting it to a fracture-driving history is a modeling choice that requires a consistent energy interpretation.

**Critical caution:** Replacing \\(W_p\\) with stored energy may improve orientation discrimination, but it is not established that stored dislocation energy alone supplies the energy of fracture. The relationship between storage, damage, dissipation, and material separation must be specified.

### **3.3 Dislocation-density criteria are particularly testable in your framework**

Your existing Kocks–Mecking hardening law already evolves dislocation-density-related variables.

A simplified dislocation energy measure can be expressed schematically as

\\[ U\_\rho \sim C\_\rho\mu b^2\sum\_{\alpha=1}^{12}\rho^\alpha. \tag{5} \\]

Here \\(\mu b^2\rho\\) has the correct dimensions of energy per volume, and \\(C\_\rho\\) accounts for dislocation line-energy and related assumptions. Its numerical value should not be assumed without an appropriate model.

A relevant quantity for fracture would ordinarily be an energy change relative to the initial microstructure, not an arbitrary absolute energy reference.

An especially informative comparison is

\\[ r\_{\rm store}(t)= \frac{\Delta U\_\rho(t)}{W_p(t)}. \tag{6} \\]

If the two prematurely failing orientations have high \\(W_p\\) but substantially lower \\(\Delta U\_\rho/W_p\\), that supports your basic hypothesis.

If they instead accumulate both greater plastic work and greater stored energy at corresponding strains, simply replacing \\(W_p\\) by \\(U\_\rho\\) may not solve the problem. A different fracture threshold or stress-state mechanism would then be needed.

Dislocation density is also not necessarily a unique proxy for stored energy. Dislocation interactions, character, arrangement, recovery, and gradients may be important.

### **3.4 Plastic-strain gates and fracture-energy degradation are not interchangeable**

Consider two distinct changes.

**Strain-gated initiation:**

\\[ D_p(t)=\int_0^t \mathcal I[p(s)\geq p_0]\\,\mathrm dW_p(s). \tag{7} \\]

Plastic-work-driven damage cannot activate through this term until a chosen plastic strain is accumulated.

**Degrading toughness:**

\\[ G_c(p)=G\_{c0}F(p),\qquad F'(p)<0. \tag{8} \\]

The latter reduces material resistance as plastic strain grows. It does not establish a protected initial plastic stage unless some additional initiation structure does so.

For two orientations that accumulate similar plastic strain but different work, equation (7) can delay the effect of the high flow stress before activation. Equation (8), by contrast, may fail to provide the required selectivity because both specimens experience similar resistance reduction.

One further implementation detail is important: if \\(G_c\\) depends on a spatially varying field, the simple constant-\\(G_c\\) Laplacian form of equation (1) is no longer generally the complete variational derivative. Additional terms arise from the spatial variation of resistance, and potentially from its coupling to plasticity.

### **3.5 AT1 or PF-CZM formulations address a separate defect**

The quadratic AT2 crack-density function has zero slope at zero damage. Under conventional energetic driving, that is associated with no strictly finite undamaged activation threshold.

AT1 instead uses a crack-density function with nonzero slope at zero damage, which permits a finite elastic damage-initiation threshold in a suitable variational formulation.

A phase-field cohesive-zone model can introduce an independently controlled material strength and a cohesive separation response, with details depending on the chosen degradation function.

These formulations may be valuable if your model develops damage much earlier than observed. But they do not, by themselves, explain why a higher-strength orientation should tolerate substantially more plastic deformation before fracture.

They improve the definition of initiation; they do not automatically introduce the necessary crystallographic selectivity.

## **4. Experimental evidence: FCC single-crystal tensile ductility**

### **4.1 What can be inferred from single-crystal plasticity**

Crystal orientation affects the resolved shear stresses on the 12 FCC slip systems:

\\[ \tau^\alpha=\boldsymbol{\sigma}:\mathbf P^\alpha, \\]

where \\(\mathbf P^\alpha\\) is the appropriate Schmid tensor in the current configuration.

Different initial orientations can therefore produce different combinations of slip activity, work hardening, crystal rotation, and strain localization.

The following mechanisms should be distinguished.

| Mechanism                             | Expected mechanical signature                                                                                 | Relevant evidence in your experiments                        |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| Predominantly single-slip deformation | Concentration of shear on a limited number of systems; potentially strongly directional surface slip patterns | Slip traces, local shear displacement, EBSD rotation         |
| Multiple-slip deformation             | More distributed crystallographic activity and potentially different hardening                                | Slip-system activity, hardening rate, lattice rotation       |
| Geometric softening                   | Rotation changes Schmid factors and may increase slip activity under the applied loading                      | Evolution of crystal orientation and tangent hardening       |
| Diffuse or localized necking          | Loss of stability of the tensile deformation; strain concentration over a finite region                       | Load maximum, DIC strain field, thickness or area reduction  |
| Shear-band localization               | Concentrated shear and rotation within an inclined zone                                                       | DIC displacement discontinuity tendency, EBSD gradients      |
| Crystallographic material separation  | Fracture morphology follows preferred crystallographic planes                                                 | Fractography and crystal-frame indexing of fracture surfaces |
| Void growth and coalescence           | Progressive separation through cavity growth                                                                  | Fractography, microscopy, porosity characterization          |

These are mechanically distinct mechanisms, even though several can occur in the same specimen.

I cannot verify from accessible literature whether a particular ordering of tensile ductility among nickel single-crystal orientations is universal, nor whether Rezwan et al. observed each mechanism. Consequently, I would not assign a fixed *single-slip = ductile* or *multiple-slip = brittle* rule without experimental support.

### **4.2 Is failure in thin specimens a plastic instability or material separation?**

This is not necessarily an either-or question.

A possible sequence is

**Homogeneous crystal-plastic deformation**

Slip activity, hardening and lattice rotation

**Loss of deformation stability**

Geometric softening, necking or shear localization

**Severe local plastic deformation**

Large strain and stress-state changes in a localized band

**Material separation**

Crystallographic cracking, void coalescence, or another separation mechanism

This is a possible mechanistic sequence, not a claim that it has been established experimentally for your samples.

The particularly important observation is the local axial strain approaching 1.0 in your ductile-like specimens.

That observation suggests substantial deformation concentration, but does not tell us whether the local strain reached 1.0 before damage initiation or partly as a result of an already developing fracture process.

To answer the instability-versus-separation question, I would examine three events separately:

1. First sustained deviation from distributed plastic flow.
2. Load maximum and rapid growth of a localization band.
3. First observable crack followed by complete separation.

If the first two events are accurately predicted by crystal plasticity without damage, while the coupled model separates prematurely, the fracture law is the likely primary problem.

If even the undamaged crystal-plasticity model fails to reproduce the band, the flow rule, hardening, lattice rotation or geometric boundary-value problem must be revisited first.

**This distinction is critical:** the global engineering strain at final fracture is not automatically a transferable constitutive material parameter.

### **4.3 What orientation-dependent validation should demonstrate**

A convincing validation against single-crystal tensile specimens of different orientations should reproduce more than final strain:

- orientation-dependent force–displacement or stress–strain curves;
- onset and spatial structure of localization;
- lattice rotation and crystallographic deformation bands;
- failure location and observed fracture-plane orientation;
- failure strain across the relevant thickness profiles.

A paper matching only macroscopic failure strains would not establish that its orientation-dependent crack-driving mechanism is physically correct.

I cannot identify a verified publication in this session that demonstrates all these requirements for your particular nickel benchmark. That remains an open literature-verification task rather than a finding that such work does not exist.

## **5. Ranked recommendations for your case**

These recommendations are based on the structure of your constitutive model and the observations you supplied, not on an unverified claim of published superiority.

## **1**

**Stored dislocation energy or slip-state-dependent stored-work fraction**

First constitutive test

This directly addresses your central concern: high flow stress produces high accumulated plastic work, but need not produce a correspondingly high amount of energy stored in the microstructure.

## **2**

**Plastic-strain-controlled initiation with plastic work retained as the subsequent driver**

This decouples the initiation criterion from flow stress while retaining plastic-work sensitivity after activation. It is the simplest controlled test of whether the wrong variable is triggering damage.

## **3**

**Crystallographic fracture resistance or slip-plane-resolved damage**

This is attractive if the experimental fracture surfaces are strongly crystallographic. It can distinguish orientations through the relationship between loading direction, active slip systems and admissible fracture planes.

I would replace recommendation 3 with a triaxiality-/porosity-dependent model if experimental fractography instead demonstrates void-driven separation following strong necking.

### **5.1 Recommendation 1: Test dislocation stored energy**

**Reasoning.** Your model already contains dislocation-density evolution, making this a relatively economical modification to investigate.

The key question is whether the ductile-like orientation converts a smaller fraction of plastic work into stored energy.

#### **Minimal two-specimen test**

Choose one brittle-like specimen \\(B\\) and one ductile-like specimen \\(D\\), preferably with closely matched thickness profiles.

First run both with damage disabled, or with fracture strongly suppressed. Save the following local histories within the eventual failure region:

\\[ W_p(t),\quad p(t),\quad \rho^\alpha(t),\quad \gamma^\alpha(t),\quad \psi_e^+(t). \\]

Compute a consistent dislocation stored-energy proxy and the ratios

\\[ R_B=\frac{\Delta U\_{\rho,B}}{W\_{p,B}}, \qquad R_D=\frac{\Delta U\_{\rho,D}}{W\_{p,D}}. \\]

Compare both at common global strain and at corresponding stages of localization.

**Support for the hypothesis:** \\(D\\) accumulates high total plastic work but substantially less stored-energy-related crack-driving energy than \\(B\\), and a single energy-based criterion reproduces the observed failure ordering.

**Evidence against it:** The stored-energy histories have essentially the same orientation ordering as the original plastic-work histories.

A successful model should use one common constitutive rule, with differences arising from the simulated state, not separately fitted storage coefficients for each orientation.

### **5.2 Recommendation 2: Test a plastic-strain initiation gate**

**Reasoning.** A high-strength orientation should not automatically be penalized for accumulating more work at the same plastic strain.

For a simple numerical experiment, consider the strain-gated accumulated work defined in equation (7).

This means that plastic work accumulated before \\(p_0\\) is not retroactively added to the driving term. That distinction matters: a model that switches on all previously accumulated work at \\(p=p_0\\) could produce an abrupt damage-driving jump.

#### **Minimal two-specimen test**

Use the same two specimens and preserve all existing crystal-plasticity parameters.

Set one common \\(p_0\\). Sweep it from zero through a range that includes the plastic strains at which premature damage develops in the baseline.

For each candidate threshold, record:

- force–displacement response;
- first strongly localized phase-field growth;
- strain at peak load;
- global strain at final loss of load capacity;
- local plastic strain at failure.

**Support:** One common \\(p_0\\) preserves the already acceptable \\(B\\) result while delaying \\(D\\) toward the observed failure range, without an artificial load jump.

**Evidence against:** Thresholds large enough to correct \\(D\\) consistently postpone \\(B\\) beyond its experimental failure strain.

That would demonstrate that a single global strain threshold cannot adequately separate the two classes.

An important additional diagnostic is whether the plastic strain field at the premature predicted failure is already larger in \\(D\\) than in \\(B\\). If so, a common plastic-strain threshold is less likely to help.

### **5.3 Recommendation 3: Test crystallographic fracture resistance**

**Reasoning.** A scalar isotropic \\(G_c\\) treats every crack orientation equally in the crystal frame. It can produce directionally different failure indirectly through anisotropic crystal plasticity, but has no explicit crystallographic crack-plane resistance.

If the two groups fail through different crystallographic mechanisms, this omission could matter.

#### **Minimal two-specimen test**

Start with experimental fracture-plane indexing.

For each specimen, identify the orientation of the observed fracture plane relative to the crystal lattice, accounting for lattice rotation if substantial.

Then:

1. Keep the crystal-plasticity parameters unchanged.
2. Replace the isotropic resistance with a common crystal-frame directional resistance law.
3. Rotate that law with the evolving crystal orientation.
4. Use the same directional material parameters for both specimens.
5. Compare fracture onset, crack path, and final strain.

**Support:** The same crystal-frame law predicts distinct fracture behavior and the correct fracture-plane orientations in both specimens.

**Evidence against:** The observed fracture surfaces do not exhibit a consistent crystallographic relationship, or correcting the failure strains requires separate arbitrary toughness values for the specimens.

This experiment may require a directional crack-surface functional, not merely substituting a direction-dependent scalar into equation (1). A variationally consistent anisotropic gradient term is important.

## **6. A prerequisite test: Is&#x20;**\\(l_0=0.6\\)**&#x20;mm too large?**

Regularization length

# **0.60 mm**

Gauge thickness

# **0.35 mm**

To scale. The phase-field length is approximately 1.71 times the gauge thickness.

### **6.1 Why this ratio matters**

The regularization length influences the spatial distribution of damage and, in AT2 formulations, the intrinsic strength scale.

For a conventional one-dimensional elastic AT2 reference problem, the characteristic tensile strength scales as

\\[ \sigma_c\propto\sqrt{\frac{E G_c}{l_0}}. \tag{9} \\]

This is a reference scaling, not an exact strength formula for your finite-strain crystal-plasticity problem.

It shows why a change in \\(l_0\\) can change predicted fracture strain even if all crystal-plasticity parameters remain unchanged.

When the regularization length exceeds the physical thickness, the diffuse damage region may be strongly constrained by the specimen geometry.

That can influence:

- interaction of damage with free surfaces;
- the spatial extent of the damaged band;
- sensitivity to thickness-profile changes;
- competition between plastic localization and fracture;
- the relationship between nominal \\(G_c\\) and measured structural failure.

It is therefore possible that the two problematic orientations interact with the regularized damage process differently because of their distinct localization patterns.

### **6.2 Mesh convergence is not the same as length-scale convergence**

Two independent checks are needed.

**Mesh convergence at fixed&#x20;**\\(l_0\\)**.**

Reduce the element size while holding material parameters and the phase-field length fixed. Check convergence of force–displacement response, damage distribution, failure location, and failure strain.

**Length-scale sensitivity.**

Vary \\(l_0\\) using sufficiently refined meshes for each chosen value. Compare the resulting orientation-dependent failure strains.

There is a subtle identifiability issue here: reducing \\(l_0\\) without recalibration changes the intrinsic phase-field strength. For a meaningful assessment, report both the direct sensitivity and a comparison in which the calibration specimen is refitted consistently.

### **6.3 Minimal numerical study**

I would use the following diagnostic matrix for the same \\(B\\) and \\(D\\) specimens:

| Run | Length scale | Calibration strategy        | Main diagnostic                          |
| --- | ------------ | --------------------------- | ---------------------------------------- |
| A   | 0.60 mm      | Current baseline            | Reproduce original mismatch              |
| B   | 0.30 mm      | Existing parameters         | Direct length sensitivity                |
| C   | 0.15 mm      | Existing parameters         | Damage-zone and localization sensitivity |
| D   | 0.30 mm      | Recalibrate against \\(B\\) | Transferability to \\(D\\)               |
| E   | 0.15 mm      | Recalibrate against \\(B\\) | Transferability to \\(D\\)               |

These are **proposed test values**, not literature-recommended material length scales or claimed physical parameters.

All runs need appropriate mesh resolution; the listed ratios alone do not establish it.

A particularly informative outcome would be:

- the recalibrated smaller-length models still reproduce \\(B\\);
- their predictions for \\(D\\) move significantly toward the experimental failure strain;
- the location and width of localization agree better with the observed band.

That would implicate regularization as an important contributor to the original mismatch.

Conversely, if the orientation gap remains nearly unchanged over a converged length-scale study, a constitutive crack-driving modification becomes more compelling.

**My priority would be to perform this study before fitting orientation-dependent toughness or damage parameters.**

## **7. Does one universal parameter explain both orientation classes?**

There are three different levels of universality.

**Level 1: One universal scalar parameter.** A single \\(\zeta\\), \\(W_0\\), \\(p_0\\), or \\(G_c\\) is shared by every specimen.

Such a parameter does not necessarily create a new orientation-discriminating mechanism. In particular, a common monotone function of accumulated plastic work cannot reverse strictly ordered local work histories under otherwise identical conditions.

**Level 2: One universal constitutive function.** For example,

\\[ \chi=\chi(\boldsymbol\rho,\boldsymbol\gamma,\mathbf Q) \\]

could use one parameter set for all crystals while evaluating different values as their microstructures evolve.

This is compatible with material parameter transferability. The important distinction is that the orientation dependence comes from constitutive state variables rather than five separately assigned fracture parameters.

**Level 3: Explicitly orientation-dependent calibrated constants.** For example,

\\[ G_c=G_c(\mathbf Q_0), \\]

with an independently fitted toughness for each initial orientation.

This can improve agreement but has weaker predictive significance unless those values are independently measured or linked to a defensible crystallographic fracture-energy law.

For your benchmark, **Level 2 is the preferred modeling target**.

I cannot confirm whether a specific published CP-PFD model has achieved your requested two-class agreement using one global parameter, or whether published single-crystal studies rely on slip-system-resolved variables. That requires paper-level verification.

## **8. An efficient model-discrimination workflow**

Before adding another constitutive parameter, extract the following quantities from the existing simulations.

**Diagnostic extraction checklist**

9 of 9 selected

Accumulated plastic work, Wp

Accumulated plastic strain, p

Slip-system dislocation densities, ρα

Distribution of slip activity over the 12 systems

Lattice rotation history

Triaxiality and Lode parameter histories

Localization onset and local/global strain

Damage onset, peak load, and final separation

Experimental fracture-plane orientation

&#x20;Copy selected diagnostics

Use these to answer the following decision questions.

| Diagnostic finding                                            | Interpretation                                                 | Next model                                    |
| ------------------------------------------------------------- | -------------------------------------------------------------- | --------------------------------------------- |
| High \\(W_p\\), relatively low stored-energy proxy in \\(D\\) | Plastic-work conversion is plausible                           | Stored-energy driver                          |
| Similar \\(p\\), strongly different \\(W_p\\)                 | Damage is overly sensitive to flow stress                      | Plastic-strain gate                           |
| Different triaxiality or Lode histories before failure        | Local stress state may discriminate failure                    | Stress-state-sensitive damage                 |
| Strongly crystallographic fracture surfaces                   | Directional material separation may matter                     | Anisotropic resistance / slip-resolved damage |
| Localization captured but separation premature                | Fracture constitutive law or regularization likely problematic | Initiation criterion and length-scale study   |
| Localization itself predicted incorrectly                     | Fracture parameter fitting is premature                        | Revisit crystal plasticity and geometry       |
| Large change in orientation error with \\(l_0\\)              | Regularization is materially affecting predictions             | Length-scale calibration                      |

The comparisons must use the same definition of strain and the same event in both specimens. Comparing local strain at final failure in one specimen against nominal strain at damage onset in another would be misleading.

## **9. Direct answers to the six questions**

**1. Which modifications are available conceptually?**

The requested classes can all be represented through changes to the damage-driving energy, damage-initiation condition, fracture resistance, or additional internal damage variables. Section 3 presents their mathematical structures and physical meanings. Their exact published forms and attribution remain unverified.

**2. Which have been applied to crystal plasticity and validated against multi-orientation single-crystal experiments?**

I cannot verify the publication-level answer. Mathematical compatibility with CP is not evidence of published implementation or experimental validation. I have therefore marked the application and validation columns as unverified rather than supplying speculative yes/no answers.

**3. Is thin single-crystal nickel fracture controlled by instability or separation?**

The two mechanisms can operate sequentially. Your large observed local strain makes plastic localization especially important to investigate, but it does not establish the microscopic separation mechanism. For the nickel orientations in question, this requires experimental fracture and deformation evidence.

**4. Which change is most likely to delay high-flow-stress orientations?**

A physically justified stored-energy fraction or slip-state-dependent driving quantity is the most direct way to decouple crack driving from total plastic work.

A common monotonic work threshold cannot reverse a strict ordering in accumulated work. A common plastic-strain gate can help if plastic-strain histories have a different ordering. Decreasing \\(G_c\\) with accumulated plastic strain is unlikely to help if the ductile-like specimens also accumulate greater plastic strain, although differing local histories can change that conclusion.

**5. Is there a universal scalar parameter?**

No universal scalar is established by the information available. A universal state-dependent law is a more promising target than independently calibrated orientation constants, but its success must be demonstrated.

**6. Could the length scale be responsible?**

Yes. The ratio \\(l_0/t=1.71\\) makes this a substantial and testable possibility. It does not prove that regularization causes the early failure. A mesh-converged, recalibrated length-scale study is needed.

## **10. Reference status and claims not verified**

### **Reference list**

**No independently verified bibliographic entries or DOI records can be supplied in this session.**

The benchmark identified in your question is:

Rezwan et al., *“Predicting ...”*, *Computational Materials Science*, volume 272 (2026), article 114872.

That bibliographic identification is user-supplied and remains unverified here. I cannot provide a trustworthy DOI for it without checking the source record.

I have intentionally omitted author/year/DOI combinations recalled from memory because your instruction is to cite only opened and verified sources.

### **Claims requiring verification**

1. The exact fracture formulations, parameter definitions, and thermodynamic treatments in the published isotropic and crystal-plasticity ductile phase-field literature.
2. Which published CP-PFD models have been validated against tensile failure of multiple orientations of the same single-crystal material.
3. The measured ordering of tensile ductility and fracture mechanisms across FCC nickel single-crystal orientations, including the role of single versus multiple slip.
4. Whether the Rezwan et al. benchmark independently identifies regularization-length effects, stored plastic energy, or orientation-dependent fracture resistance.
5. Whether any verified study demonstrates prediction of both your ductile-like and brittle-like orientation groups with one state-dependent fracture parameterization.

These gaps prevent this response from being considered a completed, source-verified literature review.

## **Final modeling conclusion**

The most important result from your equations is that **total accumulated plastic work can be a poor discriminator of orientation-dependent ductility when flow stress and resistance to material separation do not follow the same ordering**.

Your hypothesis is therefore mechanically plausible, but it is not uniquely supported by the current failure-strain mismatch. A stored-energy or slip-state-sensitive driver could resolve it; a threshold may also help; and the unusually large regularization length could independently bias the structural failure response.

The strongest minimal investigation is to compare one specimen from each class, extract their undamaged CP histories, determine whether plastic work and stored-energy proxies actually separate the groups, and then repeat the two-specimen fracture calculation at smaller, mesh-converged phase-field lengths.

That would establish whether the missing physics is primarily **energy storage, damage initiation, crystallographic separation, or the regularization of a plastic instability**—before introducing parameters that risk fitting the observed failure strains without explaining them.