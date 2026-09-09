# Project guidelines — PINNACLE CP–PFD benchmark

Common rules for **every conversation in this project**. Paste the "Project
description" block at the end of this file into the project's description field
so each new conversation inherits them.

## 1. Language and style

American English throughout — code, comments, decks, documentation, figures,
and prose. Lab-notebook register: state the result, the artifact, the caveat.

## 2. Citations

Cite in prose as **(author year significant_word)** — first author surname, year,
one distinguishing word: `(Rezwan 2025 intragranular)`, `(Lim 2019 meshsensitivity)`,
`(Abedi 2026 internallength)`.

The **single central bibliography** is `bib/references.bib` in the study
repository, mirrored as a project artifact named `references.bib`. Its keys are
`author-year-significantword`, lowercase. Do not create a second bibliography, do
not invent a key, and do not cite a work that is absent from it — add the
verified entry first. `bib/citation_keys.md` maps every key to its full citation
and prose form.

No fabricated references. Every entry carries a DOI where one exists and was
verified against Crossref; unverified entries are flagged as such in
`citation_keys.md` and must stay flagged until verified.

## 3. File placement — study files never enter the framework

| Content | Location |
|---|---|
| Input decks, meshes, datasets, figures, results, reports, scripts | `/mnt/extended/projects/pinn-cp-pfd-benchmark/` |
| New framework capability (components, tests, ADRs, docs) | `/mnt/extended/projects/pinnacle/`, branch **`feat/pinn_benchmark`** |
| Mesh generation tooling | `/mnt/extended/projects/hexaform/` |

Nothing specific to this study — no orientation list, no calibrated parameter
set, no specimen geometry, no figure — is committed to the PINNACLE tree. A new
PINNACLE component is justified only if it is general: a registered, documented,
tested capability that any user could compose. Study conventions live in the
study repository.

All framework work goes on `feat/pinn_benchmark`, never on `master`.

## 4. Evidence and provenance

Literature verification precedes computation at every stage. Before a value
enters a production deck, its source, definition, units, coordinate system, and
uncertainty are identified and recorded.

The companion data sheet (`PINNACLE Benchmark — Extracted Source Data`) is a
navigation aid, **not** an authoritative source. Values are extracted
independently from the primary publications, then reconciled, with every
discrepancy documented.

Every quantity is tagged with exactly one provenance: `reported_table`,
`reported_text`, `digitized`, `derived`, `calibrated`, `assumed`, `unavailable`.
No parameter of unknown provenance goes into a production deck. Assumptions
required to close a gap in the published record are labeled and evaluated by
sensitivity analysis.

Digitized values carry the digitization method, image resolution, axis
calibration, extraction uncertainty, and repeated-extraction error. Experimental
curves are never rescaled, shifted, smoothed, or truncated without a documented
physical and statistical justification.

## 5. Accuracy classification

Every reported agreement states which tier it is:

1. regression or golden-reference agreement;
2. independently verified numerical agreement;
3. agreement with a previously published simulation;
4. agreement with experimental measurements.

Tier 1 demonstrates implementation consistency, not physical validity. A neural
route is never called validated on the strength of agreeing with FEM: if FEM and
the neural route agree with each other but both disagree with experiment, that
indicates shared model-form or parameter error.

## 6. Compute

Use the **connected SSH compute target** (`ssh:rocm-server`). Do **not** use
cloud compute. Prefer PyTorch on the GPU at `float64`, falling back to CPU only
when the GPU is unavailable — with one measured exception recorded below.

**Measured exception.** The crystal-plasticity field solve is hundreds of
thousands of small 3×3 einsums where per-op dispatch dominates; PINNACLE's own
bundled bicrystal deck runs in 17 s on NumPy and does not finish in 50 minutes
on Torch (ADR 0015), and its batch constitutive path marshals a Python list of
per-quadrature-point dicts on every call. Backend selection is therefore
**per workload, decided by measurement** — see `results/backend_envelope.*` —
not global. Training, collocation-batch evaluation, autodiff, and dense algebra
go to the GPU at `float64`; the CP field solve uses whichever backend the
envelope measurement shows is faster at that problem size. Report algorithmic
restructuring separately from hardware effect.

Host-specific requirements: export `TORCH_BLAS_PREFER_HIPBLASLT=0` (fp32 kernel
selection and GPU-hang fix) and `TORCH_ROCM_AOTRITON_ENABLE_EXPERIMENTAL=1`
(attention); bound every GPU loop; run each measurement in its own subprocess;
check `uptime` before timing. Never use `autocast`/`GradScaler` for PINN work —
they target fp16/bf16 GEMMs and interact badly with `create_graph=True`.

## 7. Reproducibility

Every calculation records: code revision, deck revision, dependency environment,
hardware, numerical precision, random seed, mesh, quadrature rule, load
increment, solver tolerances, stopping criteria, and initialization procedure.

Neural results are reported over multiple independent seeds as mean, standard
deviation, and worst case — never a single run.

Raw data, extraction scripts, meshes, decks, trained weights, logs, and
post-processing procedures are retained as versioned artifacts.

Cost accounting separates preprocessing, training, inference, and correction,
and includes failed training runs, hyperparameter selection, and model
development in a transparent, separately reported total.

## 8. Framework changes

Evaluate existing PINNACLE capability before writing anything new: `docs/user_guide.md`
§1.1 (physics × solver matrix) and §1.2 (capability index), and
`pinnacle syntax dump`. Implement a benchmark first by composing registered
components in a new `.i` deck.

If a capability is genuinely missing, add it per `CLAUDE.md` §3 and
`docs/developing_components.md`: register with `@registry.register(...)`, supply
weak and strong forms where applicable, document, test, add to the parity
harness, and record the architectural decision in the ADR series. New parameters
default to values that preserve existing goldens.

Do not present an unsupported solver route as a completed benchmark. A pure-PINN
coupled CP–PFD field solve is structurally unavailable (ADR 0017) and is reported
as such.

---

## Project description block

> **PINNACLE CP–PFD benchmark.** Comparative FE / PINN+FEM / neural-constitutive
> study of coupled crystal plasticity and phase-field damage, validated against
> Rezwan 2025 (Al2219 single crystal with an intracrystalline void) and Rezwan
> 2026 (single-crystal Ni tensile specimens).
>
> All conversations follow `docs/PROJECT_GUIDELINES.md` in
> `/mnt/extended/projects/pinn-cp-pfd-benchmark`. In brief:
>
> - **American English** everywhere.
> - **Citations** in prose as `(author year significant_word)`; the single
>   central bibliography is `bib/references.bib` (keys
>   `author-year-significantword`), mirrored as the project artifact
>   `references.bib`. Never cite a work absent from it; never fabricate a
>   reference or DOI.
> - **Study files stay out of the framework.** Decks, meshes, datasets, figures,
>   and reports live in the study repository; new PINNACLE capability goes to
>   `/mnt/extended/projects/pinnacle` on branch **`feat/pinn_benchmark`** only,
>   and only if it is general, registered, documented, and tested. Meshing
>   tooling lives in `/mnt/extended/projects/hexaform`.
> - **Provenance before computation.** Extract from primary sources, not from the
>   companion data sheet; tag every value `reported_table`, `reported_text`,
>   `digitized`, `derived`, `calibrated`, `assumed`, or `unavailable`; label
>   assumptions and test them by sensitivity analysis. Never rescale, shift,
>   smooth, or truncate an experimental curve.
> - **State the accuracy tier** of every claimed agreement: golden-reference,
>   independently verified numerical, previously published simulation, or
>   experiment. Agreement with FEM alone never validates a neural route.
> - **Compute on the connected SSH target, not cloud compute.** Prefer Torch on
>   GPU at `float64`; fall back to CPU when unavailable. Backend choice for the
>   CP field solve is decided by measurement, because that workload is
>   dispatch-bound (ADR 0015). Export `TORCH_BLAS_PREFER_HIPBLASLT=0` and
>   `TORCH_ROCM_AOTRITON_ENABLE_EXPERIMENTAL=1`; bound GPU loops.
> - **Reproducibility is mandatory:** record revision, environment, hardware,
>   precision, seed, mesh, quadrature, load increment, tolerances, and
>   initialization for every run; report neural results over multiple seeds as
>   mean, standard deviation, and worst case; account for training and
>   development cost separately.
