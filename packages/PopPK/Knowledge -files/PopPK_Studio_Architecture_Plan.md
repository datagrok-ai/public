# Building a browser-based PopPK/PD engine on Datagrok

**Architecture, algorithms, and roadmap**

---

## Datagrok provides a surprisingly strong foundation

Datagrok's architecture is well-suited to this task. The platform core is written in Dart (shared between client and server), with a JavaScript/TypeScript API (`datagrok-api`) as the primary developer interface. Its proprietary in-memory columnar database can handle tens of millions of rows in the browser — far exceeding typical pharmacometric dataset sizes. The plugin system supports packages written in TypeScript that can register custom viewers, applications, and computational functions.

Three existing Datagrok components are directly relevant. First, **diff-grok** is a purpose-built TypeScript ODE solver supporting stiff equations via Rosenbrock methods (MRT, ROS3PRw, ROS34PRw), with `Float64Array` performance and WebWorker support. This is critical because PK/PD models frequently exhibit stiffness. Second, the **compute-utils** library provides FittingView with in-WebWorker parameter optimization (currently using Nelder-Mead), goodness-of-fit viewers, and outlier selection. Third, the **PipelineView/TreeWizard** framework enables hierarchical multi-step computational workflows with dependency tracking — exactly the pattern needed for a stepwise covariate modeling workflow.

The platform also supports server-side compute via R, Python, and Julia through a message queue architecture. This means nlmixr2 (R) or Pumas (Julia) could serve as fallback estimation engines while the browser-native engine is built, and for computationally intensive tasks like bootstrap or full Bayesian estimation that may exceed browser capacity.

What Datagrok lacks entirely: NLME estimation algorithms, population-level random effects modeling, NMTRAN data handling, VPC generation, bootstrap infrastructure, and a pharmacometric model specification language.

---

## The NLME estimation engine requires two core algorithms

The central challenge in population pharmacometrics is computing the marginal likelihood by integrating out individual random effects (η) from the joint distribution. For N subjects, the population likelihood is **L(θ, Ω, Σ) = ∏ᵢ ∫ p(yᵢ | θ, ηᵢ, Σ) · p(ηᵢ | Ω) dηᵢ**. This integral has no closed-form solution for nonlinear models. Every estimation method is an approximation strategy.

**FOCE-I** (the industry workhorse) uses a Laplacian approximation with first-order linearization. It operates as a nested optimization: the **inner loop** finds each individual's MAP estimate η̂ᵢ by minimizing the individual negative log-posterior via Newton-Raphson or L-BFGS, while the **outer loop** minimizes the approximate population -2LL over all population parameters (θ, Ω, Σ) using quasi-Newton methods. The "interaction" means the residual variance matrix R is evaluated at η̂ᵢ rather than η=0, capturing ε-η interaction — essential for proportional and combined error models. The individual contribution to the objective function is: **OFVᵢ = η̂ᵢᵀΩ⁻¹η̂ᵢ + ln|Cᵢ| + εᵢᵀCᵢ⁻¹εᵢ + nᵢ·ln(2π)**, where Cᵢ is the linearized covariance matrix incorporating both between-subject and residual variability. The computational bottleneck is that inner optimization runs for every subject at every outer iteration, with each evaluation requiring an ODE solve.

**SAEM** (Stochastic Approximation EM) takes a fundamentally different approach. Instead of optimizing, it treats individual parameters as missing data and iterates between an E-step (MCMC sampling of individual parameters via Metropolis-Hastings) and an M-step (closed-form updates of population parameters from sufficient statistics). Monolix's implementation uses two phases: an **exploratory phase** with step-size γₖ = 1 and simulated annealing for global exploration, followed by a **smoothing phase** with γₖ = 1/(k−K₁) guaranteeing almost-sure convergence to the MLE. SAEM requires no gradient computation and no nested optimization, making it structurally simpler to implement. Each iteration needs N × L × m model evaluations (L chains per subject, m Metropolis-Hastings transitions), but the M-step is just matrix algebra: μₖ = (1/N)·sₖ,₁ and Ωₖ = (1/N)·(sₖ,₂ − sₖ,₁·sₖ,₁ᵀ).

For a Datagrok implementation, **build SAEM first, then FOCE-I**. SAEM is structurally simpler (no gradients, no nested optimization), more robust to starting values (thanks to simulated annealing), and mathematically guaranteed to converge. FOCE-I is needed for NONMEM compatibility and for computing the standard -2LL objective function (SAEM does not directly produce an OFV — it requires a separate likelihood computation step via importance sampling or linearization).

**MAP/EBE estimation** should be implemented before either algorithm. It is the inner problem of FOCE-I and is independently valuable for model-informed precision dosing (MIPD). For individual i, it minimizes φᵢ(ηᵢ) = Σⱼ[(yᵢⱼ − f(θ, ηᵢ, tᵢⱼ))²/σ²ᵢⱼ + ln(σ²ᵢⱼ)] + ηᵢᵀΩ⁻¹ηᵢ. This is a low-dimensional optimization (2–5 parameters) solvable with L-BFGS or even Nelder-Mead.

Full Bayesian estimation via HMC/NUTS is **not recommended for browser implementation**. It requires mature automatic differentiation through the ODE solver (nothing comparable to Stan Math exists in JS), and computation is 10–100× more expensive than SAEM. Delegate this to server-side Stan/Torsten or Pumas via Datagrok's compute queue.

---

## What the JS/TS numerical ecosystem offers and where it falls short

The JavaScript ecosystem has matured significantly for scientific computing, though critical gaps remain.

**Available and production-ready.** The **stdlib-js** project provides the most comprehensive scientific computing foundation: full BLAS Level 1–3 implementations with WASM fallbacks, 35+ probability distributions, and 40+ seedable PRNGs. For dense linear algebra, **ml-matrix** (~84 Snyk health score) offers Cholesky decomposition, LU, QR, SVD, eigenvalue decomposition, matrix inversion, and determinants — all operations needed for the NLME covariance step. The **ml-levenberg-marquardt** package (~17,000 weekly downloads) provides a mature Levenberg-Marquardt implementation. For L-BFGS, **optimization-js** offers a working implementation, though a custom implementation may be needed for bounded optimization (L-BFGS-B). Statistical distributions from **jStat** or stdlib-js cover all pharmacometric needs (multivariate normal, chi-squared, etc.).

**Available but needs WASM compilation.** The most critical performance dependency is ODE solving. While diff-grok provides Rosenbrock methods, compiling **SUNDIALS CVODE/CVODES** to WebAssembly via Emscripten would deliver BDF methods for stiff systems plus — critically — **forward sensitivity analysis** (CVODES), which enables gradient computation without finite differences. SUNDIALS is pure C with CMake and minimal dependencies; the serial N_Vector implementation ports cleanly to WASM. **Eigen** compiles to WASM successfully at ~40–60% of native speed, providing high-performance dense linear algebra. **NLopt** (C, self-contained) should compile straightforwardly to WASM, providing access to 20+ optimization algorithms including L-BFGS, BOBYQA, and COBYLA.

**Missing and must be built.** There is no production-quality automatic differentiation library for JS/TS suitable for pharmacometrics. Forward sensitivity equations (augmenting the ODE system to compute ∂y/∂η alongside the state) are the practical alternative — they can be implemented by symbolically differentiating the model equations or using finite differences within the ODE RHS. No JS/TS MCMC library provides HMC or NUTS; a basic Metropolis-Hastings sampler (needed for SAEM) must be written from scratch, though the algorithm is straightforward.

**Performance is viable.** WASM achieves **1.45–1.55× native speed** on average, and with SIMD enabled (`-msimd128`), numerical kernels run within 1.2–1.8× of native. For a two-compartment oral PK model with 100 subjects, estimated FOCE-I wall time in the browser is 5–30 seconds using WASM ODE solver + 8 WebWorkers for per-subject parallelism. This is slower than native C++ (1–5 seconds) but entirely usable for interactive workflows.

---

## Architecture blueprint for the computational core

The recommended architecture has three tiers, matching the pattern used successfully by nlmixr2 (R frontend → C++ core) and Pumas (Julia all-the-way-down), adapted for the browser.

**Tier 1: WASM numerical engine.** Compile SUNDIALS CVODES + Eigen + NLopt to a single WASM module exposed via a clean TypeScript interface. This module handles: ODE integration (BDF/Adams), forward sensitivity equations, dense matrix factorization (Cholesky, eigendecomposition), and constrained optimization (L-BFGS-B, BOBYQA). All performance-critical inner-loop operations run here. Alternatively, if compilation complexity is too high initially, use diff-grok (Rosenbrock) for ODE solving and ml-matrix for linear algebra in pure TS — this is slower but functional and allows faster initial development.

**Tier 2: TypeScript estimation layer.** The NLME algorithms (SAEM, FOCE-I, MAP) are implemented in TypeScript, calling Tier 1 for numerical primitives. This layer manages the estimation loop, convergence monitoring, parameter transformations, and result aggregation. Key components:

- **Model compiler**: Parses the model DSL, generates optimized evaluation functions (or WASM modules) for the structural model, its Jacobian, and derived quantities
- **Data manager**: Handles NMTRAN-format datasets — ADDL/II expansion, sorting, BLQ handling, covariate interpolation, and splitting data by individual
- **FOCE-I engine**: Outer loop (quasi-Newton via L-BFGS-B), inner loop (MAP estimation per subject), objective function assembly, covariance step (R-matrix via finite differences, S-matrix from individual gradients, sandwich estimator)
- **SAEM engine**: Metropolis-Hastings sampler with multiple transition kernels, sufficient statistics accumulator, simulated annealing schedule, convergence monitoring
- **Post-estimation**: CWRES computation, shrinkage calculation, NPDE via simulation, individual fits

**Tier 3: WebWorker parallelism.** The per-subject computations (ODE solving, MAP estimation, MCMC sampling) are embarrassingly parallel. Deploy a WebWorker pool sized to `navigator.hardwareConcurrency`. Each worker loads the WASM module and processes a batch of subjects. SharedArrayBuffer enables zero-copy sharing of the population parameters and dataset. The main thread orchestrates the estimation loop and updates the UI with convergence diagnostics in real time.

For tasks exceeding browser capacity — bootstrap with 1000 replicates, full Bayesian estimation, very large datasets (50,000+ observations) — route computation to Datagrok's server-side compute queue running nlmixr2 (R) or a containerized native NLME engine. The same model specification should compile to both browser and server targets.

---

## How existing platforms inform model specification and data design

Learning from NONMEM, Monolix, nlmixr2, and Pumas reveals convergent patterns that should guide Datagrok's design.

**Model specification** across all platforms separates into five blocks: parameter definitions, random effects structure, pre-computation (covariate-parameter relationships), dynamics (ODEs or analytical solutions), and observation model (error structure). The nlmixr2 DSL provides the best template for a browser implementation: an `ini` block for parameter initialization (θ values with bounds, Ω matrix via `eta ~ variance` syntax, Σ parameters) and a `model` block for the structural model with `d/dt()` ODE notation and `DV ~ add(σ) + prop(σ)` error specification. Pumas adds the useful pattern of explicit `@param` domains with constraints (`RealDomain(lower=0)`) and the `:=` operator for intermediate variables that won't appear in output.

For Datagrok, implement a **dual-interface model specification**: a text-based DSL (targeting the pharmacometrician persona) and a graphical model builder (targeting the clinical pharmacologist). Both should compile to the same internal JSON-like model representation:

```json
{
  "structural": { "type": "ode|analytical", "compartments": [...], "equations": {...}, "output": "..." },
  "parameters": { "theta": [{"name": "", "init": 0, "lower": 0, "upper": 0, "transform": ""}], "omega": {"structure": "", "values": []}, "sigma": {...} },
  "covariates": [{"parameter": "", "covariate": "", "type": "", "reference": 0}],
  "estimation": { "method": "focei|saem", "settings": {...} }
}
```

**The PREDPP-equivalent model library** is essential for performance. Standard analytical solutions (1-cmt IV, 1-cmt oral, 2-cmt IV, 2-cmt oral, 3-cmt variants) avoid ODE solving entirely, providing 100–1000× speedup. Implement these as optimized TypeScript functions with the standard TRANS parameterizations: microconstants (K, K12, K21), CL/V (TRANS2), CL/V1/Q/V2 (TRANS4). For a two-compartment oral model, the analytical solution is a triexponential involving eigenvalues of the transfer rate matrix — precompute these from CL, V1, Q, V2, and KA. Include superposition for multiple-dose scenarios and steady-state computation.

**NMTRAN data handling** is non-negotiable for interoperability. The data manager must parse the standard columns (ID, TIME, DV, AMT, EVID, CMT, MDV, RATE, ADDL, II, SS), handle EVID codes (0=observation, 1=dose, 2=other, 3=reset, 4=reset+dose), expand ADDL/II into explicit dose records, enforce sorting (by ID then TIME, doses before observations at same time), and carry forward dose-related covariates. Datagrok's powerful data grid can render these datasets with domain-aware formatting — color-coding dosing vs. observation records, flagging BQL observations, showing time-after-dose alongside clock time.

---

## Diagnostics, VPC, and bootstrap complete the pipeline

**Standard diagnostic plots** should be generated automatically post-estimation. The core set requires: DV vs PRED (population predictions), DV vs IPRED (individual predictions), CWRES vs TIME, CWRES vs PRED, |IWRES| vs IPRED, QQ plot of CWRES, ETAs vs covariates, and individual fits. CWRES computation uses the FOCE linearization: CWRESᵢ = Cᵢ^(-1/2) · (yᵢ − E[yᵢ]), where Cᵢ = Gᵢ·Ω·Gᵢᵀ + Rᵢ. A critical companion metric is **shrinkage**: η-shrinkage = 1 − SD(η̂ₖ)/ωₖ and ε-shrinkage = 1 − SD(IWRES). When shrinkage exceeds 20–30%, EBE-based diagnostics become unreliable — the platform should warn users prominently. Datagrok's synchronized viewers (scatter plots, histograms, box plots sharing the same data frame) are ideal for this — clicking a subject in one plot highlights them everywhere.

**NPDE (Normalized Prediction Distribution Errors)** provides simulation-based, model-agnostic diagnostics. The algorithm simulates 300–1000 datasets from the model, computes the empirical CDF of simulated values for each observation, then transforms to standard normal via Φ⁻¹. If the model is correct, NPDEs follow N(0,1). This is computationally intensive (hundreds of simulation runs) but parallelizes perfectly across WebWorkers.

**Visual Predictive Checks** require simulating 500–1000 replicate datasets using the final estimates, binning the independent variable (TIME or TAD), computing percentiles (5th, 50th, 95th) within each bin for both observed and simulated data, then overlaying observed percentile lines against simulated prediction interval bands. **Prediction-corrected VPCs** (pcVPC) normalize data within bins by the ratio PRED_bin_median / PREDᵢⱼ, removing covariate-driven variability. Binning strategies range from equal-width intervals to K-means clustering on the independent variable. With Datagrok's custom viewer API, the VPC plot can be implemented as a specialized viewer with interactive bin adjustment.

**Bootstrap resampling** is embarrassingly parallel: resample individuals with replacement, re-estimate the model, collect parameters. With 1000 replicates and a 10-second per-replicate browser-based estimation, total time would be ~2.7 hours serially — impractical in-browser. Two solutions: (1) run bootstrap replicates in parallel across server-side compute nodes via Datagrok's queue architecture, or (2) use **case-deletion diagnostics** and the **SIR (Sampling Importance Resampling)** method as faster alternatives that require only one model estimation plus simulation.

---

## Serving two personas with one architecture

The pharmacometrician and clinical pharmacologist need fundamentally different interfaces to the same engine.

**For the pharmacometrician**, the platform must provide a code editor with syntax highlighting for the model DSL, a console showing estimation progress (OFV by iteration, gradient norms, convergence status), direct access to the objective function value and parameter estimates, scripting capability for custom workflows (stepwise covariate modeling, sensitivity analyses), and the ability to export/import NONMEM-compatible control streams via a translation layer (analogous to nlmixr2's babelmixr2). The Diff Studio package already demonstrates this pattern — a code editor for model specification with real-time parameter exploration. Extend it with NONMEM-like syntax support and estimation method selection. Critical: parameter estimates must match NONMEM to 3+ significant figures for regulatory acceptance.

**For the clinical pharmacologist**, the platform should offer a wizard-driven workflow: (1) upload/connect data → auto-detect NMTRAN columns, (2) select a pre-built model from the library (with visual compartmental diagrams), (3) click "Estimate" with sensible defaults, (4) review auto-generated diagnostic report, (5) point-and-click covariate selection with forest plots showing effect sizes, (6) one-click VPC and bootstrap, (7) automated Word/PDF reporting. Datagrok's RichFunctionView (auto-generated UI from function annotations) and PipelineView (multi-step workflow) are purpose-built for this pattern. The model library should include the **30+ most common pharmacometric models**: 1/2/3-compartment PK with various absorption models (first-order, zero-order, transit compartment, Weibull), Michaelis-Menten elimination, target-mediated drug disposition (TMDD), Emax PD models, indirect response models, and tumor growth inhibition models.

The two interfaces should share the same underlying model representation and estimation engine. A pharmacometrician's hand-coded model should appear identically in the clinical pharmacologist's GUI, and vice versa. Monolix demonstrates this pattern well: every GUI action produces an editable Mlxtran text file.

---

## A phased roadmap from MVP to regulatory-grade platform

**Phase 1 (Months 1–4): MAP estimation + analytical models.** Build the NMTRAN data parser, implement analytical solutions for 1-cmt and 2-cmt models (with TRANS2/TRANS4 parameterizations), implement MAP/EBE estimation using Nelder-Mead (already available) and L-BFGS (from optimization-js or custom), create basic individual fit plots. This delivers immediate value for **model-informed precision dosing** (MIPD/TDM) — a high-demand clinical application. Use diff-grok for ODE-based models. Deliverables: data import, model library (6 analytical models), individual parameter estimation, individual fit visualization.

**Phase 2 (Months 4–8): SAEM estimation engine.** Implement Metropolis-Hastings sampler, sufficient statistics accumulator, simulated annealing schedule, and the two-phase SAEM algorithm. Add population-level parameter estimation, between-subject variability, and residual error model estimation. Implement post-hoc likelihood computation (via linearization). Add standard diagnostics (GOF plots, CWRES, shrinkage). This is the first true **population PK** capability. WebWorker parallelism for per-subject MCMC sampling. Deliverables: SAEM estimation, population parameter estimates with standard errors, core diagnostic suite.

**Phase 3 (Months 8–12): FOCE-I and model qualification.** Implement the FOCE-I objective function (requires forward sensitivity equations or finite-difference gradients), nested optimization architecture, and the covariance step (R-matrix, S-matrix, sandwich estimator). Add VPC simulation engine, bootstrap infrastructure (server-side), stepwise covariate modeling, and the model DSL with code editor. Compile SUNDIALS CVODES to WASM for performance-critical stiff ODE solving with sensitivity analysis. Deliverables: FOCE-I estimation matching NONMEM OFV, VPC viewer, bootstrap, covariate modeling workflow.

**Phase 4 (Months 12–18): Full platform maturity.** NONMEM control stream import/export (via translation layer), NPDE diagnostics, FREM covariate analysis, automated reporting (Word/PDF), model comparison tools (AIC/BIC tables), simulation workspace with dosing optimization, server-side Bayesian estimation via Stan/Torsten, and regulatory-grade validation documentation. Deliverables: complete pharmacometrics workbench competitive with Monolix's GUI capabilities while offering NONMEM-level computational rigor.

**Performance strategy**: for datasets under 500 subjects with analytical models, run everything in the browser. For ODE models or >500 subjects, offer a "compute on server" toggle that routes estimation to Datagrok's compute nodes running the same algorithms compiled natively (or nlmixr2/Pumas). The user experience should be identical regardless of where computation runs — only wall time changes.

---

## Key technical decisions that will determine success

**ODE solver strategy**: Start with diff-grok (Rosenbrock, already integrated) for immediate functionality, but invest in WASM-compiling SUNDIALS CVODES for Phase 3. The forward sensitivity analysis capability of CVODES eliminates the need for finite-difference gradient computation in FOCE-I, providing both better accuracy and **10–100× speedup** in gradient evaluation. Without sensitivity equations, FOCE-I requires ~2p additional full model evaluations per outer iteration (p = number of parameters), each requiring all N inner optimizations.

**Gradient computation**: Forward sensitivity equations (augmented ODE system computing ∂y/∂η alongside the state) are strongly preferred over finite differences. They can be implemented by symbolically differentiating the model's ODE right-hand side — nlmixr2/rxode2 does this via symengine, generating C code for the augmented system. In TypeScript, a lightweight symbolic differentiation pass during model compilation can generate the sensitivity ODEs. For analytical models, gradients are available in closed form.

**Model compilation**: Follow rxode2's successful pattern — transpile the model DSL into optimized evaluation functions. In the JS/TS context, this means: parse the model → generate a JavaScript function (or WASM module) that computes model predictions, Jacobians, and derived quantities for a given set of parameters and time points. JIT compilation by V8 will optimize these functions after a few calls. For even higher performance, generate C code and compile to WASM via the Emscripten toolchain, analogous to how Pumas leverages Julia's JIT.

**Numerical validation**: The FOCE-I implementation must reproduce NONMEM's OFV to at least 3 significant figures on benchmark datasets (theophylline, warfarin, phenobarbital) to be credible in the pharmacometrics community. nlmixr2 achieved this and published validation results — follow their test suite. The Wang (2007) paper provides the definitive reference for objective function derivations.
