# PopPK Studio — Algorithms and Numerical Methods

## Datagrok current capabilities (DiffStudio / diff-grok)

Based on the [DiffStudio documentation](https://datagrok.ai/help/compute/diff-studio) and the [diff-grok npm package](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio), Datagrok currently provides:

**ODE solvers:**
- Implicit / stiff: ROS34PRw (default), ROS3PRw, MRT — Rosenbrock-Wanner methods
- Explicit / non-stiff: Bogacki-Shampine RK3(2), Runge-Kutta-Fehlberg RK4(5), Dormand-Prince RK5(4), Adams multistep AB4, AB5
- Adaptive: **LSODA** (automatic Adams ↔ BDF switching), **CVODE** (SUNDIALS variable-order Adams + BDF)

**Model optimization & analysis:**
- Parameter optimization (fitting to experimental data) via Nelder-Mead in WebWorkers
- Sensitivity analysis: Monte Carlo, Sobol variance decomposition, Grid
- Multi-stage simulation (#loop, #update blocks)
- Model DSL (IVP format) with auto-generated UI
- Model Hub integration and script export

---

## 1. Minimal set

These algorithms are required to deliver MAP estimation, SAEM-based population estimation, and core diagnostics.

### 1.1 ODE Solvers

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Rosenbrock-Wanner methods (ROS34PRw, ROS3PRw, MRT) | Stiff ODE integration for compartmental PK/PD models; A-stable implicit methods | **Yes** — diff-grok | Rang & Angermann (2005) "New Rosenbrock W-Methods of Order 3 for PDE"; Hairer & Wanner (1996) "Solving ODEs II" | **JS/TS**: diff-grok; **C/C++**: SUNDIALS, Boost.Odeint; **Fortran**: LSODA; **Julia**: OrdinaryDiffEq.jl; **R**: deSolve, rxode2; **Python**: scipy.integrate.solve_ivp |
| Runge-Kutta methods (RK3, RK4, Dormand-Prince) | Non-stiff ODE integration; adaptive step-size control | **Yes** — diff-grok (rk3, rk4, rkdp) | Dormand & Prince (1980); Bogacki & Shampine (1989); Hairer, Nørsett & Wanner (1993) "Solving ODEs I" | **JS/TS**: diff-grok; **C/C++**: SUNDIALS CVODE (Adams), Boost.Odeint; **Fortran**: DVODE; **Julia**: Tsit5(); **R**: deSolve::ode; **Python**: scipy RK45 |
| Adams multistep methods (AB4, AB5) | Non-stiff predictor-corrector methods; efficient for smooth problems | **Yes** — diff-grok (ab4, ab5) | Bashforth & Adams (1883); Hairer, Nørsett & Wanner (1993) | **JS/TS**: diff-grok; **C/C++**: SUNDIALS CVODE; **Fortran**: DVODE; **Julia**: VCABM(); **R**: deSolve; **Python**: scipy |
| LSODA (automatic stiff/non-stiff switching) | Adaptive solver switching between Adams (non-stiff) and BDF (stiff) automatically; the classic pharmacometric solver | **Yes** — diff-grok (lsoda) | Hindmarsh (1983) "ODEPACK"; Petzold (1983) "Automatic selection of methods for solving stiff and nonstiff systems" | **JS/TS**: diff-grok; **C**: ODEPACK; **Fortran**: LSODA (original); **Julia**: lsoda() via LSODA.jl; **R**: deSolve::lsoda; **Python**: scipy.integrate.odeint |
| CVODE (SUNDIALS variable-order Adams + BDF) | Production-grade stiff/non-stiff solver; used by NONMEM (ADVAN14) and Monolix; variable-order BDF for stiff systems | **Yes** — diff-grok (cvode) | Hindmarsh et al. (2005) "SUNDIALS: Suite of Nonlinear and Differential/Algebraic Equation Solvers"; Cohen & Hindmarsh (1996) | **JS/TS**: diff-grok; **C**: SUNDIALS CVODE/CVODES (Lawrence Livermore, BSD); **Julia**: CVODE_BDF() via Sundials.jl; **R**: rxode2; **Python**: Assimulo |

### 1.2 Analytical PK Solutions

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Matrix exponential for linear compartmental systems | Analytical solution for 1/2/3-compartment models avoiding ODE integration; 100–1000× faster than numerical ODE | **No** | Teorell (1937); Wagner (1993) "Pharmacokinetics for the Pharmaceutical Scientist"; D'Argenio et al. (2009) | **JS/TS**: ml-matrix (eigendecomposition); **C++**: Eigen (MatrixExp); **Fortran**: NONMEM ADVAN1-4, ADVAN11-12; **Julia**: LinearAlgebra exp(A*t); **R**: rxode2 linCmt(), nlmixr2; **Python**: scipy.linalg.expm |
| Superposition principle for multiple dosing | Compute steady-state and multi-dose profiles from single-dose analytical solutions | **No** | Gibaldi & Perrier (1982) "Pharmacokinetics, 2nd ed."; Rowland & Tozer (2011) "Clinical PK and PD" | **Fortran**: NONMEM SS system; **Julia**: Pumas DosageRegimen; **R**: rxode2; **C++**: Custom |

### 1.3 Optimization Algorithms

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Nelder-Mead simplex | Derivative-free optimization; MAP estimation inner loop; robust but slow for high dimensions | **Yes** — compute-utils, DiffStudio fitting | Nelder & Mead (1965); Lagarias et al. (1998) "Convergence properties" | **JS/TS**: compute-utils (Datagrok); **C**: NLopt; **Fortran**: NONMEM; **Julia**: Optim.jl; **R**: optim(); **Python**: scipy.optimize.minimize |
| L-BFGS (Limited-memory BFGS) | Quasi-Newton gradient-based optimization; outer loop of FOCE-I; fast convergence for smooth objectives | **No** | Nocedal (1980); Liu & Nocedal (1989) "On the limited memory BFGS method" | **JS/TS**: optimization-js, lbfgs-js; **C**: NLopt, L-BFGS-B (Zhu et al., 1997); **Fortran**: L-BFGS-B.f; **Julia**: Optim.jl; **R**: optim("L-BFGS-B"); **Python**: scipy L-BFGS-B |
| L-BFGS-B (box-constrained) | L-BFGS with bound constraints; essential for constrained PK parameters (CL > 0, V > 0) | **No** | Byrd, Lu, Nocedal & Zhu (1995); Zhu et al. (1997); Morales & Nocedal (2011) | **C**: L-BFGS-B 3.0, NLopt; **Fortran**: L-BFGS-B.f; **Julia**: Optim.jl; **R**: optim("L-BFGS-B"); **Python**: scipy L-BFGS-B |
| Levenberg-Marquardt | Nonlinear least-squares; individual curve fitting; hybrid Gauss-Newton / gradient descent | **No** | Levenberg (1944); Marquardt (1963); Moré (1978) | **JS/TS**: ml-levenberg-marquardt; **C**: cminpack, GSL; **Fortran**: MINPACK; **Julia**: LsqFit.jl; **R**: minpack.lm; **Python**: scipy.optimize.least_squares |

### 1.4 NLME Estimation Methods

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| MAP / EBE (Maximum A Posteriori) | Individual parameter estimation given population priors; inner problem of FOCE-I; standalone value for precision dosing (TDM) | **No** | Sheiner & Beal (1981); Beal & Sheiner (1982) "Estimating population kinetics" | **Fortran**: NONMEM POSTHOC; **C++**: Phoenix NLME; **Julia**: Pumas fit(MAP()); **R**: mapbayr, nlmixr2; **Python**: No mature library |
| SAEM (Stochastic Approximation EM) | Population NLME estimation; no gradient required; robust convergence; primary estimation method for PopPK Studio | **No** | Delyon, Lavielle & Moulines (1999); Kuhn & Lavielle (2005); Lavielle (2014) "Mixed Effects Models for the Population Approach" | **C++**: Monolix (proprietary); **Julia**: Pumas SAEM(); **R**: nlmixr2 saem(), saemix; **C**: NONMEM $EST METHOD=SAEM; **Python**: No mature library |
| Metropolis-Hastings sampler | MCMC sampling of individual random effects within SAEM E-step | **No** | Metropolis et al. (1953); Hastings (1970); Robert & Casella (2004) "Monte Carlo Statistical Methods" | **JS/TS**: Custom (straightforward); **C++**: Stan Math; **Julia**: AdvancedMH.jl; **R**: MCMCpack; **Python**: PyMC, emcee |

### 1.5 Linear Algebra

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Cholesky decomposition | Decompose Ω, Σ positive-definite matrices; sampling multivariate normals; covariance step | **No** | Golub & Van Loan (2013) "Matrix Computations, 4th ed." | **JS/TS**: ml-matrix, stdlib-js; **C**: LAPACK dpotrf, Eigen; **Julia**: cholesky(); **R**: chol(); **Python**: numpy.linalg.cholesky |
| LU decomposition | Solving linear systems; matrix inversion for FOCE-I covariance matrix Cᵢ | **No** | Golub & Van Loan (2013) | **JS/TS**: ml-matrix; **C**: LAPACK, Eigen; **Julia**: lu(); **R**: solve(); **Python**: scipy.linalg.lu |
| Eigendecomposition | Matrix exponential for analytical PK solutions; spectral decomposition of rate matrices | **No** | Golub & Van Loan (2013) | **JS/TS**: ml-matrix; **C**: LAPACK, Eigen; **Julia**: eigen(); **R**: eigen(); **Python**: numpy.linalg.eig |
| Log-determinant | Computing ln|Cᵢ| in FOCE-I objective function; numerically stable via Cholesky | **No** | Via Cholesky: ln|A| = 2·Σ ln(Lᵢᵢ) | **JS/TS**: Via ml-matrix Cholesky; **C**: Eigen; **Julia**: logdet(); **R**: determinant(logarithm=TRUE); **Python**: numpy.linalg.slogdet |

### 1.6 Statistical Methods

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Multivariate normal sampling | Generating random η vectors in SAEM and simulation; requires Cholesky of Ω | **No** | Box & Muller (1958); Gentle (2003) "Random Number Generation and Monte Carlo Methods" | **JS/TS**: stdlib-js + Cholesky; **C**: GSL; **Julia**: Distributions.jl MvNormal; **R**: MASS::mvrnorm; **Python**: numpy.random.multivariate_normal |
| Multivariate normal density | Evaluating p(ηᵢ|Ω) in MAP objective function and likelihood | **No** | Standard textbook | **JS/TS**: Custom from Cholesky; **Julia**: logpdf(MvNormal()); **R**: mvtnorm::dmvnorm; **Python**: scipy.stats.multivariate_normal |
| Seedable PRNG | Reproducible simulation and SAEM sampling; regulatory reproducibility requirement | **No** | Matsumoto & Nishimura (1998) "Mersenne Twister"; Blackman & Vigna (2021) "Xoshiro256" | **JS/TS**: stdlib-js mt19937, seedrandom; **C**: PCG; **Julia**: Random.Xoshiro; **R**: set.seed(); **Python**: numpy.random.Generator |
| Monte Carlo sensitivity analysis | Explore function behavior at randomly sampled points | **Yes** — DiffStudio | Saltelli et al. (2008) "Global Sensitivity Analysis" | **JS/TS**: DiffStudio (Datagrok); **R**: sensitivity; **Julia**: GlobalSensitivity.jl; **Python**: SALib |
| Sobol variance decomposition | Decompose output variance into input contributions; global sensitivity method | **Yes** — DiffStudio | Sobol (1993); Saltelli et al. (2010) "Variance based sensitivity analysis" | **JS/TS**: DiffStudio (Datagrok); **R**: sensitivity::sobol; **Julia**: GlobalSensitivity.jl; **Python**: SALib |
| CWRES (Conditional Weighted Residuals) | Primary diagnostic metric for NLME model adequacy; should follow N(0,1) | **No** | Hooker, Staatz & Karlsson (2007) "Conditional Weighted Residuals (CWRES)" | **Fortran**: NONMEM; **R**: xpose, nlmixr2; **Julia**: Pumas; **Python**: No mature library |
| Shrinkage calculation | Quantify information content of individual parameter estimates; η-shrinkage = 1 − SD(η̂)/ω | **No** | Savic & Karlsson (2009) "Importance of shrinkage in empirical Bayes estimates" | **R**: nlmixr2, xpose; **Julia**: Pumas; **Fortran**: NONMEM |

### 1.7 Diagnostics and Visualization

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| GOF plots | DV vs PRED, DV vs IPRED, residuals vs time/predictions | **Partially** — scatter plot viewers exist | Nguyen et al. (2017) "Model evaluation of continuous data pharmacometric models" | **R**: xpose, ggPMX, nlmixr2; **Julia**: PumasPlots; **Python**: pharmpy |
| QQ-plot | Assess normality of CWRES and ETAs | **Partially** — generic QQ available | Standard statistical reference | **R**: qqnorm(), ggplot2; **Julia**: StatsPlots; **Python**: scipy.stats.probplot |
| Individual fit plots | Observed vs. model-predicted concentration-time profiles per subject | **Partially** — line/scatter viewers | Standard pharmacometric practice | **R**: xpose, nlmixr2; **Julia**: Pumas; **Python**: Custom |

---

## 2. Extended set

These algorithms extend PopPK Studio to competitive parity with NONMEM/Monolix and beyond.

### 2.1 Advanced ODE Solving

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| CVODES (CVODE with forward sensitivity analysis) | CVODE + simultaneous ∂y/∂θ and ∂y/∂η computation; eliminates finite-difference gradients for FOCE-I; 10–100× speedup | **Partially** — CVODE available; forward sensitivity not yet exposed | Serban & Hindmarsh (2005) "CVODES"; Almquist et al. (2015) "Sensitivity equations for FOCE/FOCEI gradients" | **C**: SUNDIALS CVODES (BSD); **Julia**: Sundials.jl; **R**: rxode2; **Python**: Assimulo; **C++**: CasADi; **WASM**: Compilable via Emscripten |
| DDE solvers (Delay Differential Equations) | Models with discrete delays (maturation, transit compartments with explicit delay) | **No** | Baker, Paul & Willé (1995); Bellen & Zennaro (2003) "Numerical Methods for DDEs" | **C**: SUNDIALS IDAS; **Fortran**: RADAR5, NONMEM ADVAN16-18; **Julia**: DDE.jl; **R**: deSolve::dede; **Python**: pydelay |

### 2.2 FOCE-I Infrastructure

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| FOCE-I objective function | Industry-standard population likelihood approximation; -2LL via Laplace + first-order conditional linearization with η-ε interaction | **No** | Wang (2007) "Derivation of various NONMEM estimation methods"; Almquist et al. (2015) | **Fortran**: NONMEM; **C++**: Phoenix NLME; **R**: nlmixr2 focei(); **Julia**: Pumas FOCE() |
| Forward sensitivity equations | Augmented ODE computing ∂y/∂η alongside state y; avoids finite differences; essential for efficient FOCE-I | **No** | Almquist et al. (2015); Dickinson & Gelinas (1976) | **C**: SUNDIALS CVODES; **Julia**: ForwardDiff.jl + DifferentialEquations.jl; **R**: rxode2 (symengine); **Python**: jax.experimental.ode; **C++**: CasADi |
| Sandwich estimator (covariance step) | Robust standard errors: Cov(θ̂) = R⁻¹·S·R⁻¹ where R = Hessian, S = score outer product | **No** | White (1980); Royall (1986); NONMEM $COV | **Fortran**: NONMEM; **R**: nlmixr2; **Julia**: Pumas vcov(); **C++**: Phoenix NLME |
| Finite-difference Hessian | Numerical R-matrix when analytical Hessian unavailable | **No** | Dennis & Schnabel (1996); Gill, Murray & Wright (1981) | **JS/TS**: Custom (central differences); **C**: NLopt; **Julia**: FiniteDifferences.jl; **R**: numDeriv::hessian; **Python**: numdifftools |

### 2.3 Model Qualification

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Non-parametric bootstrap | Resample subjects, re-estimate → parameter CIs; gold standard for uncertainty | **No** | Efron (1979); Parke et al. (1999) "Population PK bootstrap" | **Perl**: PsN bootstrap; **R**: nlmixr2 bootstrapFit(); **Julia**: Pumas bootstrap() |
| VPC (Visual Predictive Check) | Simulate 500–1000 datasets → compare percentile bands; primary graphical model qualification | **No** | Bergstrand et al. (2011) "Prediction-corrected visual predictive checks"; Holford (2005) | **R**: vpc, tidyvpc, xpose; **Julia**: Pumas vpc(); **Perl**: PsN vpc |
| pcVPC (Prediction-corrected VPC) | VPC normalized by population prediction; removes covariate-driven variability | **No** | Bergstrand et al. (2011) CPT:PSP | **R**: vpc::vpc(pred_corr=TRUE), tidyvpc; **Julia**: Pumas; **Perl**: PsN vpc -predcorr |
| NPDE (Normalized Prediction Distribution Errors) | Simulation-based diagnostic; model-agnostic; NPDE ~ N(0,1) if correct | **No** | Brendel et al. (2006); Comets et al. (2008) | **R**: npde, nlmixr2; **Julia**: Pumas |
| SIR (Sampling Importance Resampling) | Fast alternative to bootstrap; single estimation + importance sampling for parameter uncertainty | **No** | Dosne et al. (2016) "Improving estimation of parameter uncertainty" | **Perl**: PsN sir; **R**: nlmixr2 |

### 2.4 Covariate Analysis

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| SCM (Stepwise Covariate Modeling) | Automated forward inclusion (ΔOFV > 3.84) / backward elimination (ΔOFV > 6.63) | **No** | Jonsson & Karlsson (1998); Ribbing et al. (2007) | **Perl**: PsN scm; **R**: nlmixr2; **Julia**: Pumas covariate_select() |
| FREM (Full Random Effects Model) | Joint modeling of parameters and covariates; avoids stepwise selection bias | **No** | Karlsson (2012); Jonsson et al. (2024) "FREM practical usage guide" | **Perl**: PsN frem; **R**: nlmixr2 (planned); **Fortran**: NONMEM |
| GAM/CART screening | Generalized Additive Models or CART trees for initial covariate screening | **No** | Jonsson & Karlsson (1998); Mandema et al. (1992) | **R**: mgcv::gam, rpart; **Julia**: GAM.jl; **Python**: pygam, sklearn |
| Forest plot generation | Visualize covariate effects with CIs; regulatory-standard figure | **No** | FDA Guidance on Population PK (2022); ISoP best practices | **R**: ggplot2, forestly, nlmixr2; **Julia**: Custom; **Python**: forestplot |

### 2.5 Simulation Methods

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Monte Carlo simulation with variability | Sample θ, η ~ N(0,Ω), ε ~ N(0,Σ) → simulate virtual populations | **Partially** — MC sensitivity exists; population-level η/ε sampling not yet | Mould & Upton (2012); Holford et al. (2010) | **R**: rxode2, mrgsolve, nlmixr2; **Julia**: Pumas simobs(); **Fortran**: NONMEM $SIMULATION |
| Monte Carlo simulation with uncertainty | Additionally sample θ from N(θ̂, Cov(θ̂)) and Ω from Wishart; captures parameter uncertainty | **No** | Dosne et al. (2016); Plan et al. (2012) | **R**: rxode2; **Julia**: Pumas; **Perl**: PsN sse |
| Target attainment probability | P(AUC > MIC target) or P(Cmax < safety threshold) from simulated distributions | **No** | Drusano (2004); Mouton et al. (2005) | **R**: Custom; **Julia**: Custom; **Python**: Custom |

### 2.6 Bayesian Methods (server-side)

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| NUTS (No-U-Turn Sampler) | Adaptive HMC; gold standard for full Bayesian NLME; requires AD through ODE solver | **No** | Hoffman & Gelman (2014); Margossian et al. (2022) "Bayesian pharmacometrics with Stan and Torsten" | **C++**: Stan; **Fortran**: NONMEM $EST METHOD=NUTS; **Julia**: Pumas BayesMCMC(), Turing.jl; **R**: rstan; **Python**: PyMC |
| HMC (Hamiltonian Monte Carlo) | Efficient MCMC using gradient information; core of modern Bayesian inference | **No** | Neal (2011); Betancourt (2017) "A Conceptual Introduction to HMC" | **C++**: Stan Math; **Julia**: AdvancedHMC.jl; **R**: rstan; **Python**: PyMC, BlackJAX |
| Torsten PK functions | Stan-compatible PK solvers with NONMEM-style event handling | **No** | Margossian et al. (2022); Weber et al. (2018) "Stan functions for pharmacometrics" | **C++**: Torsten (Metrum Research Group, BSD) |
| Automatic differentiation | Exact gradients through ODE solvers; required for HMC/NUTS | **No** | Griewank & Walther (2008) "Evaluating Derivatives" | **C++**: Stan Math, Adept, CoDiPack, CasADi; **Julia**: ForwardDiff.jl, Zygote.jl, Enzyme.jl; **Python**: JAX, PyTorch; **JS**: No mature library |

### 2.7 Advanced Estimation

| Algorithm | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| Importance sampling (post-SAEM likelihood) | Compute -2LL after SAEM; needed for model comparison (AIC/BIC) and FOCE-I OFV equivalence | **No** | Geweke (1989); NONMEM $EST METHOD=IMP; Monolix linearization | **Fortran**: NONMEM IMP, IMPMAP; **C++**: Monolix; **R**: nlmixr2; **Julia**: Pumas |
| MCEM (Monte Carlo EM) | Alternative to SAEM; MC integration in E-step with increasing sample sizes | **No** | Wei & Tanner (1990) "A Monte Carlo implementation of the EM algorithm" | **Julia**: Pumas MCEM(); **R**: Custom |
| BOBYQA | Derivative-free constrained optimization; useful when gradients are unreliable | **No** | Powell (2009) "BOBYQA algorithm" | **C**: NLopt; **Julia**: NLopt.jl; **R**: nloptr; **Python**: nlopt, pybobyqa |
| Stochastic gradient descent (Adam) | Potential acceleration of outer loop in large-dataset scenarios | **No** | Kingma & Ba (2014) "Adam" | **JS/TS**: various; **Julia**: Flux.jl; **Python**: PyTorch, TensorFlow |

### 2.8 Interoperability

| Algorithm / Method | Purpose | Datagrok | Key references | Reliable libraries |
|-----------|---------|----------|----------------|-------------------|
| NMTRAN parser | Read/write NONMEM-format datasets and control streams | **No** | Beal et al. NONMEM Users Guide; Bauer (2019) "NONMEM Tutorial Part I" | **R**: nonmem2rx, babelmixr2, NMdata; **Julia**: Pumas; **Python**: pharmpy; **Perl**: PsN |
| Mlxtran parser | Import Monolix model files | **No** | Monolix documentation (Lixoft) | **R**: babelmixr2 (planned); **Python**: pharmpy (partial) |
| SBML (Systems Biology Markup Language) | Import/export models in standard format; QSP interoperability | **No** | Hucka et al. (2003); Keating et al. (2020) "SBML Level 3" | **JS**: libsbml.js (WASM); **C++**: libSBML; **Julia**: SBMLToolkit.jl; **R**: rsbml; **Python**: python-libsbml |

---

## Summary matrix

| Category | Minimal | Extended | Total | Already in Datagrok |
|----------|---------|----------|-------|-------------------|
| ODE Solvers | 5 | 2 | 7 | **5** (all minimal) |
| Analytical PK | 2 | — | 2 | 0 |
| Optimization | 4 | 1 | 5 | **1** (Nelder-Mead) |
| NLME Estimation | 3 | 4 | 7 | 0 |
| Linear Algebra | 4 | — | 4 | 0 |
| Statistics & Sensitivity | 7 | 1 | 8 | **3** (MC, Sobol, Grid) |
| Model Qualification | — | 5 | 5 | 0 |
| Covariate Analysis | — | 4 | 4 | 0 |
| Simulation | — | 3 | 3 | 0 (partial MC) |
| Bayesian | — | 4 | 4 | 0 |
| Diagnostics | 3 | — | 3 | **~1.5** (partial viewers) |
| Interoperability | — | 3 | 3 | 0 |
| **Total** | **28** | **27** | **55** | **~10.5** |

**Coverage: Datagrok provides ~19% of the required algorithms.** The strongest coverage is in ODE solving (100% of minimal set — including LSODA and CVODE from SUNDIALS) and sensitivity analysis (3 of 3 methods). The primary gaps are NLME estimation (0/7), linear algebra (0/4), model qualification (0/5), and interoperability (0/3).
