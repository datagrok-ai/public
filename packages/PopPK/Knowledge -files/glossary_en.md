# PopPK Studio — Glossary

| Term | Definition |
|------|-----------|
| ADVAN | NONMEM subroutine defining the PK model structure (e.g., ADVAN2 = one-compartment oral) |
| AIC (Akaike Information Criterion) | Model selection metric balancing goodness of fit and complexity; lower is better |
| AMT | Amount of dose administered (NMTRAN data column) |
| AUC (Area Under the Curve) | Total drug exposure over time, calculated as the integral of the concentration-time curve |
| BID | Twice daily dosing (bis in die) |
| BLQ / BQL (Below Quantification Limit) | Observed concentration below the assay's limit of detection |
| Bootstrap | Resampling method: repeatedly re-estimate the model on resampled datasets to obtain parameter confidence intervals |
| BPV (Between-Patient Variability) | See IIV |
| CDISC (Clinical Data Interchange Standards Consortium) | Standards for clinical trial data formatting (SDTM, ADaM) |
| CL (Clearance) | Volume of plasma completely cleared of drug per unit time |
| Cmax | Maximum observed or predicted drug concentration |
| CMT | Compartment number in NMTRAN data (identifies which compartment a dose or observation applies to) |
| Compartmental model | Mathematical model describing the body as interconnected compartments (e.g., central, peripheral) |
| Control stream | NONMEM input file specifying model, data, and estimation settings in NM-TRAN syntax |
| Covariate | Patient characteristic (weight, age, GFR, genotype) that may explain variability in PK/PD parameters |
| CRO (Contract Research Organization) | Company providing outsourced pharmaceutical research services |
| CSR (Clinical Study Report) | Comprehensive regulatory document summarizing a clinical trial's design, conduct, and results |
| Ctrough | Drug concentration at the end of a dosing interval (just before next dose) |
| CWRES (Conditional Weighted Residuals) | Diagnostic metric: standardized residuals conditional on individual parameter estimates; should follow N(0,1) |
| DSL (Domain-Specific Language) | Specialized programming language designed for a specific problem domain |
| DV (Dependent Variable) | Observed data value in NMTRAN format (typically drug concentration) |
| E-R (Exposure-Response) | Analysis linking drug exposure metrics to efficacy or safety endpoints |
| EBE (Empirical Bayes Estimate) | Individual parameter estimate obtained by combining population information with individual data |
| eCTD (Electronic Common Technical Document) | Standardized format for regulatory submissions to FDA/EMA/PMDA |
| eGFR (Estimated Glomerular Filtration Rate) | Kidney function metric used as a covariate in PK models |
| Emax | Maximum drug effect in the Emax pharmacodynamic model |
| EOP2 (End of Phase 2) | FDA meeting to discuss Phase 2 results and Phase 3 design |
| Eta (η) | Random effect representing individual deviation from the population parameter |
| EVID | Event ID in NMTRAN data (0=observation, 1=dose, 2=other, 3=reset, 4=reset+dose) |
| Fixed effects (θ) | Population-average model parameters estimated from data |
| FO (First-Order) | Simplest NLME approximation method; linearizes around η=0 |
| FOCE (First-Order Conditional Estimation) | NLME method linearizing around individual conditional η estimates |
| FOCE-I (FOCE with Interaction) | FOCE accounting for η–ε interaction in the residual variance model |
| Forest plot | Visualization showing covariate effects on PK parameters with confidence intervals |
| GOF (Goodness of Fit) | Diagnostic plots assessing how well the model describes observed data |
| GxP | Umbrella term for regulatory quality practices (GLP, GCP, GMP, etc.) |
| ICH M15 | International guideline on general principles for Model-Informed Drug Development |
| IIV (Inter-Individual Variability) | Between-subject variability in PK/PD parameters, described by the Ω matrix |
| IOV (Inter-Occasion Variability) | Within-subject, between-occasion variability in PK/PD parameters |
| IPRED (Individual Prediction) | Model prediction using individual (post-hoc) parameter estimates |
| ITT (Intent-To-Treat) | Analysis population including all randomized subjects regardless of compliance |
| IWRES (Individual Weighted Residuals) | Residuals scaled by individual-level residual variance estimate |
| Ka (Absorption rate constant) | First-order rate constant for drug absorption from the GI tract |
| MAP (Maximum A Posteriori) | Bayesian estimation of individual parameters using population priors |
| MCMC (Markov Chain Monte Carlo) | Sampling-based Bayesian inference method |
| MDV (Missing Dependent Variable) | Flag in NMTRAN data indicating the observation should not contribute to the likelihood |
| MIDD (Model-Informed Drug Development) | FDA framework for using quantitative models to inform drug development decisions |
| MMRM (Mixed Model for Repeated Measures) | Linear mixed-effects model for longitudinal continuous endpoints in clinical trials |
| M&S (Modeling and Simulation) | Umbrella term for quantitative model-based analysis in drug development |
| NCA (Non-Compartmental Analysis) | Model-independent PK analysis calculating AUC, Cmax, t½ directly from data |
| NDA (New Drug Application) | Regulatory submission to FDA for marketing approval of a new drug |
| NLME (Nonlinear Mixed-Effects Model) | Statistical model combining fixed and random effects with a nonlinear structural model |
| NMTRAN | NONMEM translator: converts user-friendly model code to NONMEM's internal format; also refers to the standard data format |
| NONMEM | NONlinear Mixed-Effects Modeling — the original and most widely used pharmacometrics software |
| NPDE (Normalized Prediction Distribution Errors) | Simulation-based diagnostic; if model is correct, NPDE ~ N(0,1) |
| ODE (Ordinary Differential Equation) | Mathematical equation describing rate of change; used to define compartmental PK/PD models |
| OFV (Objective Function Value) | -2 × log-likelihood; minimized during NLME estimation; used for model comparison (lower = better fit) |
| Omega (Ω) | Variance-covariance matrix of inter-individual random effects |
| pcVPC (Prediction-Corrected VPC) | VPC variant normalizing data by population prediction to remove covariate-driven variability |
| PD (Pharmacodynamics) | Study of drug effects on the body (what the drug does to the body) |
| PK (Pharmacokinetics) | Study of drug absorption, distribution, metabolism, and excretion (what the body does to the drug) |
| PopPK (Population Pharmacokinetics) | Application of NLME models to characterize PK variability across a patient population |
| PRED (Population Prediction) | Model prediction using population-average (typical) parameters (η=0) |
| PREDPP | NONMEM library of pre-programmed PK model subroutines |
| PsN (Perl-speaks-NONMEM) | Automation toolkit for NONMEM providing bootstrap, VPC, SCM, and cluster computing |
| QD | Once daily dosing (quaque die) |
| QSP (Quantitative Systems Pharmacology) | Mechanistic modeling approach integrating systems biology with pharmacology |
| Random effects (η, ε) | Model components describing unexplained variability (between-subject and residual) |
| REML (Restricted Maximum Likelihood) | Estimation method for variance parameters that accounts for degrees of freedom lost in estimating fixed effects |
| SAEM (Stochastic Approximation EM) | NLME estimation via iterating between MCMC sampling of individual parameters and closed-form population updates |
| SAP (Statistical Analysis Plan) | Pre-specified document describing all planned analyses for a clinical trial |
| SCM (Stepwise Covariate Modeling) | Automated covariate selection via forward inclusion then backward elimination based on OFV changes |
| SDTM (Study Data Tabulation Model) | CDISC standard for organizing clinical trial raw data |
| Shrinkage | Degree to which individual parameter estimates are pulled toward the population mean; high shrinkage (>30%) reduces reliability of individual-level diagnostics |
| Sigma (Σ) | Parameter(s) describing residual unexplained variability |
| Stiff ODE | ODE system with widely separated time scales requiring specialized solvers (BDF, Rosenbrock) |
| TDM (Therapeutic Drug Monitoring) | Clinical practice of measuring drug concentrations to optimize individual dosing |
| Theta (θ) | Fixed-effect population parameter (e.g., typical CL, typical V) |
| TMDD (Target-Mediated Drug Disposition) | PK model for drugs whose disposition is significantly affected by binding to their pharmacological target |
| TRANS | NONMEM parameterization translator (e.g., TRANS2: CL,V; TRANS4: CL,V1,Q,V2) |
| V (Volume of distribution) | Apparent volume relating drug amount in the body to plasma concentration |
| VPC (Visual Predictive Check) | Model qualification plot comparing observed data percentiles against simulated prediction intervals |
| WASM (WebAssembly) | Binary instruction format enabling near-native performance in web browsers |
