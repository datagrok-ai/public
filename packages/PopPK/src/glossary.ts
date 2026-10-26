export const glossary = new Map<string, string>([
  // NMTRAN data columns
  ['AMT', 'Amount of dose administered (NMTRAN data column)'],
  ['CMT', 'Compartment number in NMTRAN data (identifies which compartment a dose or observation applies to)'],
  ['DV', 'Observed data value in NMTRAN format (typically drug concentration)'],
  ['EVID', 'Event ID in NMTRAN data (0=observation, 1=dose, 2=other, 3=reset, 4=reset+dose)'],
  ['MDV', 'Flag in NMTRAN data indicating the observation should not contribute to the likelihood'],

  // PK parameters
  ['AUC', 'Total drug exposure over time, calculated as the integral of the concentration-time curve'],
  ['CL', 'Volume of plasma completely cleared of drug per unit time'],
  ['Cmax', 'Maximum observed or predicted drug concentration'],
  ['Ctrough', 'Drug concentration at the end of a dosing interval (just before next dose)'],
  ['Ka', 'First-order rate constant for drug absorption from the GI tract'],
  ['V', 'Apparent volume relating drug amount in the body to plasma concentration'],

  // Estimation methods
  ['FO', 'Simplest NLME approximation method; linearizes around \u03B7=0'],
  ['FOCE', 'NLME method linearizing around individual conditional \u03B7 estimates'],
  ['FOCE-I', 'FOCE accounting for \u03B7\u2013\u03B5 interaction in the residual variance model'],
  ['MAP', 'Bayesian estimation of individual parameters using population priors'],
  ['SAEM', 'NLME estimation via iterating between MCMC sampling of individual parameters and closed-form population updates'],
  ['EBE', 'Individual parameter estimate obtained by combining population information with individual data'],
  ['NLME', 'Statistical model combining fixed and random effects with a nonlinear structural model'],
  ['OFV', '-2 \u00D7 log-likelihood; minimized during NLME estimation; used for model comparison (lower = better fit)'],
  ['REML', 'Estimation method for variance parameters that accounts for degrees of freedom lost in estimating fixed effects'],
  ['MCMC', 'Sampling-based Bayesian inference method'],

  // Random effects & variability
  ['Eta (\u03B7)', 'Random effect representing individual deviation from the population parameter'],
  ['Theta (\u03B8)', 'Fixed-effect population parameter (e.g., typical CL, typical V)'],
  ['Omega (\u03A9)', 'Variance-covariance matrix of inter-individual random effects'],
  ['Sigma (\u03A3)', 'Parameter(s) describing residual unexplained variability'],
  ['Fixed effects (\u03B8)', 'Population-average model parameters estimated from data'],
  ['Random effects (\u03B7, \u03B5)', 'Model components describing unexplained variability (between-subject and residual)'],
  ['IIV', 'Between-subject variability in PK/PD parameters, described by the \u03A9 matrix'],
  ['IOV', 'Within-subject, between-occasion variability in PK/PD parameters'],
  ['BPV', 'See IIV'],
  ['Shrinkage', 'Degree to which individual parameter estimates are pulled toward the population mean; high shrinkage (>30%) reduces reliability of individual-level diagnostics'],

  // Diagnostics
  ['CWRES', 'Diagnostic metric: standardized residuals conditional on individual parameter estimates; should follow N(0,1)'],
  ['IWRES', 'Residuals scaled by individual-level residual variance estimate'],
  ['GOF', 'Diagnostic plots assessing how well the model describes observed data'],
  ['IPRED', 'Model prediction using individual (post-hoc) parameter estimates'],
  ['PRED', 'Model prediction using population-average (typical) parameters (\u03B7=0)'],
  ['VPC', 'Model qualification plot comparing observed data percentiles against simulated prediction intervals'],
  ['pcVPC', 'VPC variant normalizing data by population prediction to remove covariate-driven variability'],
  ['NPDE', 'Simulation-based diagnostic; if model is correct, NPDE ~ N(0,1)'],
  ['Forest plot', 'Visualization showing covariate effects on PK parameters with confidence intervals'],

  // Software & tools
  ['NONMEM', 'NONlinear Mixed-Effects Modeling \u2014 the original and most widely used pharmacometrics software'],
  ['NMTRAN', 'NONMEM translator: converts user-friendly model code to NONMEM\'s internal format; also refers to the standard data format'],
  ['ADVAN', 'NONMEM subroutine defining the PK model structure (e.g., ADVAN2 = one-compartment oral)'],
  ['TRANS', 'NONMEM parameterization translator (e.g., TRANS2: CL,V; TRANS4: CL,V1,Q,V2)'],
  ['PREDPP', 'NONMEM library of pre-programmed PK model subroutines'],
  ['PsN', 'Automation toolkit for NONMEM providing bootstrap, VPC, SCM, and cluster computing'],
  ['WASM', 'Binary instruction format enabling near-native performance in web browsers'],

  // Modeling concepts
  ['Compartmental model', 'Mathematical model describing the body as interconnected compartments (e.g., central, peripheral)'],
  ['Control stream', 'NONMEM input file specifying model, data, and estimation settings in NM-TRAN syntax'],
  ['Covariate', 'Patient characteristic (weight, age, GFR, genotype) that may explain variability in PK/PD parameters'],
  ['Bootstrap', 'Resampling method: repeatedly re-estimate the model on resampled datasets to obtain parameter confidence intervals'],
  ['SCM', 'Automated covariate selection via forward inclusion then backward elimination based on OFV changes'],
  ['NCA', 'Model-independent PK analysis calculating AUC, Cmax, t\u00BD directly from data'],
  ['ODE', 'Mathematical equation describing rate of change; used to define compartmental PK/PD models'],
  ['Stiff ODE', 'ODE system with widely separated time scales requiring specialized solvers (BDF, Rosenbrock)'],
  ['DSL', 'Specialized programming language designed for a specific problem domain'],
  ['PopPK', 'Application of NLME models to characterize PK variability across a patient population'],
  ['QSP', 'Mechanistic modeling approach integrating systems biology with pharmacology'],
  ['TMDD', 'PK model for drugs whose disposition is significantly affected by binding to their pharmacological target'],
  ['Emax', 'Maximum drug effect in the Emax pharmacodynamic model'],
  ['E-R', 'Analysis linking drug exposure metrics to efficacy or safety endpoints'],
  ['M&S', 'Umbrella term for quantitative model-based analysis in drug development'],
  ['MMRM', 'Linear mixed-effects model for longitudinal continuous endpoints in clinical trials'],
  ['TDM', 'Clinical practice of measuring drug concentrations to optimize individual dosing'],
  ['PK', 'Study of drug absorption, distribution, metabolism, and excretion (what the body does to the drug)'],
  ['PD', 'Study of drug effects on the body (what the drug does to the body)'],

  // Regulatory & clinical
  ['AIC', 'Model selection metric balancing goodness of fit and complexity; lower is better'],
  ['BLQ / BQL', 'Observed concentration below the assay\'s limit of detection'],
  ['CDISC', 'Standards for clinical trial data formatting (SDTM, ADaM)'],
  ['CRO', 'Company providing outsourced pharmaceutical research services'],
  ['CSR', 'Comprehensive regulatory document summarizing a clinical trial\'s design, conduct, and results'],
  ['eCTD', 'Standardized format for regulatory submissions to FDA/EMA/PMDA'],
  ['eGFR', 'Kidney function metric used as a covariate in PK models'],
  ['EOP2', 'FDA meeting to discuss Phase 2 results and Phase 3 design'],
  ['GxP', 'Umbrella term for regulatory quality practices (GLP, GCP, GMP, etc.)'],
  ['ICH M15', 'International guideline on general principles for Model-Informed Drug Development'],
  ['ITT', 'Analysis population including all randomized subjects regardless of compliance'],
  ['MIDD', 'FDA framework for using quantitative models to inform drug development decisions'],
  ['NDA', 'Regulatory submission to FDA for marketing approval of a new drug'],
  ['SAP', 'Pre-specified document describing all planned analyses for a clinical trial'],
  ['SDTM', 'CDISC standard for organizing clinical trial raw data'],

  // Dosing
  ['BID', 'Twice daily dosing (bis in die)'],
  ['QD', 'Once daily dosing (quaque die)'],

  // Personas
  ['PMx', 'Pharmacometrician\n"I don\'t need another GUI that hides the math. I need the same power as NONMEM but without the Fortran-era workflow. Show me the OFV, show me the gradients, let me script everything \u2014 but in a modern environment."'],
  ['CP', 'Clinical Pharmacologist\n"I don\'t want to learn NONMEM syntax. I want to pick a two-compartment model, point it at my data, and explore what happens when I change the dose for a 40 kg pediatric patient \u2014 all before my next meeting."'],
  ['BS', 'Biostatistician\n"I spend two days every E-R analysis just merging the pharmacometrician\'s individual exposures with my efficacy dataset and making sure the joins are correct. If it were all in one platform with traceable data lineage, I\'d save a week per study."'],
  ['PM', 'Program Leadership\n"I don\'t care about the OFV or the omega matrix. I care about one thing: if I pick 200 mg for Phase 3, what\'s the probability that we hit the efficacy target in the intent-to-treat population? Show me that \u2014 live, now, in this meeting."'],
]);
