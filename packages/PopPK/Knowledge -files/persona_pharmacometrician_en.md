# Persona: Pharmacometrician (PMx)

## Persona Card

| Attribute | Details |
|-----------|---------|
| **Role** | Pharmacometrician / Population PK-PD Modeler |
| **Department** | Clinical Pharmacology & Pharmacometrics |
| **Organization** | Big Pharma, Mid-size Biotech, CRO |
| **Team size** | 5–30 per company (Big Pharma); 2–5 (Biotech); 10–50 (large CRO) |
| **Education** | PhD in Pharmaceutical Sciences, Pharmacology, Applied Mathematics, or Statistics |
| **Experience** | 3–15 years in pharmacometrics |
| **Reports to** | Head of Clinical Pharmacology or VP Quantitative Pharmacology |
| **Budget authority** | Recommends software purchases; final approval at director/VP level |

### Profile

The pharmacometrician is the primary power user of PopPK Studio. They build population PK/PD models end-to-end: from structural model selection and random effects specification through covariate analysis to model qualification and regulatory reporting. They think in terms of OFV drops, eta-shrinkage percentages, and CWRES distributions. They spend 60–80% of their working time in NONMEM + R, switching between 4–6 disconnected tools per analysis.

### Goals

- Build scientifically defensible PopPK/PD models that withstand regulatory scrutiny
- Deliver analyses on time for regulatory submissions and clinical decision milestones
- Reproduce and extend models from prior studies or literature
- Mentor junior modelers and establish best practices within the team

### Current tools

NONMEM, Monolix, R (nlmixr2, xpose, vpc, ggPMX), PsN, Pirana, SAS (occasionally), Git, Word/PowerPoint for reports

### Pain points

- **Fragmented workflow**: juggling 4–6 tools (NONMEM → R → PsN → Excel → Word) with manual data transfer between each
- **Debugging hell**: hours spent troubleshooting NONMEM control stream syntax errors with cryptic Fortran error messages
- **No collaboration**: sharing models via email/shared drives; no version control by default; "which control stream is the final one?"
- **Slow iteration**: waiting for cluster job queues; no interactive exploration of model behavior
- **Reproducibility burden**: manually documenting every step for regulatory audit trails

### What they need from PopPK Studio

- Code-level control over model specification (DSL with syntax highlighting, autocomplete)
- Full estimation suite: FOCE-I, SAEM, MAP — matching NONMEM OFV to 3+ significant figures
- Integrated diagnostics (GOF, VPC, CWRES) without switching to R
- NONMEM interoperability: import/export control streams
- Version-controlled, auditable analysis workflow
- Performance: estimation in seconds, not minutes

### Quote

> "I don't need another GUI that hides the math. I need the same power as NONMEM but without the Fortran-era workflow. Show me the OFV, show me the gradients, let me script everything — but in a modern environment."

---

## User Journey Map

### Scenario: PopPK analysis for a Phase 2b oral drug, regulatory submission in 8 weeks

#### Stage 1 — Data Assembly (Day 1–3)

| Aspect | Details |
|--------|---------|
| **Task** | Combine dosing records, PK samples, and covariates into analysis-ready dataset |
| **Actions** | Import SDTM/ADaM datasets → map to NMTRAN format → handle BQL records (M3 method) → QC checks (duplicates, time sequences, missing covariates) → generate data summary statistics |
| **Current experience** | Manual R scripts or SAS programs; error-prone merging of multiple CSV files; hours spent on EVID/CMT coding |
| **PopPK Studio experience** | Drag-and-drop data import → auto-detect NMTRAN columns → visual BQL flagging → built-in QC rules with one-click fix suggestions → data summary dashboard |
| **Emotion** | Current: frustration (tedious, error-prone) → PopPK Studio: confidence (validated, visual) |
| **Key metric** | Time from raw data to analysis-ready dataset |

#### Stage 2 — Base Model Development (Day 3–10)

| Aspect | Details |
|--------|---------|
| **Task** | Select structural PK model, define random effects hierarchy, choose residual error model |
| **Actions** | Test 1-cmt vs 2-cmt vs 3-cmt → first-order vs transit absorption → IIV on CL, V, Ka → proportional vs combined error → compare OFV, AIC, diagnostics |
| **Current experience** | Edit NONMEM control stream in text editor → submit to cluster → wait → check output → edit → resubmit (10–30 iterations) |
| **PopPK Studio experience** | Select model from library or write in DSL → run estimation in browser (seconds for analytical models) → instant GOF plots → compare models in side-by-side view → iterate interactively |
| **Emotion** | Current: patience-testing (slow loop) → PopPK Studio: flow state (instant feedback) |
| **Key metric** | Number of models explored per day |

#### Stage 3 — Covariate Analysis (Day 10–20)

| Aspect | Details |
|--------|---------|
| **Task** | Identify patient factors explaining inter-individual variability |
| **Actions** | ETA-vs-covariate plots → GAM/SCM screening → forward inclusion (ΔOFV > 3.84) → backward elimination (ΔOFV > 6.63) → forest plots of covariate effects |
| **Current experience** | PsN SCM automation (hours to days) or manual stepwise in NONMEM; separate R scripts for forest plots |
| **PopPK Studio experience** | Interactive ETA-covariate explorer → automated stepwise with real-time progress → integrated forest plot generator → model comparison table auto-updates |
| **Emotion** | Current: anxiety (did I miss a covariate? Is the order correct?) → PopPK Studio: systematic confidence |
| **Key metric** | Time for complete covariate analysis cycle |

#### Stage 4 — Model Qualification (Day 20–30)

| Aspect | Details |
|--------|---------|
| **Task** | Demonstrate model adequacy through bootstrap, VPC, and sensitivity analysis |
| **Actions** | Bootstrap 500–1000 replicates → compare CIs with point estimates → generate pcVPC → run sensitivity analyses on key assumptions |
| **Current experience** | PsN bootstrap submitted to cluster (overnight to days); VPC via vpc R package; manual assembly of results |
| **PopPK Studio experience** | One-click bootstrap (server-side parallelized) → interactive VPC with adjustable binning → integrated sensitivity analysis dashboard → all results linked to model version |
| **Emotion** | Current: waiting (cluster queues, overnight runs) → PopPK Studio: productive (results by lunch) |
| **Key metric** | Wall time from final model to qualified model |

#### Stage 5 — Simulation & Dose Recommendation (Day 30–40)

| Aspect | Details |
|--------|---------|
| **Task** | Simulate alternative dosing regimens, generate dose recommendations for target populations |
| **Actions** | Define virtual populations (covariates sampled from distributions) → simulate 1000+ subjects per regimen → compute target attainment → generate exposure-response plots → compare regimens |
| **Current experience** | Custom R scripts using mrgsolve or rxode2; manual assembly of simulation scenarios |
| **PopPK Studio experience** | Visual simulation builder → define scenarios via GUI or script → parallel Monte Carlo in WebWorkers → interactive dose-response explorer → export to report template |
| **Emotion** | Current: tedious (lots of scripting for each scenario) → PopPK Studio: empowering (interactive exploration) |
| **Key metric** | Number of scenarios explored; time to dose recommendation |

#### Stage 6 — Regulatory Reporting (Day 40–56)

| Aspect | Details |
|--------|---------|
| **Task** | Prepare PopPK report for Module 2.7.2 of eCTD submission |
| **Actions** | Write methods section → generate tables (parameter estimates, bootstrap CIs) → insert diagnostic figures → document model building decisions → QC review → finalize |
| **Current experience** | Copy-paste from R/NONMEM output into Word; manual table formatting; figures as PNG exports; version control nightmare |
| **PopPK Studio experience** | Auto-generated report template with live-linked tables and figures → edit narrative in platform → one-click PDF/Word export → complete audit trail → reviewer can trace any number to its source |
| **Emotion** | Current: dread (the most hated part) → PopPK Studio: relief (automated, traceable) |
| **Key metric** | Report preparation time; number of QC errors caught |
