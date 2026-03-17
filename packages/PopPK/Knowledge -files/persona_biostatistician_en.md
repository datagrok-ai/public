# Persona: Biostatistician (BS)

## Persona Card

| Attribute | Details |
|-----------|---------|
| **Role** | Biostatistician / Statistical Scientist |
| **Department** | Biostatistics, Statistical Programming |
| **Organization** | Big Pharma, CRO |
| **Team size** | 50–200 per company (Big Pharma); 20–100 (large CRO) |
| **Education** | PhD or MS in Biostatistics, Statistics, or Mathematics |
| **Experience** | 3–20 years in clinical trial statistics |
| **Reports to** | Head of Biostatistics or VP Biometrics |
| **Budget authority** | Influences statistical software purchases (SAS, R, East) |

### Profile

The biostatistician designs clinical trials, writes statistical analysis plans (SAPs), and executes primary efficacy/safety analyses. Their world centers on MMRM, ANCOVA, survival analysis, and sample size calculations — all in SAS or R. They intersect with pharmacometrics at two points: exposure-response (E-R) analyses that link drug exposure to clinical outcomes, and MMRM analyses of longitudinal endpoints where PopPK model-derived exposure metrics are used as covariates. They are rigorous about statistical methodology but often view pharmacometric models as a parallel discipline with unfamiliar tools and terminology.

### Goals

- Deliver statistically rigorous trial analyses on regulatory timelines
- Ensure proper handling of missing data, multiplicity, and estimand frameworks (ICH E9(R1))
- Collaborate effectively with pharmacometricians on E-R analyses without tool fragmentation
- Adopt reproducible, validated workflows that satisfy GxP requirements

### Current tools

SAS (PROC MIXED, PROC GLIMMIX, PROC PHREG), R (lme4, mmrm, survival), East/nQuery (sample size), Pinnacle 21 (CDISC), Git/SVN

### Pain points

- **Two worlds, zero integration**: MMRM in SAS and PopPK in NONMEM are completely disconnected — different data formats, different software, different teams, different reports
- **Exposure-response handoff friction**: receiving PopPK-derived exposure metrics (AUC, Cmax, Ctrough) requires manual data merging with efficacy datasets
- **Reproducibility across tools**: results must be traceable from raw SDTM data through analysis to regulatory tables — cross-tool pipelines break audit trails
- **Learning curve barrier**: would benefit from understanding PopPK diagnostics but the NONMEM learning cliff is too steep for occasional use

### What they need from PopPK Studio

- Seamless import of PopPK-derived individual exposure metrics into their analysis datasets
- MMRM capability within the same platform (or clean R integration for lme4/mmrm)
- Exposure-response analysis module that bridges PopPK output and statistical endpoints
- Validated, GxP-ready computational environment with audit trail
- Familiar R/SAS-like interface — not a foreign DSL

### Quote

> "I spend two days every E-R analysis just merging the pharmacometrician's individual exposures with my efficacy dataset and making sure the joins are correct. If it were all in one platform with traceable data lineage, I'd save a week per study."

---

## User Journey Map

### Scenario: Exposure-response analysis for Phase 3 primary efficacy endpoint

#### Stage 1 — Receive Exposure Metrics (Day 1–3)

| Aspect | Details |
|--------|---------|
| **Task** | Obtain individual exposure metrics (AUC, Cmax, Ctrough) from PopPK model |
| **Actions** | Receive CSV of individual post-hoc estimates from PMx → validate subject IDs → merge with ADEFF dataset → QC the merge (missing subjects, duplicates) |
| **Current experience** | Email from PMx with CSV attachment; manual merge in SAS; hours of QC checking |
| **PopPK Studio experience** | Exposure metrics auto-generated from PopPK model → available as linked dataset in Datagrok → one-click merge with efficacy data → visual merge validation dashboard |
| **Emotion** | Current: tedious (manual data plumbing) → PopPK Studio: effortless (auto-linked) |
| **Key metric** | Time from model completion to analysis-ready E-R dataset |

#### Stage 2 — Exploratory E-R Analysis (Day 3–7)

| Aspect | Details |
|--------|---------|
| **Task** | Explore relationship between drug exposure and efficacy/safety endpoints |
| **Actions** | Scatter plots of AUC vs. response → exposure quartile analysis → logistic regression (binary endpoints) → preliminary MMRM with exposure covariate |
| **Current experience** | SAS programs or R scripts; manual generation of each plot; separate code for each endpoint |
| **PopPK Studio experience** | Integrated E-R explorer → auto-generate exposure quartile plots → point-and-click logistic/linear regression → MMRM with exposure as time-varying covariate |
| **Emotion** | Current: routine programming → PopPK Studio: interactive exploration |
| **Key metric** | Number of E-R relationships explored; time to initial signal |

#### Stage 3 — Formal E-R Modeling (Day 7–15)

| Aspect | Details |
|--------|---------|
| **Task** | Build definitive exposure-response model for regulatory submission |
| **Actions** | Fit Emax / linear / logistic E-R models → covariate effects on E-R → model selection (AIC/BIC) → sensitivity analyses → generate regulatory-ready tables and figures |
| **Current experience** | SAS PROC NLMIXED or R nls/nlme; separate code for each model; manual table generation |
| **PopPK Studio experience** | E-R model library → side-by-side model comparison → automated sensitivity analysis → one-click regulatory table generation with data traceability |
| **Emotion** | Current: careful but slow → PopPK Studio: same rigor, faster execution |
| **Key metric** | Time to final E-R model; regulatory table preparation time |

#### Stage 4 — Cross-Functional Review & Reporting (Day 15–20)

| Aspect | Details |
|--------|---------|
| **Task** | Align E-R results with PMx team and prepare joint Clinical Pharmacology/Statistics section |
| **Actions** | Review E-R findings with PMx → ensure consistency between PopPK report and E-R report → contribute to integrated summary → finalize CSR statistical tables |
| **Current experience** | Meetings to reconcile different tool outputs; manual cross-referencing of results; separate report sections |
| **PopPK Studio experience** | Shared workspace with PMx → E-R results directly linked to PopPK model → unified report generation → single source of truth for all numbers |
| **Emotion** | Current: friction (reconciling parallel workflows) → PopPK Studio: alignment (shared platform) |
| **Key metric** | Number of reconciliation cycles needed; cross-team alignment time |
