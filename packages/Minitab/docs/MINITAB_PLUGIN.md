# Minitab Plugin for Datagrok — Reference Documentation

> **Purpose:** This document is the primary reference for AI-assisted development of the Datagrok Minitab plugin. It covers everything an agent needs to understand the domain, make sound architectural decisions, and implement features with appropriate context.

---

## Table of Contents

1. [Minitab Overview](#1-minitab-overview)
2. [Pharma Industry Context](#2-pharma-industry-context)
3. [Competitive Landscape](#3-competitive-landscape)
4. [Core Feature Set](#4-core-feature-set)
5. [File Formats](#5-file-formats)
6. [Technological Fit for Datagrok](#6-technological-fit-for-datagrok)
7. [Suggested Solution Workflows](#7-suggested-solution-workflows)
8. [Implementation Roadmap](#8-implementation-roadmap)

---

## 1. Minitab Overview

### What Is Minitab?

Minitab is a statistical analysis and process improvement platform originally developed at Pennsylvania State University in 1972. It is the dominant GUI-based statistics tool for **quality engineering**, **manufacturing**, and **Six Sigma/Lean** programs. As of 2026, Minitab LLC is headquartered in State College, Pennsylvania, with subsidiaries in the UK, France, Germany, Netherlands, Hong Kong, Japan, and Australia.

The product suite has expanded significantly through acquisitions:
- **2017:** Salford Systems (CART, Random Forests, TreeNet — predictive analytics)
- **2024:** Simul8 (discrete-event simulation)
- **2025:** Prolink Software (automated quality data collection, SPC connectivity)
- **2026:** Scytec / DataXchange (real-time OEE and machine data from CNC, PLCs, robotics)

The **Minitab Solution Center** (cloud platform) consolidates: Minitab Statistical Software, Minitab Engage (project management), Minitab Workspace (visual tools), Minitab Connect (data connectors), Real-Time SPC, Salford Predictive Modeler, and the Minitab Education Hub.

### Typical Users

| Persona | Role | Primary Use |
|---|---|---|
| Quality Engineer / QE | QA/QC team | SPC charts, capability analysis, Gage R&R, CAPA |
| Process Engineer | Manufacturing/CMC | DOE, regression, process validation |
| Six Sigma Black/Green Belt | Continuous improvement | DMAIC projects, hypothesis tests, MSA |
| R&D Scientist | Pharma/biotech | Stability studies, assay development, method validation |
| Manufacturing Ops | Plant floor | Real-time SPC dashboards, OEE tracking |
| Statistician/Analyst | Support role | Advanced regression, multivariate analysis |
| Academic | University | Teaching statistics (very large installed base) |

### Typical Workflows

1. **Quality Control:** Import batch data → run control charts (Xbar-R, I-MR, P, U) → assess process capability (Cpk, Ppk) → export report to Word/PowerPoint
2. **Process Validation (FDA 3-stage):** Design stage DOE → qualification Gage R&R → continued process verification SPC
3. **Stability Studies:** Import time-series assay data → regression with batch as random effect → shelf-life estimation
4. **Method Validation:** Linearity, accuracy, precision (repeatability, reproducibility) → measurement uncertainty
5. **DMAIC Projects:** Minitab Workspace for process mapping → Minitab Stats for Measure/Analyze phases → Engage for tracking

### Strengths and Weaknesses

**Strengths:**
- Dominant brand recognition in quality/manufacturing — the "Excel of SPC"
- Exceptional ease of use for non-statisticians; "Assistant" feature walks users through analyses
- 50+ years of methodological reliability; FDA validation kit available
- Strong training ecosystem (thousands of certified Lean Six Sigma courses reference Minitab)

**Weaknesses:**
- Limited data visualization vs. modern BI tools (Tableau, Spotfire)
- Desktop-first; web app exists but is less capable
- No native connection to enterprise data sources (LIMS, ERP) without Minitab Connect add-on
- Expensive for broad deployment (≈ $1,851/user/year)
- No collaborative or multi-user project environment
- Session-based workflow — not reproducible/auditable in a modern sense
- Limited scripting (proprietary macro language, not Python/R first-class)

---

## 2. Pharma Industry Context

### Market Position

In a pharma SPC survey by ECA Academy, **Minitab was cited by 44% of respondents** as their primary SPC tool — far ahead of any other dedicated package. The "Other" category (41%, including Excel) and LIMS-integrated tools (6%) were the only significant alternatives. This makes Minitab the de facto standard for pharma QC/manufacturing statistics.

### Key Pharma Use Cases

#### FDA Process Validation (21 CFR §211.100, §211.165)
The FDA three-stage process validation lifecycle maps directly to Minitab capabilities:

| Stage | FDA Goal | Minitab Tools Used |
|---|---|---|
| Stage 1 — Process Design | Characterize CQAs and CPPs | DOE (factorial, response surface), regression, ANOVA |
| Stage 2 — Process Qualification | Demonstrate process control | Gage R&R / MSA, capability analysis (Cpk, Ppk), hypothesis tests |
| Stage 3 — Continued Process Verification | Ongoing state of control | SPC control charts, capability monitoring, trend analysis |

#### Pharmaceutical CMC Statistics
- **Assay/method validation:** linearity, accuracy, precision, specificity (ICH Q2)
- **Stability studies:** regression with batch factors, shelf-life calculation (ICH Q1E)
- **Dissolution testing:** capability analysis for tablet dissolution profiles
- **In-process controls:** real-time monitoring during manufacturing

#### FDA Computer System Validation (21 CFR Part 11)
Pharma companies must validate any software that generates regulatory data. Minitab provides a **Validation Kit** (installation qualification scripts, OQ test protocols) that companies use to satisfy CSV requirements. This is a critical procurement factor — any replacement tool must also have a CSV story.

#### GMP Data Integrity (FDA 21 CFR Part 11 / EMA Annex 11)
Data must be: attributable, legible, contemporaneous, original, accurate (ALCOA). Minitab session logs and project files provide an audit trail, though it is not a validated 21 CFR Part 11 system. Datagrok's provenance tracking is actually superior here.

### Pharma Personas Relevant to Datagrok

The Minitab user in pharma is typically in **QA/QC, manufacturing, or CMC** — not in discovery or informatics. This is a distinct persona from Datagrok's typical pharma user (cheminformatics, DMPK, HTS). The Minitab plugin creates an opportunity to expand Datagrok's reach into:
- QC labs (analytical chemistry, microbiology)
- CMC groups (formulation, process development)
- Manufacturing operations (plant floor, batch release)
- Regulatory affairs (process validation packages)

---

## 3. Competitive Landscape

### Direct Competitors

| Tool | Vendor | Positioning | Price (approx.) | Key Differentiator |
|---|---|---|---|---|
| **Minitab Statistical Software** | Minitab LLC | Quality/SPC/Six Sigma leader | ~$1,851/user/year | Ease of use, brand dominance, FDA validation kit |
| **JMP** | SAS Institute | Engineers/scientists, visual analytics | ~$1,320/user/year (perpetual-style) | Richer interactive visualization, R integration, scripting (JSL) |
| **IBM SPSS Statistics** | IBM | Social sciences, clinical research | ~$1,069/user/year | Strong for survey data, clinical trials, GLM |
| **SAS Viya / Base SAS** | SAS Institute | Enterprise analytics, pharma regs | Request quote (enterprise) | Regulatory gold standard; dominates clinical/NDA submissions |
| **SigmaXL** | SigmaXL Inc. | Excel add-in, Six Sigma training | ~$280 perpetual per user | Lowest cost; leverages Excel familiarity |
| **QI Macros** | KnowWare Intl. | Excel add-in, SPC | ~$369 perpetual per user | Instant SPC inside Excel; no new software to learn |
| **XLSTAT** | Lumivero | Excel add-in, advanced stats | ~$400-800/year | Broader statistical methods than SigmaXL; pharma module |
| **GraphPad Prism** | Dotmatics | Life sciences, biological data | ~$800/user/year | Biology-focused; curve fitting, survival, dose-response |
| **Statgraphics** | Statgraphics Tech. | SPC + DOE | ~$800-1,200/year | Strong SPC + DOE combo; less brand recognition |
| **R + RStudio/Posit** | Open source / Posit | Data science, pharma analytics | Free (Posit Workbench priced separately) | Unlimited capability; requires programming |
| **Python** | Open source | General purpose | Free | Unlimited capability; requires programming |
| **numiqo** | Numiqo | Minitab alternative, web-based | Lower cost SaaS | Emerging competitor; targets Minitab switchers |

### Market Size

- **SPC Software market:** ~$943M in 2024, growing at ~12% CAGR to ~$2.1B by 2031
- **Statistical Software market (broad):** ~$11.3B in 2024, growing at ~9% CAGR to ~$23.8B by 2033
- **Minitab** has ~3,874 active job listings (significant for a quality-focused tool), behind only the largest platforms

### Positioning Opportunity for Datagrok

Datagrok enters this space not as a direct Minitab replacement but as a **migration path and augmentation platform**:

- Minitab users are increasingly frustrated with: cost, siloed desktop files, lack of enterprise data connectivity, poor visualization, and no collaboration
- JMP is the most common upgrade path from Minitab — Datagrok can compete here too
- The key differentiator: Datagrok offers Minitab's core statistical outputs **inside a connected, collaborative, browser-native platform** with live data from LIMS, ELN, and databases

---

## 4. Core Feature Set

These are the Minitab capabilities that Minitab users most rely on, ordered by frequency of use in pharma/manufacturing contexts. The plugin should prioritize implementing or mapping these to Datagrok equivalents.

### Tier 1: Must-Have (Most Used)

#### Statistical Process Control (SPC)
- **Control charts:** Xbar-R (subgroup mean and range), I-MR (individuals and moving range), Xbar-S (mean and std dev), P chart (proportion nonconforming), NP chart (number nonconforming), C chart (count), U chart (count per unit), CUSUM, EWMA
- **Chart rules:** Western Electric (WECO) rules, Nelson rules — detect special-cause variation
- **Out-of-control signals:** Automatic flagging of rule violations on charts
- **Laney charts:** Overdispersion-adjusted P' and U' charts (important for pharma)

#### Process Capability Analysis
- **Normal capability:** Cp, Cpk, Pp, Ppk, Cpm with confidence intervals
- **Non-normal capability:** Distribution fitting, Box-Cox / Johnson transformation
- **Automated capability distribution** (recent Minitab feature) — auto-selects best distribution
- **Capability sixpack:** Combined I-MR, run chart, histogram, normal probability plot, and capability indices in one panel — very frequently used in pharma QC

#### Measurement System Analysis (MSA / Gage R&R)
- Crossed Gage R&R (operators × parts × replicates)
- Nested Gage R&R
- Attribute MSA (Cohen's kappa, effectiveness)
- EMP method (Wheeler's preferred criterion)
- Bias, linearity, stability studies

#### Hypothesis Tests
- 1-sample and 2-sample t-tests, paired t-test
- 1-sample and 2-sample proportions
- Chi-square test of association
- 1-sample and 2-sample Poisson rate
- Equivalence tests (TOST) — critical for pharma method comparisons
- Normality tests (Anderson-Darling, Ryan-Joiner, Kolmogorov-Smirnov)

### Tier 2: Highly Important

#### Regression and ANOVA
- Simple and multiple linear regression with diagnostics
- Binary, ordinal, and nominal logistic regression
- One-way, two-way, and general linear model (GLM) ANOVA
- Analysis of means (ANOM)
- Fitted-line plots

#### Design of Experiments (DOE)
- Full and fractional factorial designs (2k and 2k-p)
- General full factorial (multi-level)
- Response surface designs (central composite, Box-Behnken)
- Mixture designs
- Definitive Screening Designs (DSD)
- Optimal custom designs
- Contour plots, response optimizer, overlaid contour plots

#### Reliability / Survival Analysis
- Distribution ID plot (fit multiple distributions)
- Probability plots (Weibull, lognormal, exponential)
- Parametric and nonparametric survival analysis
- Accelerated life testing

#### Time Series
- Time series plots with trend lines
- Moving averages, exponential smoothing
- Autocorrelation (ACF, PACF)
- ARIMA modeling and forecasting

### Tier 3: Supplementary

#### Multivariate Analysis
- Principal component analysis (PCA)
- Cluster analysis (hierarchical, k-means)
- Discriminant analysis
- Factor analysis (with pharma applications in spectroscopy and formulation)

#### Power and Sample Size
- Sample size calculators for all hypothesis test types
- Power curves

#### Stability Studies (Pharma-specific)
- Regression with batch as random/fixed effect
- Shelf-life estimation at specified confidence level
- Poolability tests (test whether batches can be pooled)

#### The Assistant
- Guided analysis wizard — walks users through choosing and running the correct statistical test
- Generates summary cards, diagnostic cards, and report cards in plain language
- Critical for non-statistician users (most Minitab pharma users are engineers, not statisticians)

---

## 5. File Formats

### Minitab Native Formats

#### `.mpx` — Minitab Project (current, Minitab 19+)
- **Structure:** ZIP archive containing JSON and/or XML files
- **Contents:**
  - `/project_metadata.json` — project name, Minitab version, creation date
  - `/sheets/0/sheet.json` — worksheet data (columns, types, values, missing value coding)
  - `/commands/1/command.json` — analysis commands run in session
  - `/commands/1/groups/2/group.json` — output groups (graphs, tables, session text)
  - Embedded graph data (SVG or proprietary format)
- **Parse approach:** Rename to `.zip`, extract, parse JSON. No official public spec but format is reverse-engineerable. Average file: 110 KB; range: 2 KB – several MB for large projects.
- **Priority:** HIGH — this is the primary export format users will have

#### `.mwx` — Minitab Worksheet (current, Minitab 19+)
- **Structure:** ZIP-based (same family as MPX), contains worksheet data only — no graphs or session output
- **Contents:** Column data (numeric, text, date/time), column names, column descriptions, value labels, missing values (`*` for numeric, blank for text)
- **Priority:** HIGH — simplest path to data import; users often share worksheets without full projects

#### `.mpj` — Minitab Project (legacy, Minitab 15–18)
- **Structure:** OLE2 Compound Document (same format as legacy `.doc`, `.xls`)
- **Magic bytes:** `D0 CF 11 E0 A1 B1 1A E1`
- **Internal streams:** `UI`, `Worksheets`, `Compressed`, `ProjectInfo`, `Charts`
- **CLSID:** `de9268e3-30f8-11d0-9a73-00a024c92ff4`
- **Parse approach:** Use `cfb` (Compound File Binary) library. Complex; lower priority than MPX.
- **Priority:** MEDIUM — many legacy files exist in pharma (long-lived validation packages)

#### `.mtw` — Minitab Worksheet (legacy, up to Minitab 18)
- **Structure:** OLE2 (same as MPJ)
- **Contents:** Worksheet data only — equivalent to `.mwx` but older format
- **Priority:** MEDIUM

#### `.mgf` — Minitab Graph File (all versions)
- **Structure:** Binary/proprietary
- **Contents:** A single saved graph
- **Priority:** LOW — graphs are better re-created than imported

#### `.mtb` — Minitab Exec Script (all versions)
- **Structure:** Plain text
- **Contents:** Series of Minitab session commands; used for automation
- **Priority:** LOW initially, but HIGH value for "bring your macros" migration story

#### `.mac` — Minitab Macro (all versions)
- **Structure:** Plain text, begins with `GMACRO` or `MACRO`
- **Contents:** Minitab command language with control flow (IF/DO/WHILE/CALL)
- **Priority:** LOW for initial release; MEDIUM for advanced migration scenarios

### Minitab Session Command Language

All Minitab GUI operations generate session commands (visible in the History pane) that can be saved and replayed. This is the key to macro automation. Examples:

```
HISTOGRAM C1;
  NINTERVAL 10;
  BAR.

CAPABILITY C1 C2;
  USPEC 100;
  LSPEC 80;
  TARGET 90.
```

Understanding this language enables potential future features: importing Minitab analysis history and recreating analyses natively in Datagrok.

### Other Relevant Import Formats

These are the formats Minitab users commonly work with alongside `.mpx`/`.mwx`:

| Format | Notes |
|---|---|
| Excel (`.xlsx`, `.xls`) | Most common data exchange format; Minitab imports/exports natively |
| CSV / TXT / DAT | Simple column-per-variable text files; Minitab imports with delimiter detection |
| XML | Minitab can export worksheets to XML |
| HTML | Minitab can export session output as HTML tables |
| ODBC | Minitab Connect provides ODBC connectivity to databases/LIMS |

---

## 6. Technological Fit for Datagrok

### Alignment

Datagrok is exceptionally well-positioned to absorb Minitab workflows because:

1. **Column-oriented data model:** Minitab's worksheet is fundamentally columns of typed data — identical to Datagrok's DataFrame. Column types (numeric, text, date/time) map 1:1. Missing value conventions (`*` → `null`) are trivial to handle.

2. **Browser-native architecture:** Minitab's main weakness — desktop-first, file-based, no collaboration — is Datagrok's core strength. The plugin delivers an immediate "why switch" story.

3. **Statistical functions:** Datagrok's JS/R/Python script execution environment can host any statistical computation. R packages (`qcc`, `SixSigma`, `capability`, `MASS`, `lme4`) cover all Minitab statistical methods. The plugin can either wrap R scripts or implement statistical functions natively in JavaScript using WebAssembly (faster, no server round-trip).

4. **Visualization infrastructure:** Datagrok's viewer framework can implement all Minitab chart types (control charts, capability histograms, probability plots, residual plots) as custom viewers with proper statistical annotations.

5. **Scripting bridge:** Minitab's `.mtb`/`.mac` command language is a simple imperative language. A transpiler from Minitab session commands to Datagrok/R/Python scripts is architecturally feasible.

6. **Regulatory readiness:** Datagrok's audit trail, data provenance, and access control are superior to Minitab's native capabilities — a selling point for FDA 21 CFR Part 11 compliance.

### Gaps to Address

| Gap | Minitab Capability | Datagrok Path |
|---|---|---|
| Statistical output formatting | Minitab produces rich formatted output cards (summary, diagnostic, report) | Custom output panels for each analysis; implement "report card" pattern |
| "Assistant" guided workflow | Step-by-step analysis wizard for non-statisticians | AI-powered analysis assistant (leverages Datagrok's LLM integration) |
| Capability sixpack | Single-panel combining 6 plots + indices | Composite viewer layout |
| Control chart rules engine | Western Electric / Nelson rules, visual flagging | Implement rules engine in JS or R; flag points in viewer |
| DOE design generation | Optimal design algorithms, design catalog | R packages (`AlgDesign`, `rsm`, `FrF2`) or integrate Minitab's recently acquired Effex catalog |
| Gage R&R ANOVA output | Structured ANOVA tables with %contribution, %tolerance | R `lme4` or `qcc` package; output formatter |
| Shelf-life calculation | Regression + ICH Q1E method | R `stability` or custom ICH-compliant implementation |
| FDA validation kit | IQ/OQ scripts for CSV | Generate Datagrok-equivalent qualification documentation |

---

## 7. Suggested Solution Workflows

### 7.1 Data Import Workflows

#### Import `.mpx` Project File
**User story:** "I have a Minitab project from our QC lab. I want to continue working with the data in Datagrok."

**Flow:**
1. User drags `.mpx` file onto Datagrok or uses File > Open
2. Plugin detects `.mpx` signature (ZIP + `project_metadata.json`)
3. Parser extracts all worksheets as Datagrok DataFrames (one per Minitab worksheet)
4. Column metadata (name, type, description, value labels) is preserved
5. User sees a "Minitab Import Summary" panel:
   - Worksheets imported (with row/column counts)
   - Session commands found (analysis history — read-only, shown for reference)
   - Graphs found (count; offer to re-create natively or skip)
6. Each worksheet opens in a Datagrok table view
7. Optional: Session history shown as a read-only audit log panel

**Parser implementation notes:**
- Use `JSZip` to unzip MPX in the browser
- Parse `sheet.json` for column definitions and data
- Map Minitab types: `Numeric` → `float64`, `Text` → `string`, `Date/Time` → `datetime`
- Minitab missing value `*` → Datagrok `null`

#### Import `.mwx` Worksheet File
- Simpler path: same as MPX but only one worksheet, no session history or graphs
- Fast path for data exchange between Minitab and Datagrok users

#### Import Legacy `.mpj` / `.mtw` Files
- Use `cfb` library to parse OLE2 compound document
- Extract `Worksheets` stream; parse binary worksheet data
- Medium complexity; prioritize after MPX/MWX support

### 7.2 Analysis Workflows

#### SPC / Control Charts Workflow
**User story:** "Our QC analyst wants to monitor a critical quality attribute across batches."

**Flow:**
1. User selects a numeric column + optional subgroup column in Datagrok
2. Right-click → Quality Tools → Control Chart (or from the ribbon)
3. Dialog: select chart type (I-MR, Xbar-R, P, U, etc.), configure subgroup size, limits
4. Chart renders as a Datagrok viewer with:
   - Data points plotted over time/batch sequence
   - UCL / CL / LCL lines
   - Out-of-control points highlighted with rule labels (e.g., "Rule 1: Point beyond 3σ")
   - WECO / Nelson rule selector in chart properties
5. Summary panel alongside: process statistics, number of violations
6. Chart is live — updates when new data is added to the table

**Implementation:** Custom Datagrok viewer using D3 or Canvas; R `qcc` package for limit calculation and rule testing.

#### Process Capability Workflow
**User story:** "Formulation team wants to show Cpk > 1.33 for a critical parameter before submitting the process validation report."

**Flow:**
1. Select column with process measurements
2. Quality Tools → Capability Analysis
3. Dialog: enter USL, LSL, target; select distribution (Normal / non-normal / auto-detect)
4. Output panel (replicates Minitab's "Capability Sixpack"):
   - I-MR chart (top left)
   - Capability histogram with spec limits (top center)
   - Normal probability plot (top right)
   - Last 25 subgroup ranges (bottom left)
   - Run chart (bottom center)
   - Summary statistics panel with Cp, Cpk, Pp, Ppk, % out of spec
5. One-click export: copy to clipboard / export to Word report

**Implementation:** Composite multi-viewer layout; capability indices computed via R `SixSigma` or native JS.

#### Gage R&R / MSA Workflow
**User story:** "Before validating an HPLC method, we need to qualify the measurement system."

**Flow:**
1. User imports a Gage R&R data worksheet (columns: Part, Operator, Measurement)
2. Quality Tools → Measurement System Analysis → Gage R&R (Crossed)
3. Dialog: identify Part, Operator, Measurement columns; enter study limits
4. Output:
   - ANOVA table: Source, DF, SS, MS, F, P
   - Variance components: Gage (repeatability), Reproducibility (operator), Part-to-Part, Total
   - %Study Variation (%SV) and %Tolerance for each component
   - Six-panel graphical output: Components of Variation, R chart by operator, Xbar chart by operator, By-Part plot, By-Operator plot, Operator-Part interaction plot
5. Pass/fail summary: %GRR < 10% = acceptable; 10–30% = marginal; >30% = unacceptable

**Implementation:** R `lme4` for ANOVA; custom multi-panel viewer.

#### DOE (Design of Experiments) Workflow
**User story:** "Process development team is optimizing a fermentation process — 3 factors at 2 levels."

**Flow:**
1. Quality Tools → Design of Experiments → Create Design
2. Wizard:
   - Step 1: Choose design type (2-level factorial, RSM, etc.)
   - Step 2: Enter factor names, levels, units
   - Step 3: Select full/fractional factorial, number of centerpoints, blocks, replicates
   - Step 4: Preview design summary and generate
3. Datagrok creates a run-order DataFrame as a structured worksheet
4. User conducts experiments, enters response values into the Datagrok table
5. Quality Tools → DOE → Analyze Design
6. Output: ANOVA, effects plots (main effects, interaction plots), Pareto chart of effects, residual diagnostic plots, response optimizer

**Implementation:** R packages `FrF2` (factorial), `rsm` (response surface), `AlgDesign` (optimal).

#### Stability Study Workflow (Pharma-specific)
**User story:** "Regulatory affairs needs a shelf-life calculation from our stability data per ICH Q1E."

**Flow:**
1. Import stability data: columns = Batch, Time_months, Assay_pct
2. Quality Tools → Stability Analysis
3. Dialog: specify assay column, time column, batch column, acceptance criteria (e.g., 90% of label claim)
4. Tests: test for poolability of batches (ANCOVA), select pooled vs. unpooled model
5. Regression with batch as fixed effect; extrapolate to 90% lower confidence bound
6. Output: regression plots per batch + combined, shelf-life estimate with 95% confidence
7. Generate ICH Q1E-formatted summary table for regulatory submission

### 7.3 Migration and Compatibility Workflows

#### Bring Your Macros (`.mtb` / `.mac` Conversion)
**User story:** "Our QC lab has 30 Minitab macros that run monthly reports. We want to migrate them to Datagrok."

**Flow:**
1. User uploads `.mtb` / `.mac` files to Datagrok
2. Plugin parses Minitab command language
3. Generates equivalent Datagrok/R/Python script (best-effort transpilation)
4. Shows diff view: "Minitab command → Datagrok equivalent"
5. User reviews, adjusts, saves as Datagrok script function
6. Macros become Datagrok Functions visible in the function browser

**Scope note:** Full transpilation is ambitious. V1 can produce an annotated script with comments explaining each Minitab command and the equivalent Datagrok/R approach, leaving the user to finalize.

#### Side-by-Side Validation Mode
**User story:** "Regulatory team needs to verify that Datagrok produces the same numerical results as Minitab before decommissioning Minitab."

**Flow:**
1. User imports `.mpx` file containing analysis results
2. Plugin runs the same analyses natively in Datagrok
3. Comparison table: Minitab result vs. Datagrok result for each statistic
4. Tolerance: flag any difference > 0.001 (numerical precision)
5. Generate validation report (PDF) for IQ/OQ documentation

This feature directly addresses the CSV concern and positions Datagrok as a validated Minitab replacement.

#### Export Back to Minitab
**User story:** "I've cleaned data in Datagrok but my colleague still uses Minitab."

**Flow:**
1. Any Datagrok DataFrame → File > Save As > Minitab Worksheet (.mwx)
2. Generate valid MWX ZIP structure with correct JSON schema
3. Preserves column names, types, and values
4. Column metadata (descriptions, value labels) included where supported

### 7.4 Visualization Workflows

#### Re-creating Minitab Graphs as Datagrok Viewers

When a `.mpx` file contains saved graphs, the plugin can offer to **re-create** them as live Datagrok viewers rather than importing static images:

| Minitab Graph | Datagrok Implementation |
|---|---|
| Histogram | Native histogram viewer + distribution fit overlay |
| Scatterplot + regression line | Scatter viewer + R regression script |
| Time series plot | Line viewer with time axis |
| Boxplot by group | Box plot viewer |
| Individual value plot | Strip/dot plot viewer |
| Control chart | Custom SPC viewer (new) |
| Capability histogram | Custom capability viewer (new) |
| Normal probability plot | Custom probability plot viewer (new) |
| Residual plots | Composite residuals viewer (new) |
| Pareto chart | Custom Pareto viewer (new) |
| Main effects / interaction plots | Custom DOE viewers (new) |

The advantage over static images: these viewers are **live** — they update when data changes, support filtering and selection, and can be shared via Datagrok links.

### 7.5 Reporting Workflows

#### Export Analysis Results to Word / PDF

Minitab users are accustomed to pasting results and graphs into Word documents for SOPs, validation reports, and regulatory submissions.

**Flow:**
1. Any Datagrok analysis panel → right-click → Export to Word / Export to PDF
2. Include: analysis title, configuration, statistical output table, graph(s), interpretation notes
3. For pharma regulatory use: include data source, version, timestamp, analyst name
4. Template system: pre-configured Word templates for common report types (process validation, stability, Gage R&R)

---

## 8. Implementation Roadmap

### Phase 1 — Data Bridge (MVP)

**Goal:** Make it frictionless to get Minitab data into Datagrok.

| Feature | Priority | Effort |
|---|---|---|
| `.mpx` file parser → DataFrames | P0 | Medium |
| `.mwx` worksheet parser | P0 | Low |
| Column type and metadata preservation | P0 | Low |
| Session history viewer (read-only) | P1 | Low |
| Export DataFrame as `.mwx` | P1 | Low |
| Legacy `.mpj` / `.mtw` parser | P2 | High |

### Phase 2 — Core Statistical Analyses

**Goal:** Replicate the 80% of Minitab analyses that pharma users run most.

| Feature | Priority |
|---|---|
| Control charts (I-MR, Xbar-R, P, NP, C, U, CUSUM, EWMA) | P0 |
| WECO / Nelson rules engine | P0 |
| Capability analysis (normal + auto-distribution) | P0 |
| Capability Sixpack composite viewer | P0 |
| Hypothesis tests (t-tests, proportions, equivalence) | P1 |
| Normality tests | P1 |
| Gage R&R / MSA | P1 |
| Basic regression + diagnostics | P1 |
| One-way / Two-way ANOVA | P2 |

### Phase 3 — Advanced Workflows

| Feature | Priority |
|---|---|
| DOE design generation + analysis | P1 |
| Stability study / shelf-life (ICH Q1E) | P1 |
| Multivariate (PCA, cluster) | P2 |
| Power and sample size | P2 |
| Reliability / distribution analysis | P2 |
| Report templates (Word/PDF) | P2 |

### Phase 4 — Migration Tooling

| Feature | Priority |
|---|---|
| Side-by-side validation mode | P1 |
| Macro transpiler (`.mtb` → R/Python) | P2 |
| CSV/IQ/OQ documentation generator | P2 |
| Guided "Assistant"-style analysis wizard (AI-powered) | P3 |

---

*Last updated: 2026. See also: `MINITAB_FILE_FORMATS.md` for detailed format specifications and `MINITAB_STATISTICAL_METHODS.md` for algorithm references.*
