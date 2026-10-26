# Persona: Clinical Pharmacologist (CP)

## Persona Card

| Attribute | Details |
|-----------|---------|
| **Role** | Clinical Pharmacologist |
| **Department** | Clinical Pharmacology, Clinical Development |
| **Organization** | Big Pharma, Mid-size Biotech |
| **Team size** | 10–50 per company (Big Pharma); 3–10 (Biotech) |
| **Education** | PharmD, MD, or PhD in Pharmacology/Clinical Pharmacology |
| **Experience** | 5–20 years in drug development |
| **Reports to** | VP Clinical Pharmacology or Chief Medical Officer |
| **Budget authority** | Influences software and CRO budgets at the program level |

### Profile

The clinical pharmacologist is the strategic decision-maker who translates modeling results into clinical actions. They own the Clinical Pharmacology strategy for a drug program: dose selection rationale, special population recommendations (renal/hepatic impairment, pediatrics, elderly), drug-drug interaction assessments, and the Clinical Pharmacology sections of the drug label. They commission analyses from pharmacometricians but increasingly want to run basic models and simulations themselves — especially "what if?" scenarios during team meetings.

### Goals

- Make data-driven dose recommendations that survive regulatory review
- Anticipate and preemptively answer FDA/EMA questions about PK variability and dose justification
- Communicate complex PK/PD concepts to clinical teams, medical affairs, and executives
- Reduce dependency on pharmacometricians for exploratory analyses

### Current tools

Phoenix WinNonlin (primarily NCA), Excel, R (basic), PowerPoint, forest plots from pharmacometrician's output, clinical trial databases

### Pain points

- **Bottleneck dependency**: every "what if?" question requires a pharmacometrician to modify code, re-run, and deliver results — turnaround days to weeks
- **Black box frustration**: receives model results but cannot explore the model interactively, change a covariate assumption, or run a quick simulation
- **Translation gap**: struggles to connect OFV-level model details to clinical questions ("so what does this mean for a 90 kg patient with CKD stage 3?")
- **NCA-to-PopPK cliff**: comfortable with NCA in Phoenix but the jump to NONMEM-based PopPK requires entirely different skills and tools
- **Report assembly**: manually combines pharmacometrician output with clinical interpretation for regulatory documents

### What they need from PopPK Studio

- GUI-driven workflow: point-and-click model selection, estimation, and diagnostics
- Pre-built model library with visual compartmental diagrams and plain-language descriptions
- Interactive simulation sandbox: adjust dose, body weight, renal function → see predicted PK profile instantly
- Forest plot generator for covariate effects with clinical interpretation aids
- Seamless NCA-to-PopPK transition within the same platform
- Exportable figures and tables formatted for regulatory documents

### Quote

> "I don't want to learn NONMEM syntax. I want to pick a two-compartment model, point it at my data, and explore what happens when I change the dose for a 40 kg pediatric patient — all before my next meeting."

---

## User Journey Map

### Scenario: Preparing dose justification for an End-of-Phase-2 meeting with FDA

#### Stage 1 — Review Available Data (Day 1–2)

| Aspect | Details |
|--------|---------|
| **Task** | Understand the PK data landscape across Phase 1 and Phase 2a studies |
| **Actions** | Review NCA summary from Phase 1 → examine concentration-time profiles → identify variability patterns → check covariate coverage (weight, age, renal function ranges) |
| **Current experience** | Request data summaries from biostatistician; review Phoenix NCA output; manually browse clinical database |
| **PopPK Studio experience** | Open integrated data explorer → visual concentration-time overlay by study/dose → covariate distribution dashboard → automatic identification of sparse vs. rich sampling subjects |
| **Emotion** | Current: passive (waiting for others' output) → PopPK Studio: engaged (exploring data directly) |
| **Key metric** | Time to form initial PK assessment |

#### Stage 2 — Commission or Run PopPK Analysis (Day 2–10)

| Aspect | Details |
|--------|---------|
| **Task** | Obtain population PK model for the compound |
| **Actions** | If PMx available: brief them on objectives, wait for results. If self-service: select model from library → run estimation → review diagnostics |
| **Current experience** | Write analysis request → wait 1–3 weeks for PMx team → receive results package → back-and-forth clarifications |
| **PopPK Studio experience** | Open model library → select "2-compartment oral, first-order absorption" → auto-configure from data → click "Estimate (SAEM)" → review auto-generated diagnostics → refine if needed or escalate to PMx for complex issues |
| **Emotion** | Current: frustration (dependency, delays) → PopPK Studio: autonomy (self-service for 70% of cases) |
| **Key metric** | Days from analysis need to first results |

#### Stage 3 — Explore Dose-Exposure-Response (Day 10–18)

| Aspect | Details |
|--------|---------|
| **Task** | Evaluate dose-exposure relationships and identify optimal dose for Phase 2b |
| **Actions** | Simulate multiple dosing regimens → overlay predicted AUC/Cmax distributions → map to efficacy and safety thresholds → assess impact of covariates on target attainment |
| **Current experience** | Request specific simulation scenarios from PMx → wait → receive static plots → request modifications → wait again |
| **PopPK Studio experience** | Open simulation sandbox → drag sliders (dose: 100→200→400 mg; weight: 50→120 kg; eGFR: 30→120 mL/min) → real-time PK profile updates → toggle efficacy/safety thresholds → capture scenarios |
| **Emotion** | Current: impatience (can't iterate fast enough) → PopPK Studio: discovery (exploring the design space) |
| **Key metric** | Number of dose scenarios evaluated; time to dose recommendation |

#### Stage 4 — Build Dose Justification Narrative (Day 18–25)

| Aspect | Details |
|--------|---------|
| **Task** | Construct the clinical pharmacology argument for proposed Phase 2b dosing |
| **Actions** | Generate forest plots of covariate effects → create dose-exposure-response summary → prepare special population recommendations → draft Clinical Pharmacology section |
| **Current experience** | Collect figures from PMx, NCA tables from Phoenix, clinical data from database → manually assemble in PowerPoint/Word |
| **PopPK Studio experience** | Generate forest plot with one click → export dose-response figures directly → auto-populate regulatory template with model-linked tables → edit narrative in platform |
| **Emotion** | Current: tedious assembly work → PopPK Studio: streamlined storytelling |
| **Key metric** | Report preparation time; consistency between figures and narrative |

#### Stage 5 — Present to Cross-Functional Team & FDA (Day 25–30)

| Aspect | Details |
|--------|---------|
| **Task** | Present dose rationale to internal governance and prepare for FDA EOP2 meeting |
| **Actions** | Present to internal review board → address questions with live model exploration → finalize briefing document → present at FDA meeting |
| **Current experience** | Static PowerPoint; if questioned on "what about dose X for population Y?" → "we'll get back to you" |
| **PopPK Studio experience** | Present from interactive dashboard → respond to ad-hoc questions with live simulations during the meeting → export updated scenarios immediately |
| **Emotion** | Current: anxiety (unprepared for ad-hoc questions) → PopPK Studio: confidence (can answer in real time) |
| **Key metric** | Ability to answer ad-hoc questions; meeting outcomes |
