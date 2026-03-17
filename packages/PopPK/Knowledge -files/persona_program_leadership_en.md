# Persona: Program Leadership (PM)

## Persona Card

| Attribute | Details |
|-----------|---------|
| **Role** | VP Clinical Development / Program Lead / Head of Clinical Pharmacology |
| **Department** | Clinical Development, R&D Leadership |
| **Organization** | Big Pharma, Mid-size Biotech |
| **Team size** | Oversees 20–100+ across Clinical Pharmacology, Pharmacometrics, and Medical functions |
| **Education** | MD, PharmD, or PhD; often with MBA or management training |
| **Experience** | 15–25+ years in drug development; 5+ in leadership |
| **Reports to** | Chief Medical Officer, SVP R&D, or CEO (Biotech) |
| **Budget authority** | Approves $500K–$5M+ annual software/services budgets; signs CRO contracts |

### Profile

The program leader makes high-stakes decisions that depend on modeling and simulation results: which dose to advance into Phase 3, whether to conduct a dedicated renal impairment study or rely on PopPK extrapolation, how to respond to an FDA Complete Response Letter. They do not build or run models — they consume results in the form of forest plots, simulation dashboards, and executive summaries. Their frustration is the turnaround time and opacity of the modeling process: asking a question today and receiving an answer next week, with no ability to interactively explore the implications.

### Goals

- Make timely, data-driven decisions on dose selection and clinical development strategy
- Minimize time and cost of the development program while maximizing probability of regulatory success
- Communicate the value of M&S to the executive committee and board
- Ensure compliance with regulatory expectations (ICH M15, FDA MIDD)

### Current tools

PowerPoint, PDF reports, email, dashboards (Tableau/Spotfire — rarely for PMx), clinical trial management systems

### Pain points

- **Week-long turnaround**: "What if we dose 200 mg BID instead of 300 mg QD?" takes a week to answer because PMx must re-run simulations
- **Opacity**: cannot see how the model works, cannot challenge assumptions, must trust the PMx black box
- **Static deliverables**: receives PowerPoint slides with frozen figures — cannot drill down, change parameters, or explore scenarios
- **Cost of indecision**: every week of delay in dose selection costs $500K–$2M in development timeline; poor dose choice can kill a $1B program
- **Board communication**: struggles to explain model-based evidence to non-technical executives and investors

### What they need from PopPK Studio

- Interactive dashboard: select dose, population, scenario → see predicted outcomes in real time
- Scenario comparison view: side-by-side dosing regimens with probability of target attainment
- Executive summary auto-generation with key findings highlighted
- Audit trail for decision documentation (what was decided, based on what evidence, when)
- Access from any device (browser-native = no IT installation requests)

### Quote

> "I don't care about the OFV or the omega matrix. I care about one thing: if I pick 200 mg for Phase 3, what's the probability that we hit the efficacy target in the intent-to-treat population? Show me that — live, now, in this meeting."

---

## User Journey Map

### Scenario: Dose selection decision for Phase 3 go/no-go at portfolio review

#### Stage 1 — Receive M&S Results Package (Day 1)

| Aspect | Details |
|--------|---------|
| **Task** | Review modeling and simulation results from the PMx team |
| **Actions** | Receive final PopPK/PD report → review executive summary → examine dose-response simulations → identify questions for the team |
| **Current experience** | 40-page PDF report arrives by email; skims executive summary; flags questions for follow-up meeting next week |
| **PopPK Studio experience** | Notification that analysis is complete → open interactive dashboard → browse key results → drill into specific figures → bookmark questions for PMx |
| **Emotion** | Current: overwhelmed (dense report) → PopPK Studio: oriented (interactive overview) |
| **Key metric** | Time to grasp key findings |

#### Stage 2 — Explore Scenarios (Day 1–3)

| Aspect | Details |
|--------|---------|
| **Task** | Evaluate alternative dosing strategies before the governance meeting |
| **Actions** | Compare 150 mg BID vs. 200 mg QD vs. 300 mg QD → check predicted safety margins → evaluate special populations (elderly, renal impairment) → assess dose-response plateau |
| **Current experience** | Send email to PMx: "Can you also simulate 150 mg BID?" → wait 3–5 days → receive updated slides |
| **PopPK Studio experience** | Open scenario comparison dashboard → toggle between regimens → adjust population filters (age > 65, eGFR < 60) → see predicted Cmax distributions and target attainment → save preferred scenario |
| **Emotion** | Current: frustration (waiting, can't explore) → PopPK Studio: empowerment (exploring in real time) |
| **Key metric** | Number of scenarios explored before decision; time from question to answer |

#### Stage 3 — Governance Meeting (Day 5)

| Aspect | Details |
|--------|---------|
| **Task** | Present dose recommendation to the clinical development governance committee |
| **Actions** | Present M&S-supported dose rationale → answer committee questions → document decision and rationale |
| **Current experience** | Present static slides; committee asks "what about lower dose in elderly?" → "We'll get back to you after the next modeling cycle" → decision delayed 2 weeks |
| **PopPK Studio experience** | Present from live dashboard → committee member asks about elderly → filter population on screen → show predicted exposure → decision made in-meeting |
| **Emotion** | Current: frustration (can't answer, delays decision) → PopPK Studio: decisive (real-time evidence) |
| **Key metric** | Decisions made per meeting; delay days avoided |

#### Stage 4 — Document Decision & Communicate (Day 5–7)

| Aspect | Details |
|--------|---------|
| **Task** | Document the dose selection decision with supporting evidence for regulatory trail and internal communication |
| **Actions** | Capture decision rationale → link to supporting analyses → communicate to broader team → archive for regulatory reference |
| **Current experience** | Write meeting minutes in Word; attach PDF report; email to stakeholders |
| **PopPK Studio experience** | One-click decision capture with linked model version, scenarios explored, and committee notes → auto-generate executive brief → share dashboard link with broader team → permanent auditable record |
| **Emotion** | Current: administrative burden → PopPK Studio: seamless documentation |
| **Key metric** | Decision documentation time; traceability quality |
