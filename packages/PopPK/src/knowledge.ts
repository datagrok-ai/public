export const knowledgeFiles: Record<string, string> = {
  'PopPK_Studio_Overview.md': `# PopPK Studio

## Overview

PopPK Studio (part of PharmaFlow) — a browser-native population PK/PD modeling module for the Datagrok platform. Implements the full pharmacometrics pipeline: NMTRAN data handling, compartmental models (analytical + ODE), NLME estimation (MAP → SAEM → FOCE-I), diagnostics (GOF, VPC, bootstrap), covariate modeling, Monte Carlo simulation, and regulatory-grade reporting. Serves four pharma personas:
- **Pharmacometricians** — code editor, DSL, full estimation control
- **Clinical pharmacologists** — GUI wizard, pre-built model library, self-service simulations
- **Biostatisticians** — exposure-response integration, MMRM bridge
- **Program leadership** — interactive dashboards, scenario comparison, executive reporting

No browser-native competitor exists in this space.

---

## Existing software

| Platform | Owner | Est. methods | Key weakness for us |
|----------|-------|-------------|-------------------|
| [NONMEM 7.6](https://www.iconplc.com/solutions/technologies/nonmem) | ICON plc | FO, FOCE-I, SAEM, NUTS | No GUI, no visualization, Fortran syntax, fragmented 4–6 tool workflow |
| [Monolix 2024R1](https://www.simulations-plus.com/software/monolix/) | Simulations Plus | SAEM only | No FOCE (complicates NONMEM comparison), desktop-only, no collaboration |
| [Phoenix NLME 8.7](https://www.certara.com/software/phoenix-nlme/) | Certara | FOCE, QRPEM | No SAEM, no Bayesian, Windows desktop only |
| [nlmixr2 4.0](https://nlmixr2.org/) | Open-source (R) | FOCE-I, SAEM | No GUI, steep learning curve, R-only |
| [Pumas 2.8](https://pumas.ai/) | PumasAI (Julia) | FOCE, SAEM, NUTS | Julia unfamiliar to most PMx, small community, opaque pricing |

**Gap**: zero browser-native platforms → zero unified workflows → 6–8 week analysis cycles.

→ *Full competitive analysis: [Pharmacometrics Software Competitive Landscape](software.md)*

---

## Target users

| Persona | Role | What they do | Primary pain |
|---------|------|-------------|-------------|
| **PMx** — [Pharmacometrician](persona_pharmacometrician_en.md) | Builds PopPK/PD models end-to-end | Writes control streams, runs estimation, VPC, bootstrap, regulatory reports | 4–6 disconnected tools; hours debugging Fortran syntax; no collaboration |
| **CP** — [Clinical Pharmacologist](persona_clinical_pharmacologist_en.md) | Translates models into clinical decisions | Dose justification, special populations, drug label sections | Every "what if?" takes days; black-box dependency on PMx |
| **BS** — [Biostatistician](persona_biostatistician_en.md) | Designs trials, runs E-R analyses | MMRM, exposure-response, sample size calculations | 2 days merging PMx exposure metrics with efficacy datasets |
| **PM** — [Program Leadership](persona_program_leadership_en.md) | Makes dose selection decisions | Reviews M&S results, presents to governance and FDA | Static slides; can't answer ad-hoc questions; $500K–$2M/week delay cost |

→ *Full persona cards and user journey maps:*
- *[Pharmacometrician](persona_pharmacometrician_en.md)*
- *[Clinical Pharmacologist](persona_clinical_pharmacologist_en.md)*
- *[Biostatistician](persona_biostatistician_en.md)*
- *[Program Leadership](persona_program_leadership_en.md)*

---

## Pipeline

![PopPK/PD Pipeline with Personas](popPKPD_pipeline_personas_en.svg)

**Input**: dosing records, PK blood samples, PD biomarker measurements, patient covariates (weight, age, renal function, genotype) — typically in NMTRAN format (ID, TIME, DV, AMT, EVID, CMT).

**Output**: dose recommendations for target populations, Clinical Pharmacology section of drug label, design of next clinical trial phase, regulatory-grade PopPK report (Module 2.7.2 eCTD) for FDA/EMA submission.

---

## Minimal algorithm set

28 algorithms required for MVP and first population estimation capability.

| Category | Count | Already in Datagrok | Key algorithms |
|----------|-------|-------------------|----------------|
| ODE Solvers | 5 | **5 of 5** | Rosenbrock, RK3/4/DP, Adams AB4/5, **LSODA**, **CVODE** |
| Analytical PK | 2 | 0 | Matrix exponential, superposition for multi-dose |
| Optimization | 4 | **1 of 4** | **Nelder-Mead**, L-BFGS, L-BFGS-B, Levenberg-Marquardt |
| NLME Estimation | 3 | 0 | MAP/EBE, SAEM, Metropolis-Hastings sampler |
| Linear Algebra | 4 | 0 | Cholesky, LU, eigendecomposition, log-determinant |
| Statistics | 7 | **3 of 7** | MVN sampling/density, seedable PRNG, **Monte Carlo**, **Sobol**, CWRES, shrinkage |
| Diagnostics | 3 | ~1.5 | GOF plots, QQ-plot, individual fits (partial viewer support) |
| **Total** | **28** | **~10.5** | **~19% coverage** |

Strongest coverage: **ODE solving is 100% complete** (including LSODA and SUNDIALS CVODE — the same solvers used by NONMEM and Monolix). Primary gaps: NLME estimation (0/3), linear algebra (0/4), optimization (1/4).

→ *Full algorithm catalog with references and libraries: [Algorithms and Numerical Methods](PopPK_Studio_Algorithms_and_Methods.md)*

→ *Terminology: [Glossary](glossary_en.md)*`,


  'software.md': `# Pharmacometrics software: a competitive landscape in transition

**The population PK/PD modeling software market — valued at approximately $613–680 million in 2024 and growing at 10–13% CAGR — is experiencing its most significant disruption in decades.** NONMEM still commands roughly 85% of FDA PopPK submissions after 45 years of dominance, but cloud-native platforms, open-source alternatives, and AI-integrated tools are rapidly fragmenting the landscape. No comprehensive browser-native pharmacometrics platform exists today, creating a clear whitespace opportunity for a Datagrok-based solution. This analysis evaluates every major competitor across technology, regulatory positioning, pricing, and strategic trajectory to inform product strategy.

---

## The established giants: NONMEM and Certara Phoenix

### NONMEM (ICON plc) — the regulatory gold standard

NONMEM was created in the late 1970s by **Stuart Beal and Lewis Sheiner at UCSF** and remains the single most important tool in population pharmacokinetics. Ownership passed from UCSF to GloboMax (~2000), then to **ICON plc** (NASDAQ: ICLR) when ICON acquired GloboMax in 2003. Today, ICON is a $8B+ CRO with 41,000 employees; NONMEM is a tiny but strategically prestigious product maintained primarily by a single developer, Robert Bauer.

**Technology.** Written in **Fortran 95**, NONMEM is a command-line batch program with no native GUI. Users write NM-TRAN control streams — a Fortran-derived domain-specific language — which are compiled on-the-fly. Version **7.6.0 (2024)** is current, running on Windows and Linux. It offers the **most comprehensive estimation algorithm suite** in the field: FO, FOCE, FOCE-I, Laplacian, ITS, IMP, IMPMAP, SAEM, full Bayesian MCMC, and NUTS (Hamiltonian Monte Carlo). Its ODE solver library is equally deep: LSODA (ADVAN13) for automatic stiff/non-stiff switching, CVODES from SUNDIALS (ADVAN14), and uniquely, **delay differential equation solvers** (ADVAN16–18 using RADAR5). Version 7.6 added distributed delay convolution algorithms and a **$DESIGN block** for optimal experimental design.

**Regulatory position.** NONMEM is referenced in an estimated **85% of FDA-submitted PopPK analyses**. The FDA's own Division of Pharmacometrics purchases NONMEM licenses. Every regulatory reviewer worldwide is familiar with NONMEM output formats. This deep institutional entrenchment is NONMEM's most formidable competitive moat — one that has resisted challenge for three decades.

**Ecosystem.** NONMEM's weaknesses in usability are compensated by a rich third-party ecosystem: **PsN** (Perl-speaks-NONMEM, v5.7.0, Uppsala University) provides bootstrap, VPC, stepwise covariate modeling, and cluster computing; **Xpose** delivers ggplot2-based diagnostics in R; **Pirana** (now Certara-owned) offers a GUI workbench; and ICON's own **PDx-Pop** provides execution management.

**Pricing.** Academic licenses run approximately **$990 for 3 CPUs** (per 2015 reference). Industry pricing is not publicly disclosed but estimated at **$5,000–15,000+ per CPU/year**. There is no free tier.

**Critical weaknesses.** The steep learning curve (Fortran-derived syntax, command-line only), lack of built-in visualization, absence of a model library, requirement for Fortran compiler installation, no native cloud solution, and dependence on a single core developer represent significant vulnerabilities. NONMEM's market position depends almost entirely on regulatory inertia rather than technological superiority.

### Certara Phoenix WinNonlin / NLME — the enterprise platform

Certara was formed in 2008 through the Vector Capital-orchestrated merger of **Pharsight Corporation** (creator of WinNonlin) and Tripos International. The company executed an aggressive acquisition strategy — Simcyp (2012, PBPK), Synchrogenix (2014, regulatory writing), Pinnacle 21 (2021, CDISC compliance, ~$310M), Chemaxon (2024, cheminformatics) — before going public in **December 2020** (NASDAQ: CERT), raising $670M at $23/share. FY2024 revenue reached **$385 million** (software $156M, services $229M), with Q3 2025 showing 22% YoY software revenue growth to $43.8M.

**Technology.** The NLME engine is written in **C/C++** with a native Windows GUI (Phoenix Workbench). Phoenix Modeling Language (PML) provides text-based model specification alongside a **graphical model editor** and extensive pre-built model library. The platform supports three model entry methods — library selection, drag-and-drop building, or PML scripting — dramatically lowering the barrier to entry. Phoenix **8.6 (June 2025)** introduced **automatic differentiation (ADPO)** for FOCE ELS, Laplacian, and AGQ algorithms, delivering **20–50% computation time reductions**. Version 8.7 is the latest release.

**Estimation methods** include QRPEM (a unique quasi-random parametric EM method using Sobol sequences), FOCE (Lindstrom-Bates and ELS variants), Laplacian, Adaptive Gaussian Quadrature, and naïve pooled estimation. Critically, Phoenix **lacks SAEM and Bayesian MCMC/NUTS** — a significant gap versus NONMEM and newer competitors.

**Market position.** Phoenix WinNonlin is the **undisputed gold standard for NCA analysis**, used by 6,000+ researchers at 1,500+ organizations. **75+ of top 100 pharmaceutical companies** use Phoenix. For PopPK specifically, NLME is growing but remains secondary to NONMEM. Certara's broader ecosystem — Simcyp (PBPK, qualified by EMA in August 2025), Pinnacle 21 (CDISC), CoAuthor (AI regulatory writing), CertaraIQ (AI-powered QSP) — creates a comprehensive biosimulation platform unmatched in breadth.

**Cloud strategy.** Phoenix Hosted provides cloud deployment with SSO via Certara Cloud. The Integral platform offers an API-enabled GxP data repository. RsNLME packages enable R-based execution. The cloud strategy is accelerating but remains largely a hosted desktop experience rather than a browser-native application.

---

## The SAEM challenger: Monolix and Simulations Plus

MonolixSuite emerged from a **2003 research collaboration at Inria** (France) led by mathematician Marc Lavielle, who pioneered the practical implementation of the SAEM algorithm. **Lixoft** was spun off from Inria in 2011 by Jérôme Kalifa and Lavielle, commercializing Monolix. **Simulations Plus** (NASDAQ: SLP) acquired Lixoft in **April 2020 for up to $16.5 million** — a remarkably modest price for what is now a cornerstone of their platform.

Simulations Plus reported **FY2025 revenue of $79.2 million** (+13% YoY), with a portfolio spanning GastroPlus (PBPK), ADMET Predictor (AI cheminformatics), DILIsym (QST), and MonolixSuite. The June 2024 acquisition of Pro-ficiency for ~$100M doubled their addressable market.

**Core technology.** MonolixSuite 2024R1 (the latest release as of March 2026) is a **C++ engine** with Qt-based GUI running on Windows, Linux, and macOS. Models are specified in **Mlxtran**, a declarative language that is considerably more readable than NONMEM's NM-TRAN. The platform includes a graphical model editor and an extraordinary **library of 30,000+ pre-built models** covering PK, PK/PD, TMDD, tumor growth, count data, and time-to-event models.

**Estimation.** Monolix offers **only SAEM** — no FOCE. This is both its greatest strength and its most debated limitation. SAEM performs maximum likelihood estimation without linearization approximation (unlike FOCE), using a two-phase approach with burn-in, exploratory (with simulated annealing), and smoothing phases. Auto-stop convergence detection and multi-start capabilities are built in. The **COSSAC method** for automated covariate search requires **2–20x fewer runs** than traditional stepwise covariate modeling while producing equivalent results.

**ODE solvers** include Adams (non-stiff), BDF (stiff), CVODE v5.5, a DDE solver, and even a PDE solver using finite element methods — a capability unique among pharmacometrics tools.

**Regulatory acceptance** is growing but NONMEM's inertia remains the primary barrier. All MonolixSuite file types are accepted by the FDA, and regulators have published research using Monolix. **Free academic licenses** have driven strong university adoption.

**Key weaknesses** include no FOCE support (complicating direct NONMEM comparisons for regulatory reviewers), a smaller community and ecosystem than NONMEM, no native Python API, and no cloud-native offering yet (though Simulations Plus's S+ Cloud strategy targets FY2026 integration).

---

## The open-source ecosystem: nlmixr2, Stan/Torsten, and supporting tools

### nlmixr2 — the R-native NONMEM alternative

Led by **Matt Fidler** with contributors from Novartis, LAP&P, Occams, and multiple pharma companies, nlmixr2 (v5.0.0, November 2025, GPL-3 license) has matured into a serious NONMEM competitor. In November 2025, it became an **R Consortium ISC Working Group**, establishing a sustainable multi-company funding model.

The architecture is modular: **nlmixr2** (umbrella), **nlmixr2est** (estimation), **rxode2** (compiled C ODE simulation engine), plus 15+ companion packages. Models use R-like syntax compiled to C for execution. Estimation methods include **FOCEI, FOCE, SAEM, and Laplace/AGQ** (added 2025). ODE solving via LSODA and VODE with OpenMP parallelization. Published concordance studies show **population bias of 0–5.3%** versus NONMEM across multiple models.

The **babelmixr2** package enables bidirectional model conversion between nlmixr2, NONMEM, and Monolix — a unique interoperability advantage. **nonmem2rx** imports NONMEM control streams with automatic validation. Version 4.0 (July 2025) switched to WNL-type linear compartment solutions for improved NONMEM compatibility.

nlmixr2's primary limitation remains **regulatory acceptance track record**. While FDA guidance is software-agnostic, institutional comfort with NONMEM output persists. The strategy of enabling NONMEM execution through babelmixr2 provides a practical bridge.

### Stan/Torsten — Bayesian gold standard

**Stan** (Columbia University, NumFOCUS, BSD-3 license) provides state-of-the-art Bayesian inference via NUTS/HMC. **Torsten** (Metrum Research Group) extends Stan with pharmacometric-specific functions — compartmental model solvers, event handling following NONMEM data conventions, and steady-state computation. Funded by the Office of Naval Research and the Bill & Melinda Gates Foundation.

Stan/Torsten offers complete posterior uncertainty quantification and flexible hierarchical modeling that no frequentist tool can match. However, it remains computationally expensive, has a steep learning curve (requiring both Stan programming and Bayesian statistics expertise), and Torsten is still labeled as a prototype. Pfizer published a Bayesian PopPK application using Stan/Torsten in 2025, signaling growing industry interest.

### Supporting tools in the R ecosystem

The broader R pharmacometrics ecosystem includes **mrgsolve** (Metrum Research Group, fast C++ ODE simulation for QSP/PBPK), **Pmetrics** (USC, non-parametric population modeling via NPAG), **PopED** (Uppsala University, FIM-based optimal experimental design), **PFIM** (INSERM/University of Paris, design optimization), **vpc** (visual predictive checks), **xpose** (diagnostics), and **PKNCA** (NCA). Together, these packages create a comprehensive open-source toolkit — though one that requires substantial R expertise and manual integration.

### ADAPT (USC)

Developed since 1985 by **David D'Argenio** at USC's Biomedical Simulations Resource, ADAPT 5 (2009) provides individual and population estimation methods (MLEM, ITS, MAP Bayesian) and **optimal experimental design** (D-optimal, C-optimal). Written in Fortran, Windows-only, free under a Release Agreement but not truly open-source. It remains relevant primarily for optimal design and academic research but is showing its age.

---

## The next-generation disruptor: Pumas and PumasAI

**PumasAI** was founded in 2019 by **Joga Gobburu and Vijay Ivaturi** (University of Maryland) with software development led by **Chris Rackauckas**, as a joint venture with JuliaHub. Vijay Ivaturi became CEO in January 2024 and currently serves as President of ISoP.

**Why Pumas matters.** Built entirely in **Julia**, Pumas eliminates the "two-language problem" — delivering an R-like interface at Fortran-level speed. Benchmarks show a **median 81x acceleration** for NLME fitting versus competitors and **4x improvement** in ODE solving. On a 15-stiff-ODE PBPK model, Pumas demonstrated **2–4x advantage over CVODE** (used by NONMEM and Monolix). These performance gains are not incremental — they are transformational for large-scale simulations and PBPK modeling.

**Version 2.8 (February 2026)** supports FO, FOCE (with automatic interaction detection), LaplaceI, SAEM (via @emmodel), BayesMCMC (NUTS), and newly added MCEM and VEM — making it the most complete estimation suite alongside NONMEM. ODE solving leverages Julia's DifferentialEquations.jl, the most comprehensive solver ecosystem in any language, with automatic stiffness detection.

**DeepPumas** integrates Scientific Machine Learning (Neural ODEs, Universal Differential Equations, Mixed-Effect Neural Networks) directly into the pharmacometric workflow — a capability no competitor offers at this level. **PumasAide**, an AI coding agent introduced in v2.8, generates, executes, and refines models through natural language interaction.

**Regulatory traction** is building: **26 FDA submissions** reference Pumas/JuliaHub. The platform runs on a **SOC-2 and CFR Part 11 compliant** JuliaHub cloud environment. Moderna has publicly stated that "Pumas has emerged as our go-to tool for most of our analyses."

**Weaknesses** include a smaller established user base, Julia's unfamiliarity to most pharmacometricians (who know R and Fortran), opaque pricing, and the SAEM/@emmodel separation that adds complexity. The JuliaHub dependency creates vendor concentration risk.

---

## Cloud infrastructure and emerging platforms

**Metworx** (Metrum Research Group, est. 2004) is the leading cloud-native pharmacometrics **Platform-as-a-Service**, providing pre-configured AWS environments with NONMEM, R, Stan, PsN, Torsten, and mrgsolve — all GxP-validated with IQ/OQ/PQ packages. Over 200 companies use Metworx. In September 2025, Metworx integrated **Aplos NCA** — the first truly browser-based NCA tool — into its ecosystem. Metworx enables "thousands of simulations in parallel, cutting a 3-month timeline down to just weeks" but is fundamentally cloud *infrastructure* around existing tools, not a browser-native application.

**MATLAB SimBiology** (MathWorks) provides pharmacometrics capabilities including NLME estimation (FO, FOCE, SAEM), NCA, graphical model building, and sensitivity analysis. Its strength lies in MATLAB ecosystem integration but it suffers from expensive licensing, a non-pharmacometrics-first design, and a small PMx-specific community.

**InsightRX** has achieved remarkable scale in clinical precision dosing — **1,000+ hospital deployments**, 4 million+ dose advices, 150+ validated PK/PD models — with a KLAS Spotlight 2025 score of 96.1. Its Vantage Platform extends into clinical trials. However, InsightRX applies existing models at the bedside rather than building new population models.

**Berkeley Madonna** (Version 10) serves as a visualization and teaching tool for ODE-based models with drag-and-drop building and interactive parameter sliders, but lacks NLME estimation entirely.

---

## Market dynamics and the pricing landscape

The pharmacometrics software market was estimated at **$613–680 million in 2024**, projected to reach **$1.4–2.1 billion by 2033** at 10–13% CAGR. North America holds ~41% market share, with Asia Pacific growing fastest at 13.5–16.5% CAGR.

| Platform | Academic pricing | Industry pricing (est.) | Model |
|---|---|---|---|
| **NONMEM** | ~$990/3 CPUs | $5,000–15,000+/CPU/yr | Per-CPU annual subscription |
| **Phoenix WinNonlin/NLME** | Free (teaching); discounted | $5,000–15,000+/seat/yr | Per-seat or hosted |
| **MonolixSuite** | Free (non-commercial) | Contact sales | Named/floating/server |
| **Pumas** | Free (academic/research) | Contact sales (SaaS) | Cloud SaaS |
| **nlmixr2** | Free | Free | Open-source (GPL-3) |
| **Stan/Torsten** | Free | Free | Open-source (BSD-3) |
| **ADAPT** | Free | Free | Release agreement |
| **Metworx** | N/A | Enterprise cloud pricing | PaaS |

A typical pharma team running NONMEM + Phoenix + Simcyp may spend **$50,000–200,000+ annually** on software licenses alone, before computing infrastructure and CRO consulting costs ($100K–500K per pharmacometric analysis).

---

## Five trends reshaping the competitive landscape

**Cloud migration is accelerating.** 65% of pharmaceutical software was cloud-deployed in 2024, growing at 17% CAGR. Metworx, JuliaHub, and Certara Cloud lead, but most "cloud" offerings are hosted desktop experiences, not browser-native applications.

**AI/ML integration is becoming table stakes.** DeepPumas (Neural ODEs), Certara's Darwin (ML model search), PumasAide (AI coding agent), and MonolixSuite 2024R1 (Neural ODE support) all signal that AI augmentation of mechanistic models is transitioning from research novelty to commercial feature.

**Regulatory harmonization is expanding the market.** ICH M15 (Model-Informed Drug Development) is creating global standards for MIDD, while the FDA's PDUFA VII continues the MIDD Paired Meeting Program through June 2027. MIDD yields an estimated **10 months cycle time savings and $5M per program annually**, driving adoption beyond traditional pharmacometrics groups.

**Open-source credibility is growing.** nlmixr2's R Consortium status, concordance studies demonstrating NONMEM equivalence, and the babelmixr2 interoperability bridge are eroding the "only NONMEM is accepted" paradigm. The FDA's 2022 PopPK guidance is explicitly software-agnostic.

**The usability gap remains massive.** A 2025 CPT:PSP paper documented significant UX deficiencies in pharmacometrics tools. Pharmacometricians still routinely use 4–6 disconnected tools per analysis. Interactive R Shiny apps (gPKPDviz, ModVizPop, e-Campsis) are filling some gaps but lack platform-level integration.

---

## Competitive vulnerability map: where each platform falls short

| Vulnerability | Affected platforms | Opportunity |
|---|---|---|
| No browser-native interface | NONMEM, Monolix, Phoenix (desktop), ADAPT | Browser-first platform |
| Fragmented multi-tool workflows | NONMEM ecosystem (PsN + Xpose + Pirana + R) | Unified workflow |
| No built-in visualization | NONMEM | Integrated interactive VPC/GOF/diagnostics |
| No SAEM or Bayesian | Phoenix NLME | Complete estimation suite |
| No FOCE | Monolix | Both FOCE and SAEM in one tool |
| Steep learning curve | NONMEM, Stan/Torsten, ADAPT | Guided model building, templates |
| No collaboration features | All desktop tools | Multi-user real-time collaboration |
| Poor data integration | All (manual dataset creation) | Direct CDISC/SDTM/ADaM connectors |
| No cloud-native architecture | NONMEM, Monolix, ADAPT | Cloud-first scalable compute |
| High cost for small biotechs | NONMEM, Phoenix, Monolix (industry) | Affordable browser-based access |
| Julia unfamiliarity | Pumas | R/Python-native interface |
| Reproducibility challenges | All (especially NONMEM) | Version-controlled, auditable analyses |

---

## The browser-based pharmacometrics gap: Datagrok's opportunity

**No comprehensive browser-native population PK/PD modeling platform exists today.** This is the single clearest gap in the competitive landscape. Current "cloud" offerings are either hosted desktop applications (Phoenix Hosted), cloud infrastructure around existing tools (Metworx), code-based cloud environments (JuliaHub/Pumas), or single-purpose browser tools (Aplos NCA, R Shiny apps). None provides an integrated, browser-native workflow from data import through estimation to interactive reporting.

Datagrok — a web-based data analytics platform with 40+ database connectors, exceptional in-browser performance (handling millions of data points), native R/Python/Julia integration, and established pharma partnerships (TetraScience) — is architecturally positioned to fill this gap. Its plugin architecture could integrate open-source estimation engines (nlmixr2/rxode2 for FOCE/SAEM, Stan/Torsten for Bayesian) while providing what no competitor offers: **unified data management, interactive visualization, real-time collaboration, and guided workflows in a zero-installation browser interface**.

The strategic playbook for a Datagrok pharmacometrics platform should prioritize: (1) browser-native NCA as the entry wedge (competing with Aplos, validating against Phoenix WinNonlin), (2) integrated PopPK modeling powered by nlmixr2/rxode2 backend engines with a graphical model builder comparable to Monolix's, (3) interactive diagnostics (VPC, GOF, forest plots) leveraging Datagrok's visualization strengths, (4) NONMEM interoperability via babelmixr2/nonmem2rx to ease migration, and (5) collaboration and reproducibility features (version control, audit trails, shared projects) targeting the regulatory compliance requirements that CFR Part 11 demands.

## Conclusion

The pharmacometrics software market is at an inflection point. NONMEM's dominance rests on regulatory inertia, not technology — a position that erodes with every FDA submission using an alternative tool. Pumas has demonstrated that **81x performance gains** are achievable with modern computing. The open-source ecosystem has reached a maturity where nlmixr2 produces statistically concordant results with NONMEM. Yet every platform — commercial and open-source alike — remains fundamentally a desktop or code-first tool deployed with varying degrees of cloud wrapping. The pharmacometrician's daily experience still involves juggling control streams, R scripts, diagnostic packages, and PowerPoint reports across disconnected applications. A browser-native platform that unifies this workflow while integrating proven open-source engines could capture the underserved segments (small biotech, academia scaling to industry, cross-functional teams) while establishing a new paradigm for how pharmacometric analyses are conducted, shared, and submitted to regulators.`,


  'persona_pharmacometrician_en.md': `# Persona: Pharmacometrician (PMx)

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
| **Key metric** | Report preparation time; number of QC errors caught |`,


  'persona_clinical_pharmacologist_en.md': `# Persona: Clinical Pharmacologist (CP)

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
| **Key metric** | Ability to answer ad-hoc questions; meeting outcomes |`,


  'persona_biostatistician_en.md': `# Persona: Biostatistician (BS)

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
| **Key metric** | Number of reconciliation cycles needed; cross-team alignment time |`,


  'persona_program_leadership_en.md': `# Persona: Program Leadership (PM)

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
| **Key metric** | Decision documentation time; traceability quality |`,


  'PopPK_Studio_Algorithms_and_Methods.md': `# PopPK Studio — Algorithms and Numerical Methods

## Datagrok current capabilities (DiffStudio / diff-grok)

Based on the DiffStudio documentation and the diff-grok npm package, Datagrok currently provides:

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

| Algorithm | Purpose | Datagrok | Key references |
|-----------|---------|----------|----------------|
| Rosenbrock-Wanner methods | Stiff ODE integration for compartmental PK/PD models | **Yes** — diff-grok | Rang & Angermann (2005); Hairer & Wanner (1996) |
| Runge-Kutta methods (RK3, RK4, Dormand-Prince) | Non-stiff ODE integration; adaptive step-size control | **Yes** — diff-grok | Dormand & Prince (1980); Bogacki & Shampine (1989) |
| Adams multistep methods (AB4, AB5) | Non-stiff predictor-corrector methods | **Yes** — diff-grok | Bashforth & Adams (1883) |
| LSODA | Adaptive solver switching between Adams and BDF automatically | **Yes** — diff-grok | Hindmarsh (1983); Petzold (1983) |
| CVODE (SUNDIALS) | Production-grade stiff/non-stiff solver; used by NONMEM and Monolix | **Yes** — diff-grok | Hindmarsh et al. (2005) |

### 1.2 Analytical PK Solutions

| Algorithm | Purpose | Datagrok |
|-----------|---------|----------|
| Matrix exponential for linear compartmental systems | 100–1000x faster than numerical ODE | **No** |
| Superposition principle for multiple dosing | Steady-state and multi-dose profiles | **No** |

### 1.3 Optimization Algorithms

| Algorithm | Purpose | Datagrok |
|-----------|---------|----------|
| Nelder-Mead simplex | Derivative-free optimization; MAP estimation | **Yes** — compute-utils |
| L-BFGS | Quasi-Newton gradient-based; outer loop of FOCE-I | **No** |
| L-BFGS-B (box-constrained) | Essential for constrained PK parameters (CL > 0, V > 0) | **No** |
| Levenberg-Marquardt | Nonlinear least-squares; individual curve fitting | **No** |

### 1.4 NLME Estimation Methods

| Algorithm | Purpose | Datagrok |
|-----------|---------|----------|
| MAP / EBE | Individual parameter estimation given population priors | **No** |
| SAEM | Population NLME estimation; no gradient required | **No** |
| Metropolis-Hastings sampler | MCMC sampling within SAEM E-step | **No** |

### 1.5 Linear Algebra

| Algorithm | Purpose | Datagrok |
|-----------|---------|----------|
| Cholesky decomposition | Decompose positive-definite matrices; MVN sampling | **No** |
| LU decomposition | Solving linear systems; matrix inversion | **No** |
| Eigendecomposition | Matrix exponential for analytical PK solutions | **No** |
| Log-determinant | Computing ln|Ci| in FOCE-I objective function | **No** |

### 1.6 Statistical Methods

| Algorithm | Purpose | Datagrok |
|-----------|---------|----------|
| Multivariate normal sampling | Generating random eta vectors in SAEM and simulation | **No** |
| Multivariate normal density | Evaluating p(eta|Omega) in MAP objective function | **No** |
| Seedable PRNG | Reproducible simulation and SAEM sampling | **No** |
| Monte Carlo sensitivity analysis | Explore function behavior at randomly sampled points | **Yes** — DiffStudio |
| Sobol variance decomposition | Decompose output variance into input contributions | **Yes** — DiffStudio |
| CWRES | Primary diagnostic metric for NLME model adequacy | **No** |
| Shrinkage calculation | Quantify information content of individual estimates | **No** |

### 1.7 Diagnostics and Visualization

| Algorithm | Purpose | Datagrok |
|-----------|---------|----------|
| GOF plots | DV vs PRED, DV vs IPRED, residuals vs time/predictions | **Partially** |
| QQ-plot | Assess normality of CWRES and ETAs | **Partially** |
| Individual fit plots | Observed vs. predicted per subject | **Partially** |

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

**Coverage: Datagrok provides ~19% of the required algorithms.**`,


  'glossary_en.md': `# PopPK Studio — Glossary

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
| WASM (WebAssembly) | Binary instruction format enabling near-native performance in web browsers |`,
};
