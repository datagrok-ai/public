# Vision

**A customizable DMTA workflow skeleton with embedded AIDD — built on Datagrok.**

This plugin is not a finished application. It is a *skeleton*: an opinionated but open framework on which a customer builds their own Design–Make–Test–Analyze process. Every step, computation, attribute, notification, and approval is meant to be extended, replaced, or removed. Out-of-the-box components get a team productive on day one; the "bring-your-own" seams let that same team grow the system into their proprietary discovery engine without leaving the platform or waiting on us.

---

## The problem we're solving

Drug discovery runs on the DMTA cycle, and the cycle is now a chain of AI models. In a modern loop, generative design proposes candidates, property and ADMET predictors and docking/co-folding score them, retrosynthesis tools judge whether they can be made, and active-learning selects the next batch — all under multi-parameter optimization that balances potency, selectivity, ADMET, and synthesizability at once. The models are improving monthly and are increasingly open (REINVENT4, the Boltz family, DiffDock, ADMET predictors, retrosynthesis engines).

Yet the software that orchestrates this loop has not kept up:

- **Incumbent suites are locked or heavy.** The strongest DMTA cockpits are best inside a single vendor's ecosystem and physics suite, or are desktop/server-rooted workflow engines that are powerful but slow to deploy and hard to extend. Enterprise rollouts routinely take 12–18 months.
- **AI-native platforms are closed pipelines.** Most sell an end-to-end pipeline or a discovery partnership, not a framework a customer can shape to their own science.
- **Assemble-it-yourself has no floor.** The open-model ecosystem is a genuine tailwind, but stitching REINVENT, a co-folding model, an ADMET predictor, and a retrosynthesis engine into a governed, auditable, human-in-the-loop process is a large in-house build with no UI, no governance, and no DMTA structure.

The gap is an **open, browser-native, composable orchestration layer** that is vendor-neutral about which models run inside it. That is the space this plugin targets.

## The vision

A team should be able to stand up a working DMTA process in days, then spend the following months making it *theirs* — not fighting the tool, not filing feature requests, not exporting to a notebook and losing the thread.

Concretely, the skeleton lets a customer:

- **Model their own cycle.** Define the stages of their DMTA process and the transitions between them, including custom stages that don't exist anywhere else.
- **Extend every step.** Add computations, attributes, and views; wire in their own scoring functions, filters, and triage logic.
- **Own the human-in-the-loop.** Configure approvals, review gates, notifications, and audit trails so the process fits how their chemists and biologists actually make decisions.
- **Choose their models.** Use Datagrok's built-in components, or plug in proprietary and third-party models — with no expectation that they abandon the tools they already trust.

The skeleton supplies structure and defaults. The customer supplies science and judgment. Neither is locked in.

## Core principles

**1. Skeleton, not cage.** Everything shipped is a starting point. If a customer replaces every default we provide and keeps only the orchestration spine, the plugin has succeeded, not failed.

**2. Bring-your-own-model is the headline, not a footnote.** Out-of-the-box components exist so nobody starts from zero. Open seams exist so nobody hits a ceiling. A proprietary in-house predictor and a public MIT-licensed model should be equally first-class citizens in the loop.

**3. Browser-native and zero-install.** The loop runs where the scientist already is — in the browser, on real datasets, interactively. No client install, no per-seat desktop provisioning, fast to deploy.

**4. One platform, every modality.** Small molecules, peptides, oligonucleotides, antibodies, and sequences belong in the same framework, not in separate chemistry-only and biology-only silos.

**5. Governance is part of the workflow, not bolted on.** Approvals, notifications, audit trails, and access control are core skeleton features because regulated, collaborative discovery demands traceability — and because human triage of generative output is where the loop earns its value.

**6. Composition over lock-in.** We orchestrate the ecosystem's best components rather than trying to out-build all of them. The moat is the workflow, the governance, and the openness — not any single model.

## What the skeleton provides

The framework ships a reference DMTA process — Design → Make → Test → Analyze — with example steps a team can run immediately and then fork:

- **Design.** Generative proposal, enumeration, and molecule optimization; property, ADMET, and docking/affinity scoring; retrosynthesis and synthesizability assessment; medicinal-chemist triage views.
- **Make.** Synthesis planning and route candidates; registration and tracking hooks.
- **Test.** Assay-result capture and linkage back to designs; data-quality and curation steps.
- **Analyze.** SAR, activity cliffs, matched molecular pairs, model performance, and batch selection for the next iteration.

Cross-cutting, every stage exposes the same extension points:

| Extension point | What a customer does with it |
|---|---|
| **Steps** | Add, reorder, or replace stages; define custom transitions and their conditions. |
| **Computations** | Attach scoring, filtering, and derived-property functions; swap defaults for proprietary logic. |
| **Attributes** | Extend the data model with custom fields, metadata, and provenance. |
| **Notifications** | Trigger alerts on state changes, thresholds, and hand-offs. |
| **Approvals & gates** | Insert human-in-the-loop review and sign-off wherever the process requires it. |
| **Views & UI** | Auto-generate interfaces from functions; customize dashboards and triage surfaces. |

## What Datagrok already brings

The skeleton is thin because the platform underneath is not. Datagrok is a browser-native scientific analytics platform (a custom in-memory columnar engine that stays interactive on very large datasets), with an extensible plugin architecture and a mature cheminformatics stack the DMTA loop draws on directly:

- **Cheminformatics (RDKit-WASM, client-side):** substructure and similarity search, descriptors, R-group analysis, scaffold trees, MCS, matched molecular pairs, activity cliffs, notation conversion, curation.
- **Predictive & generative modeling:** QSAR/QSPR and ADMET modeling; generative chemistry (REINVENT4); no-code ML plus scripting in Python/R/Julia and custom containers with GPU support.
- **Structure-based design:** molecular docking (AutoDock) with 3D visualization; retrosynthesis and virtual synthesis.
- **Multi-modality tooling:** peptide SAR/MSA and clustering, oligonucleotide format handling and modification support, sequence and biologics workflows.
- **Data & governance:** many data connectors, query builders, fine-grained access control, and one-click plugin deployment.

The DMTA plugin's job is to *orchestrate* these into a governed, iterative loop — not to reinvent them.

## Bring-your-own-model

The design phase is a moving target: the best generative, co-folding, affinity, ADMET, and retrosynthesis models change from quarter to quarter, and much of the frontier is now open and permissively licensed. A framework that hard-codes today's models is obsolete by next year.

So the skeleton treats models as pluggable. A customer can:

- Use built-in components (RDKit, AutoDock, REINVENT4, ADMET and property models) as sensible defaults.
- Register **proprietary in-house models** through scripting and a model-registry pattern, keeping them entirely private.
- Wire in **third-party and open models** via scripts, REST endpoints, containers, and model-tracking integrations.

**Roadmap.** First-class connectors for the modern open AI stack are a priority — generative chemistry, retrosynthesis, structure/affinity prediction, and docking — including the Boltz family of co-folding and binder-design models as they mature. (Datagrok integrates AutoDock and REINVENT4 today; broader co-folding and next-generation binder-design integrations are on the roadmap and should be treated as such until shipped.) The selection is demand-driven: we prioritize the models that show up most in target customers' actual processes.

## Multi-modality

Discovery is no longer small-molecule-only, and neither is this skeleton. The same DMTA framework, extension points, and governance apply across small molecules, peptides, oligonucleotides/siRNA, and antibodies/biologics. A team working across modalities gets one process to learn, one place to govern, and one loop to close — rather than a chemistry tool and a biology tool that never meet.

## The AIDD design loop

The framework is built to close the loop, not just to run steps once:

- **Multi-parameter optimization (MPO)** balances competing objectives — potency, selectivity, ADMET, synthesizability — as first-class scoring, not afterthoughts.
- **Active learning / closed loop:** batch selection, surrogate-model retraining hooks, and next-iteration feedback let the process improve as data accumulates, mirroring the lab cycle.
- **Human-in-the-loop by design:** generative output is triaged, reviewed, and approved by scientists at configurable gates, because models are unreliable outside their applicability domain and because accountability matters.

These orchestration primitives — batch selection, retraining hooks, approval gates, notifications, audit trails — are the defensible core. Pure model providers don't have them; pure analytics tools don't have them.

## How we're different

| | This plugin (on Datagrok) | Ecosystem-locked suites | AI-native pipelines |
|---|---|---|---|
| Architecture | Browser-native, zero-install | Often desktop/server-rooted | Varies; often closed SaaS |
| Extensibility | Customer extends every step | Extensible within the vendor's world | Largely fixed pipeline |
| Models | Bring-your-own, vendor-neutral | Tied to the vendor's models | The vendor's models |
| Modalities | Unified across all | Often chemistry- or biology-first | Often single-focus |
| Time to value | Days to deploy, fork, extend | Months-long rollouts | Onboarding or partnership |

The differentiated bet, stated plainly: **the open orchestration layer for the AIDD-era DMTA cycle** — browser-native, customizable to the customer's own process, neutral about which models run inside it, and unified across modalities.

We deliberately concede depth in specialized physics (e.g., free-energy perturbation) to dedicated suites and position to *integrate* their outputs rather than replicate them. The value we add is the loop, the governance, and the openness around those components.

## Non-goals

- **Not a physics engine.** We orchestrate docking, co-folding, and free-energy outputs; we don't aim to out-compute dedicated simulation suites.
- **Not a single-vendor model store.** We won't privilege our own models over a customer's or the open ecosystem's.
- **Not a rigid application.** If the skeleton can't be reshaped to a customer's process, that's a bug in the vision.
- **Not a lab-automation controller.** The "Make" stage plans and tracks; it doesn't drive instruments (though it can integrate with systems that do).

## Direction

The phasing below is directional, not a commitment:

1. **Reference skeleton.** Ship the Design→Make→Test→Analyze spine with runnable example steps and the full set of extension points, plus templates a team forks on day one.
2. **Open-model connectors.** First-class integrations for the priority generative, ADMET, docking, retrosynthesis, and co-folding models, driven by target-customer demand.
3. **Closed-loop orchestration.** Harden active-learning, MPO, batch selection, and human-in-the-loop governance as core primitives.
4. **Multi-modality depth.** Extend reference workflows and model connectors across peptides, oligos, and biologics.

**How we'll know it's working:** if early adopters heavily customize the skeleton, we double down on the framework and SDK; if most run it as-is, we invest in turnkey, verticalized variants. Either signal is a win — it tells us which product we're actually building.

---

*This document describes intent and direction. Capabilities marked as roadmap are not yet shipped; built-in capabilities reflect the current Datagrok platform. Model and modality priorities are demand-driven and expected to evolve.*
