# Project Research Summary

**Project:** Datagrok Proteomics — v1.4 "Cross-Team Review" Milestone
**Domain:** Mass-spec proteomics analytics — read-only publishing, longitudinal sample-level SPC, screening-campaign cross-run comparison (Cytokinetics-driven)
**Researched:** 2026-06-06
**Confidence:** HIGH (stack + architecture verified against repo HEAD; features grounded in MSstatsQC/SProCoP/Spectronaut precedent + 2026-05-04 Cytokinetics meeting notes; pitfalls anchored in Phase 13 round-3 scar tissue)

## Executive Summary

v1.4 is an additive milestone on a v1.3-stable foundation. The three user stories (read-only publish to biologists, longitudinal SPC across weekly assays, and cross-run compound comparison for screening campaigns) are independently scoped, share one cross-cutting data primitive (the "run identity" / "campaign" model), and ship together as the close-out for the Cytokinetics engagement. **Zero new npm dependencies** — all three capabilities land on platform-native primitives (`DG.Project` + `grok.dapi.permissions` for publishing; ~390 LOC of in-package Shewhart/Nelson math for SPC; tag-based campaign model with AppData CSV persistence following the HitTriage precedent for cross-run comparison). The one cross-package coupling is `Chem` for SMILES rendering and (optionally) MCS substructure highlight, accessed via the existing `grok.functions.call('Chem:...')` / `DG.Func.find` pattern already proven in v1.3's Dendrogram integration — no `@datagrok/chem` npm dep.

All four researchers converge on a four-phase delivery (15→18) layered Publishing → SPC → Campaign Data Model → Cross-Run Comparison, with one material disagreement: FEATURES.md treats the "screening run" data model as a shared spine that should land before either Story 2 or Story 3 viewers, while ARCHITECTURE.md folds the data model into Phase 17 directly. **Resolution: keep the four-phase shape; Phase 17 leads with Task 0 = data model FIRST** (tags + SEMTYPE + AppData schema, before any viewer code). SPC math (Phase 16) is independent of the campaign model because SPC works on any single df with intensity columns; the campaign-index just consumes the per-run `proteomics.spc_status` tag at index time.

Key risks concentrate in publishing (stale-snapshot leak, permission inheritance bypass, tag-stripping on Project serialization, versioning ambiguity) and the campaign data model (inner-join silent loss, compound name drift, subscription leak from copy-pasting v1.2's enrichment cross-link pattern). The most material open question: **whether `proteomics.*` tags and `Proteomics-*` semTypes survive `DG.Project` save/reopen** — Phase 13 round-3 already proved layout-config tags are partially stripped (`look.filters[]` evidence at commit `e527d07ba1`). Publishing phase must lead with a round-trip test (`assertPublishedShape`) and a belt-and-braces "encode the critical tag value in a column as well" mitigation. Reviewer-group naming convention and compound column source are the two requirements-side decisions that gate Phase 15 and Phase 17 respectively.

## Key Findings

### Recommended Stack

v1.3 baseline unchanged: TypeScript + webpack + `datagrok-api ^1.25.0` + R scripts via `grok.functions.call` with TS fallbacks. v1.4 ships three new in-package modules on existing primitives. Single integration point with another Datagrok package is **Chem**, accessed via cross-package `grok.functions.call('Chem:...')` and `DG.Func.find` patterns — not as an npm dependency. See [STACK.md](./STACK.md).

**Core technologies (additions only):**
- **`DG.Project` + `grok.dapi.permissions.grant`**: read-only publishing — native, durable; permission API supports group-scoped view/edit; canonical 3-line shape verified in `packages/ApiSamples/scripts/dapi/projects.js` and `packages/Bio/src/tests/projects-tests.ts:26`.
- **DataFrame tags (`proteomics.campaign_id`, `proteomics.run_id`, `proteomics.compound`) + AppData CSV index**: campaign data model — matches v1.3 "everything is a `proteomics.*` tag" contract; HitTriage uses the same shape (`packages/HitTriage/src/app/hit-triage-app.ts:375`); avoids over-engineering vs `@datagrok-libraries/cruddy` or a custom Postgres schema.
- **In-package SPC (~390 LOC TS)**: Shewhart I-MR + Nelson rules 1–8 + per-run QC extractor + viewer + index store. `@datagrok-libraries/statistics` ships no SPC primitives (grep confirmed); no usable npm. Math is public-domain (Nelson 1984; Wheeler's d2=1.128, D4=3.267 for n=2) and deterministic enough for TS-only client path with optional R `qcc` fallback parked for v1.5.
- **`grok.functions.call('Chem:detectSmiles')` + `col.semType = 'Molecule'`**: compound structure rendering — Chem auto-registers the `'Molecule'` cell renderer; setting semType wires it up with zero coupling. Optional `Chem:convertMolNotation`, `Chem:mcsSearch` (Pattern B with `DG.Func.find` + arity guard + fallback) for canonicalization and substructure highlight.

### Expected Features

Three independently scoped user stories. Table stakes anchored in prior art (Spectronaut Direct / Scaffold Viewer for publishing; MSstatsQC / SProCoP for SPC; "Excel + PowerPoint" for cross-run comparison). Biggest competitive moat is Story 3 — no prior art ships cross-run, compound-aware proteome comparison. See [FEATURES.md](./FEATURES.md).

**Must have (table stakes):**
- **Publishing** — Trimmed published DataFrame (Protein ID, Gene, log2FC, p-value, adj.p-value, sig flag, direction); volcano renders on open with thresholds drawn; audit context panel reading `proteomics.*` tags; frozen/no-edit; target-keyed filing; reviewer-group view-only ACL; versioned/dated artifact name.
- **SPC** — Four standard metrics per run (median log2 intensity, % missing, control-rep Pearson correlation, protein count above threshold); Shewhart I-chart + MR-chart per metric; Nelson rules 1–4 default with selectable enable; run-level `proteomics.spc_status` tag; click-to-drill; configurable baseline (rolling N=10 default or fixed reference); phase-aware mode for lot/instrument changes.
- **Campaign comparison** — First-class screening-run + compound entities (tags + SMILES); Annotate Screening Campaign dialog (sibling to v1.3 Annotate Experiment); side-by-side volcano with linked selection and shared axes; diff table (sig in A only / B only / both, log2FC_A, log2FC_B, delta); Chem cell renderer for compound structures; top-impact ranking with configurable score rule; QC-status warning banner reading Story 2's tag.

**Should have (competitive differentiators):**
- Enrichment results carried with the published share (Scaffold doesn't do GO/KEGG; Datagrok already wires v1.2 enrichment with cross-DF protein highlight).
- "Request re-run with different parameters" reviewer-side button (mailto: shape).
- Pareto chart of which SPC metric flags most often.
- Enrichment-overlap panel in comparison view (terms in A only / B only / both).
- MCS / substructure highlight on the compound pair in comparison view (graceful degrade to plain SMILES).
- Comparison view itself publishable as a Story-1 artifact.

**Defer (v1.5+):**
- Reviewer comment threads; PDF report export; N-way (3+) compound comparison; CUSUM/EWMA charts; per-protein curated panel SPC; first-class target taxonomy; compound clustering / proteome-fingerprint SAR.

### Architecture Approach

Functional pipeline preserved: each step mutates a shared `DG.DataFrame` in-place and records its completion as a `proteomics.*` tag; downstream functions read those tags as preconditions. v1.4 adds two new top-level menu branches (`Campaign`, `Share`), one new entry under existing `Visualize` (`SPC Dashboard`), and **13 new `proteomics.*` tags** under three sub-prefixes (`published_*`, `campaign_*`, `spc_*`). See [ARCHITECTURE.md](./ARCHITECTURE.md).

**Major components:**
1. **`src/publishing/`** (NEW directory) — `publish-state.ts`, `trim-dataframe.ts`, `publish-project.ts`, `share-dialog.ts`. Sibling to `analysis/` because it reaches into project serialization + permissions UI.
2. **`src/analysis/spc.ts` + `src/viewers/spc-dashboard.ts`** (NEW files in existing dirs) — Shewhart I-MR math + Nelson rules engine + per-run QC extractor + tag helpers + viewer factory delegating to `DG.Viewer.lineChart` with formula lines for UCL/CL/LCL.
3. **`src/campaigns/`** (NEW directory) — `campaign-state.ts`, `campaign-types.ts`, `campaign-storage.ts`, `new-campaign-dialog.ts`, `impact-scoring.ts`, `compound-integration.ts`. Sibling to `analysis/` because it owns file I/O and multi-df coordination.
4. **`src/viewers/campaign-comparison.ts`** (NEW file) — `createSideBySideVolcano` reuses v1.3 `createVolcanoPlot` twice; `createCampaignDiffTable` builds union-joined diff DF; new `CampaignSelectionBus` abstraction with per-viewer subscription auto-eviction tied to viewer dispose (NOT a copy of v1.2's module-level `activeSubscriptions[]`).
5. **`src/panels/published-analysis-panel.ts`** (NEW file) — reviewer-side context panel registered via `@grok.decorators.panel` with semType filter on `Proteomics-ProteinId` and `isPublished(df)` first-line check.

**Convergence and conflict:**
- **All four converge** on four-phase Publishing → SPC → Campaign → Comparison ordering (15–18), zero new npm deps, tag-based campaign model, in-package SPC, and `grok.functions.call('Chem:...')` for Chem integration.
- **FEATURES vs ARCHITECTURE on the screening-run data model**: FEATURES wants a foundation pre-phase; ARCHITECTURE folds it into Phase 17. **Resolution: keep four phases; Phase 17 leads with Task 0 = data model before any viewer code.**
- **STACK vs ARCHITECTURE on SPC directory naming**: STACK proposes `src/spc/`; ARCHITECTURE keeps SPC math inside `src/analysis/spc.ts`. **Resolution: use ARCHITECTURE's shape** (consistent with v1.3 layout).

### Critical Pitfalls

Top five all in publishing or shared with it. See [PITFALLS.md](./PITFALLS.md).

1. **Pitfall 1 — Stale-snapshot leak (publish forgot to deep-clone the DataFrame).** Mitigation: `publishAnalysis(df, opts)` first line is `const frozen = df.clone(undefined, trimmedColumnNames);` — never publish the live `df`. Unit test mutates source after publish and asserts frozen is unchanged.
2. **Pitfall 3 — Project serialization drops `proteomics.*` tags / semTypes / df.name** (Phase 13 round-3 evidence at commit `e527d07ba1`). Mitigation: `assertPublishedShape(reopenedDf)` round-trip test in `src/tests/publish-roundtrip.ts`; belt-and-braces encode critical tag values in a single-row metadata column; mirror SEMTYPE constants in `detectors.js`.
3. **Pitfall 2 — Reviewer group accidentally gets Edit permission via Space-level inheritance.** Mitigation: publish ONLY the trimmed clone as a new entity into a target-named Space; post-grant verification via `grok.dapi.permissions.get(publishedProject)`; refuse to publish if Edit is found.
4. **Pitfall 13 — Subscription leak from copy-pasting v1.2 enrichment cross-link pattern.** Mitigation: build `CampaignSelectionBus` with per-viewer subscription auto-eviction tied to `viewer.onDispose`; unit test asserts subscription count returns to 0 after closing.
5. **Pitfall 9 — Mismatched protein populations across runs cause silent inner-join data loss.** Mitigation: full outer join on protein ID; "Quantified-In" column; "Run A: 4,800 / Run B: 4,200 / Both: 3,800" count panel.

Additional load-bearing pitfalls for requirements writer:
- **Pitfall 4** (versioning ambiguity): republish creates a NEW project (`publish_id`, `publish_version` tags + `superseded_by` pointer) — never an overwrite.
- **Pitfall 5** (per-protein SPC alarm flood): lock SPC metric set to RUN-LEVEL (one number per run × 4 metrics); default to Nelson rules 1+5 only.
- **Pitfall 7** (run-order ambiguity): SPC primary key is `(instrument_id, acquisition_datetime)` not import time.
- **Pitfall 11** (compound name drift): compound resolves to canonical identity at IMPORT time via explicit `compound_id` column or Chem canonicalization on SMILES.
- **Pitfall 14** (Cytokinetics demo audience contains biologists): jargon audit on every reviewer-side string; `TaskBarProgressIndicator` on every >500ms async path.

## Implications for Roadmap

Based on combined research, suggested phase structure continues v1.3 numbering:

### Phase 15: Read-Only Publishing Foundation
**Rationale:** Lowest implementation risk; highest demo value; establishes round-trip test infra that Phases 17–18 inherit. Must lead with the round-trip test as a gate — Pitfall 3 makes "tag survival across Project save" the load-bearing unknown.
**Delivers:** `src/publishing/{publish-state,trim-dataframe,publish-project,share-dialog}.ts`; `src/panels/published-analysis-panel.ts`; `Share` top-menu branch; 5 new tags (`proteomics.published`, `published_at`, `published_by`, `published_target`, `published_audit`); `src/tests/publish-roundtrip.ts`.
**Addresses:** Story 1 P1 — trimmed DataFrame, audit-context panel, frozen/no-edit, target-keyed filing, reviewer-group view-only ACL, versioned/dated artifact name, enrichment-carried-with-share.
**Avoids pitfalls:** 1, 2, 3, 4, 14.

### Phase 16: Sample-Level SPC Tracking
**Rationale:** Independent of Phases 15 and 17; can parallelize with Phase 15. Must land before Phase 17 because Phase 17 consumes `proteomics.spc_status`. SPC math is small (~390 LOC), pure TS, public-domain.
**Delivers:** `src/analysis/spc.ts`; `src/viewers/spc-dashboard.ts`; `Analyze → Compute SPC Status` + `Visualize → SPC Dashboard...` menu items; 3 new tags (`proteomics.spc_status`, `spc_metrics`, `spc_rules_tripped`); `src/tests/spc.ts`.
**Avoids pitfalls:** 5, 6, 7, 8.

### Phase 17: Campaign Data Model (with initial menu surfaces)
**Rationale:** Largest mental-model shift: file I/O under `System:AppData/Proteomics/campaigns/<id>/`, multi-df coordination, new menu branch, new SEMTYPE values, cross-package Chem coupling. Depends on Phase 16. **Task 0 = data model FIRST** before any viewer code.
**Delivers:** `src/campaigns/{campaign-state,campaign-types,campaign-storage,new-campaign-dialog,impact-scoring,compound-integration}.ts`; `Campaign` top-menu branch (`New Campaign...`, `Save Current Run to Campaign...`, `Open Campaign...`); 4 new tags (`proteomics.campaign`, `campaign_id`, `campaign_run_id`, `campaign_compound`); 2 new SEMTYPE values mirrored in `detectors.js`; AppData write/read; `src/tests/campaigns.ts`.
**Avoids pitfalls:** 11, 12, partial 9.

### Phase 18: Cross-Run Comparison Viewer
**Rationale:** Highest UX polish surface; gets spotlight in milestone demo. Needs ≥2 campaign runs from Phase 17; needs Phase 16's SPC status for QC warning banner.
**Delivers:** `src/viewers/campaign-comparison.ts`; `Campaign → Compare Runs...` menu item; `CampaignSelectionBus` abstraction; full-outer-join diff DF; Chem MCS substructure highlight via Pattern B with fallback; QC-status warning banner.
**Avoids pitfalls:** 9, 10, 13, 14.

### Phase Ordering Rationale

- **Publishing first (15)** because lowest implementation risk; isolates tag-stripping unknown to a small-scope phase; ships highest-demo-value feature for mid-milestone Cytokinetics checkpoint.
- **SPC second (16)** because math is pure-TS-deterministic; can parallelize with Phase 15; produces tag Phase 17 consumes.
- **Campaign data model third (17)** because it introduces the largest mental-model shift and depends on Phase 16's output. Task 0 = data model before viewers resolves FEATURES/ARCHITECTURE divergence.
- **Cross-run comparison last (18)** because it needs ≥2 campaign runs to be meaningful and carries highest UX polish stakes.
- **Cross-cutting reviewer-side UX hygiene** is a per-phase exit gate, not a separate phase.

### Research Flags

- **Phase 15 — Publishing**: REQUIRED research at plan time. The unknown is **which `proteomics.*` tags + `Proteomics-*` semTypes + `df.name` survive `DG.Project` save/reopen**. Plan-time research task: throw-away round-trip script publishes a fixture, reopens, enumerates surviving tags. Outcome decides whether belt-and-braces column-encoding mitigation is universal or selective.
- **Phase 17 — Campaign Data Model**: LIGHT research. AppData CSV pattern is well-precedented (HitTriage); the **compound column source convention** is a requirements-side decision (per-campaign SDF/CSV upload vs Cytokinetics-specific internal compound-DB connector vs free-text label).
- **Phase 16 — SPC**: standard patterns; research-phase skippable.
- **Phase 18 — Comparison**: standard patterns; plan-time decisions are scoping not research.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Every Datagrok API verified against repo HEAD; every Chem entrypoint verified; no-suitable-npm verdict verified; canonical publish/permissions/groups patterns verified against `packages/ApiSamples/scripts/dapi/`. |
| Features | MEDIUM-HIGH | Industry conventions verified (Spectronaut, Scaffold, MSstatsQC, SProCoP). Cytokinetics asks inferred from 2026-05-04 meeting notes via the three v1.4 todo files. Three open questions surfaced. |
| Architecture | HIGH | Existing v1.3 architecture established and well-documented. Reference packages inspected directly (HitTriage, Bio, Proteomics's Dendrogram pattern). 13-new-tag namespace verified non-colliding via grep. |
| Pitfalls | HIGH for codebase-grounded pitfalls; MEDIUM for Datagrok project-serialization surface; MEDIUM for SPC statistical traps. |

**Overall:** HIGH for four-phase ordering and stack/architecture/scope decisions. MEDIUM for two requirements-side decisions (target filing, compound source) and one technical unknown (tag survival).

### Gaps to Address

- **Tag survival across `DG.Project` save/reopen** (Phase 15 plan-time research).
- **Reviewer-group naming / target-keyed filing convention** (Phase 15 requirements decision).
- **Compound column source convention** (Phase 17 requirements decision).
- **Publish-block-on-QC-fail policy** (Phase 15 + Phase 16 cross-cut).
- **Optional R `qcc` SPC fallback path** (Phase 16, deferred to v1.5).

## Sources

See [STACK.md](./STACK.md), [FEATURES.md](./FEATURES.md), [ARCHITECTURE.md](./ARCHITECTURE.md), [PITFALLS.md](./PITFALLS.md) for full source citations.

---
*Research completed: 2026-06-06*
*Ready for roadmap: yes*
