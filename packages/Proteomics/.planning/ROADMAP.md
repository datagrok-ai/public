# Roadmap: Proteomics Package

## Milestones

- ✅ **v1.0 MVP** -- Phases 1-5 (shipped 2026-03-05)
- ✅ **v1.1 Parser & QC** -- Phases 6-7 (shipped 2026-03-07)
- ✅ **v1.2 Biological Interpretation** -- Phases 8-9 (shipped 2026-03-07)
- ✅ **v1.3 Richer Dialogs & UX** -- Phases 10-14 + 999.1 + 999.4 (shipped 2026-06-05; initial cut 2026-03-10 with Phases 10-11, extended through Phases 12/13/14 + two hotfixes May–June 2026)
- 🚧 **v1.4 Cross-Team Review** -- Phases 15-18 (started 2026-06-06)

## Phases

<details>
<summary>✅ v1.0 MVP (Phases 1-5) -- SHIPPED 2026-03-05</summary>

- [x] Phase 1: Data Import and Foundation (2/2 plans)
- [x] Phase 2: Analysis Pipeline (3/3 plans)
- [x] Phase 3: Visualization (2/2 plans)
- [x] Phase 4: Annotation and Alternatives (2/2 plans)
- [x] Phase 5: Gap Closure & Hardening (3/3 plans)

</details>

<details>
<summary>✅ v1.1 Parser & QC (Phases 6-7) -- SHIPPED 2026-03-07</summary>

- [x] Phase 6: Generic Matrix Parser (2/2 plans) -- completed 2026-03-06
- [x] Phase 7: QC Dashboard (3/3 plans) -- completed 2026-03-07

</details>

<details>
<summary>✅ v1.2 Biological Interpretation (Phases 8-9) -- SHIPPED 2026-03-07</summary>

- [x] Phase 8: Gene ID Mapping & Enrichment Analysis (2/2 plans) -- completed 2026-03-07
- [x] Phase 9: Enrichment Visualization & Volcano Integration (2/2 plans) -- completed 2026-03-07

</details>

<details>
<summary>✅ v1.3 Richer Dialogs & UX (Phases 10-14 + 999.1 + 999.4) -- SHIPPED 2026-06-05</summary>

- [x] Phase 10: Spectronaut Parser and Core Algorithms (2/2 plans) -- completed 2026-03-07
- [x] Phase 11: Dialog Expansion and UX Polish (4/4 plans) -- completed 2026-03-08
- [x] Phase 12: Spectronaut Input Coverage (4/4 plans) -- completed 2026-05-15
- [x] Phase 13: CK-omics Volcano and Enrichment Parity (6/6 + 4 gap plans) -- completed 2026-06-04
- [x] Phase 14: CK-omics Analyst-Experience Enhancements (5/5 plans) -- completed 2026-06-02
- [x] Phase 999.1 hotfix: Annotate Experiment dialog prefill (1/1) -- completed 2026-06-04
- [x] Phase 999.4 hotfix: Spectronaut malformed-counter parity (1/1) -- completed 2026-06-05

Full milestone details: [.planning/milestones/v1.3-ROADMAP.md](milestones/v1.3-ROADMAP.md)

</details>

### 🚧 v1.4 Cross-Team Review (Phases 15-18) -- ACTIVE

- [x] **Phase 15: Read-Only Publishing Foundation** -- Trimmed `DG.Project` publish flow, target-keyed filing, reviewer-group view-only ACL, round-trip survival of `proteomics.*` tags (completed 2026-06-08)
- [x] **Phase 16: Sample-Level SPC Tracking** -- Per-run Shewhart I-MR + Nelson rules engine, run identity = `(instrument_id, acquisition_datetime)`, locked baseline, `proteomics.spc_status` tag (completed 2026-06-08)
- [ ] **Phase 17: Campaign Data Model** -- AppData campaign storage, compound canonicalization, two new SEMTYPE values, `Campaign` top-menu branch (Task 0 = data model before viewers)
- [ ] **Phase 18: Cross-Run Comparison Viewer** -- Side-by-side volcano, full-outer-join diff table, `CampaignSelectionBus` with per-viewer eviction, QC-status warning, optional Chem MCS

## Phase Details

### Phase 15: Read-Only Publishing Foundation
**Goal**: A proteomics expert can publish a frozen, trimmed snapshot of a completed DE analysis as a Datagrok Project that a reviewer group can open read-only and see exactly the columns and thresholds the expert intended.
**Depends on**: v1.3 (DE pipeline + volcano + enrichment shipped)
**Requirements**: PUB-01, PUB-02, PUB-03, PUB-04, PUB-05, PUB-06, PUB-07, PUB-08, PUB-09, PUB-10, PUB-11, PUB-12, PUB-13
**Success Criteria** (what must be TRUE):
  1. Expert can run `Proteomics → Share → Share Analysis for Review...` on a DE-complete table, fill in target + reviewer group + optional note, and produce a versioned/dated `DG.Project` containing ONLY Protein ID, Gene, log2FC, p-value, adj.p-value, sig, and direction columns -- raw intensities and peptide counts are absent.
  2. A reviewer in the named group can open the published Project and see a volcano plot rendered on first paint with the stored FC and p-value thresholds drawn as formula lines, plus a context panel showing DE method, thresholds, group names, target, share date, and sharer's friendly name.
  3. The `assertPublishedShape` round-trip test passes: every required `proteomics.published*` tag, the trimmed allowlist of columns, `df.name`, and the volcano viewer's dock position survive `DG.Project.save` followed by re-open in a fresh session (Pitfall 3 mitigation -- belt-and-braces metadata column carries the critical tag values too).
  4. Mutating the live source DataFrame after publish (rerun DE, drop a column, change a tag) leaves the published clone unchanged when re-opened -- verified by `src/tests/publish-roundtrip.ts`.
  5. Re-publishing the same analysis creates a NEW project with `proteomics.superseded_by` on the previous version; post-grant `grok.dapi.permissions.get` shows the reviewer group with View access only, never Edit, on any published Project.
**Plans**: 9 plans
  - [x] 15-00-spike-publish-roundtrip-enumeration-PLAN.md — Wave 0 spike enumerating which tags / semTypes / df.name / Project.options / viewer dock positions survive DG.Project save+reopen
  - [x] 15-01-publish-state-PLAN.md — Tag helpers + PUBLISHED_TAGS / META_COLUMNS / PublishedMetadata / slugifyTarget / findPriorShare
  - [x] 15-02-trim-dataframe-PLAN.md — trimForPublish + trimEnrichmentForPublish (deep-clone + allowlist + belt-and-braces metadata columns)
  - [x] 15-03-assert-published-shape-PLAN.md — Round-trip contract assertion (consumed by Plan 04 orchestrator AND Plan 08 test)
  - [x] 15-04-publish-project-PLAN.md — publishAnalysis orchestrator (9 steps incl. 2 non-negotiable gates: verify-and-rollback + assertPublishedShape)
  - [x] 15-05-share-dialog-PLAN.md — showShareForReviewDialog (target + group ChoiceInput + note + republish banner) + buildMailtoUrl (P2)
  - [x] 15-06-published-analysis-panel-PLAN.md — Reviewer-side audit panel with isPublished guard + column-first read + supersede link + mailto (P2)
  - [x] 15-07-package-wireup-PLAN.md — Register Share menu + panel decorator + test import line in src/package.ts and src/package-test.ts
  - [x] 15-08-publish-roundtrip-test-PLAN.md — 5 regression tests covering load-bearing gates per ROADMAP success criterion 3
**UI hint**: yes

### Phase 16: Sample-Level SPC Tracking
**Goal**: A proteomics expert running weekly assays can compute and inspect per-run statistical-process-control metrics across a campaign so a downstream cross-run comparison knows which runs are in control before combining them.
**Depends on**: v1.3 (intensity columns + `proteomics.groups`); independent of Phase 15
**Requirements**: SPC-01, SPC-02, SPC-03, SPC-04, SPC-05, SPC-06, SPC-07, SPC-08
**Success Criteria** (what must be TRUE):
  1. Running `Analyze → Compute SPC Status` on any DE-complete or annotated run computes the four standard metrics (median log2 intensity, percent missing, control-replicate Pearson correlation, protein count above threshold) and writes them onto the DataFrame as `proteomics.spc_metrics` + a `proteomics.spc_status` of `pass` / `flagged` / `out_of_control`.
  2. `Visualize → SPC Dashboard...` renders a Shewhart I-chart + MR-chart per metric across a campaign's runs with UCL/CL/LCL formula lines, ordered by `(instrument_id, acquisition_datetime)` -- not file-import time, not row-insertion order.
  3. Nelson rules 1 (single point beyond 3σ) and 5 (2-of-3 beyond 2σ) trip a run by default; rules 2-4 and 6-8 are user-toggleable per metric, with `proteomics.spc_rules_tripped` JSON recording which rules fired on flagged runs.
  4. The SPC baseline is an explicit operator-annotated subset of runs locked at definition time with iterative outlier removal applied -- never silently rolling-recomputed when new runs land; rebuilding requires an explicit user action.
  5. Clicking a flagged point on the SPC chart drills into the originating run's source DataFrame in its own TableView, and a Pareto chart shows which metric trips most often across the baseline window.
**Plans**: 7 plans
  - [x] 16-01-PLAN.md — Wave 0 spike: formula-line dashed-style on datetime X + RED test scaffold (28 SPC tests under category SPC)
  - [x] 16-02-PLAN.md — `src/analysis/spc.ts` math + tag helpers (4 metrics, Nelson 1-8 rules engine, status classification, RunMeta helpers)
  - [x] 16-03-PLAN.md — `src/analysis/spc-storage.ts` AppData I/O (runs.csv upsert, baseline JSON, slugify, iterative outlier removal, schema version)
  - [x] 16-04-PLAN.md — `src/parsers/spectronaut-parser.ts` augmentation: capture R.RunDate + R.InstrumentMethod → `proteomics.spc_run_meta_seed`
  - [x] 16-05-PLAN.md — `src/analysis/experiment-setup.ts` dialog extension + `src/package.ts` Compute SPC Status orchestrator + SPC Dashboard menu stub
  - [x] 16-06-PLAN.md — `src/viewers/spc-dashboard.ts` Dashboard + Define/Rebuild baseline modal + sidebar + drill-down (P1 surface complete)
  - [x] 16-07-PLAN.md — `src/viewers/spc-dashboard.ts` Pareto chart extension (SPC-08, P2, deferrable to v1.4.1 hotfix)
**UI hint**: yes

### Phase 17: Campaign Data Model
**Goal**: A proteomics expert can persist analyzed runs into a named screening campaign with canonical compound identity so the cross-run comparison viewer (Phase 18) has a stable, indexed corpus to read from. Task 0 of this phase is the data model -- tags + SEMTYPE + AppData schema -- before any new viewer code lands.
**Depends on**: Phase 16 (campaign index reads `proteomics.spc_status` per run)
**Requirements**: CAMP-01, CAMP-02, CAMP-03, CAMP-04, CAMP-05, CAMP-06, CAMP-07
**Success Criteria** (what must be TRUE):
  1. `Proteomics → Campaign → New Campaign...` creates a campaign under `System:AppData/Proteomics/campaigns/<id>/` with target name, owner, and notes; `Save Current Run to Campaign...` saves a DE-complete run with compound identity + run metadata; both stamp `proteomics.campaign_id`, `campaign_run_id`, `campaign_compound` onto the saved DataFrame.
  2. Compound identity is canonicalized at save time -- either from an explicit `compound_id` column or via Chem SMILES canonicalization (`Chem:convertMolNotation`); fuzzy name matching is rejected so the same compound across weeks always resolves to the same identity.
  3. `Open Campaign...` opens a campaign-index DataFrame in its own TableView with one row per saved run showing compound (rendered via Chem `Molecule` cell renderer), acquisition date, sample count, protein count, SPC status, and a top-hit summary; the compound column auto-wires through `Chem:detectSmiles`.
  4. Two new SEMTYPE values (`Proteomics-CompoundVid`, `Proteomics-SpcStatus`) are registered in `src/utils/proteomics-types.ts` AND mirrored in `detectors.js`; opening a campaign-index from a clean session re-detects both column types without manual tagging.
  5. The AppData schema (campaign JSON + per-run CSV layout under `campaigns/<id>/`) is locked at Phase 17 entry so Phase 18 is purely additive against it -- no schema migration mid-milestone.
**Plans**: TBD
**UI hint**: yes

### Phase 18: Cross-Run Comparison Viewer
**Goal**: A chemist or proteomics expert can pick two saved runs from a campaign and see, side by side, how each compound moved the proteome, including which proteins shifted in only one run vs both, with a QC warning if either run is flagged.
**Depends on**: Phase 17 (≥1 campaign with ≥2 saved runs); reads Phase 16's `proteomics.spc_status` for the QC banner
**Requirements**: CAMP-08, CAMP-09, CAMP-10, CAMP-11, CAMP-12, CAMP-13, CAMP-14, CAMP-15, CAMP-16
**Success Criteria** (what must be TRUE):
  1. `Proteomics → Campaign → Compare Runs...` lets the user pick two runs from a campaign and opens a comparison view with two volcano plots side by side -- shared X-axis (log2FC), independent Y-axis (p-value scale is sample-size-dependent), shared FC/p threshold lines, distinct per-run color/title.
  2. The comparison view shows a full-outer-join diff DataFrame on protein ID (columns: log2FC_A, log2FC_B, delta_log2FC, sig_in_A, sig_in_B, sig_in_both, quantified_in) plus a visible count panel: "Run A: N₁ / Run B: N₂ / Both: N₃ / A only: N₄ / B only: N₅" -- so reviewers see the population asymmetry instead of an inner-join-style silent loss.
  3. Selecting a protein in one volcano highlights the same protein in the other volcano and the diff table via a new `CampaignSelectionBus` whose per-viewer subscriptions auto-evict when a viewer is disposed -- closing the comparison view returns the bus subscription count to 0 (not a copy of v1.2's `activeSubscriptions[]` global).
  4. A configurable top-impact ranking surfaces (significant-hit count / sum-magnitude-of-log2FC / top-N-strength rules); when either run's `proteomics.spc_status` is `flagged` or `out_of_control`, a QC-status warning banner renders above the side-by-side panel so reviewers don't compare drifting runs without knowing.
  5. An enrichment-overlap panel reuses v1.2 enrichment infrastructure to show GO/KEGG/Reactome terms enriched in A only / B only / both; the comparison header optionally highlights the MCS substructure on the compound pair via `DG.Func.find({package: 'Chem', name: 'mcsSearch'})` with arity guard and graceful degrade to plain SMILES when Chem MCS is unavailable.
**Plans**: TBD
**UI hint**: yes

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 15. Read-Only Publishing Foundation | 9/9 | Complete    | 2026-06-08 |
| 16. Sample-Level SPC Tracking | 7/7 | Complete    | 2026-06-08 |
| 17. Campaign Data Model | 0/0 | Not started | - |
| 18. Cross-Run Comparison Viewer | 0/0 | Not started | - |

All v1.0–v1.3 phases shipped. v1.4 phases 15-18 are the active milestone; plan-phase will fill the per-phase plan counts.

## Backlog

### Phase 999.2: Pre-normalized banner wording — "is" not "may be" (BACKLOG)

**Goal:** [Captured for future planning] The pre-normalized warning banner (src/analysis/normalization.ts:192) reads "This data may be pre-normalized (Spectronaut). Additional normalization may distort results." but it only renders when df.getTag('proteomics.preNormalized') === 'true' (normalization.ts:189) — a definite assertion. Change "This data may be pre-normalized (Spectronaut)." to "This data is pre-normalized (Spectronaut)." Keep the trailing "Additional normalization may distort results." Severity: COSMETIC. Surfaced by Phase 12 UAT (12-UAT.md gap obs-B); pre-existing, out of Phase-12 scope (normalization.ts untouched by Phase 12). Explicitly NOT pulled into v1.4 -- kept tight per user direction.
**Requirements:** TBD
**Plans:** 0 plans

Plans:
- [ ] TBD (promote with /gsd:review-backlog when ready)

### Phase 999.3: Publish read-only proteomics analysis for biologist review (SUPERSEDED BY v1.4)

**Status:** Promoted into v1.4 Phase 15 as `PUB-01..13`. The strategic spine described here is now the active milestone goal; the user-stories and the trim-contract design are embodied in v1.4 REQUIREMENTS.md.
**Requirements:** PUB-01..13 (Phase 15)
**Plans:** Tracked under Phase 15 once /gsd:plan-phase runs.
