# Milestone v1.4 Requirements: Cross-Team Review

**Milestone goal:** Enable proteomics experts to publish frozen analyses for biologist consumers, run longitudinal screening campaigns with weekly QC trust, and compare compound impact across runs — all within Datagrok.

**Customer driver:** Cytokinetics check-in defines the milestone close.

**Categories:**
- **PUB** — Read-only analysis publishing (Story 1)
- **SPC** — Sample-level statistical process control (Story 2)
- **CAMP** — Screening-campaign data model & cross-run comparison (Story 3)

---

## v1.4 Requirements

### PUB — Read-Only Analysis Publishing

**P1 (table stakes — required for credible demo)**

- [x] **PUB-01**: Proteomics expert can publish a frozen, read-only snapshot of a completed DE analysis as a Datagrok Project via `Proteomics → Share → Share Analysis for Review...`
- [x] **PUB-02**: Published snapshot contains only Protein ID, Gene, log2FC, p-value, adj.p-value, sig flag, and direction columns; raw intensities and peptide counts are excluded
- [x] **PUB-03**: Published snapshot is a deep clone of the source DataFrame — mutations to the source after publish do not leak to the snapshot (verified by `src/tests/publish-roundtrip.ts`)
- [x] **PUB-04**: Published snapshot is target-keyed via a free-text target identifier captured at publish time and persisted in `proteomics.published_target` tag
- [x] **PUB-05**: Reviewer group receives view-only ACL on the published Project; post-grant verification confirms no Edit permission was granted via Space inheritance
- [x] **PUB-06**: Volcano plot renders automatically on opening the published snapshot, with the stored FC and p-value thresholds drawn as formula lines
- [x] **PUB-07**: Reviewer-side context panel (`src/panels/published-analysis-panel.ts`) displays DE method, thresholds, group names, target, share date, and sharer's friendly name
- [x] **PUB-08**: Publishing dialog captures target, reviewer group, optional note, and shows a confirmation summary before creating the Project
- [x] **PUB-09**: Published Projects are named with a versioned/dated convention (e.g., `Proteomics-Review-<target>-v<N>-<date>`) that survives reopening
- [x] **PUB-10**: Re-publishing creates a new Project with a `proteomics.superseded_by` pointer on the previous version; never overwrites an existing published Project
- [x] **PUB-11**: Critical tags (`proteomics.published`, `published_at`, `published_by`, `published_target`, `published_audit`) survive `DG.Project` save/reopen via belt-and-braces encoding in a single-row metadata column

**P2 (differentiators)**

- [x] **PUB-12**: Published snapshot carries v1.2 enrichment results when present in the source DataFrame, with cross-DF protein-highlight selection wired
- [x] **PUB-13**: Reviewer-side "Request re-run with different parameters" button opens a mailto: addressed to the sharer with the published Project name in the subject

### SPC — Sample-Level Statistical Process Control

**P1 (table stakes)**

- [x] **SPC-01**: Per-run sample-level QC metrics (median log2 intensity, percent missing, control-replicate Pearson correlation, protein count above threshold) are computed when a run completes
- [x] **SPC-02**: SPC dashboard renders Shewhart I-chart and MR-chart per metric across a campaign's runs (via `DG.Viewer.lineChart` + formula lines for UCL/CL/LCL)
- [x] **SPC-03**: Nelson rules 1 (single point beyond 3σ) and 5 (2-of-3 beyond 2σ) are enabled by default; rules 2–4 and 6–8 are user-toggleable
- [x] **SPC-04**: Each run is tagged with `proteomics.spc_status` (`pass` / `flagged` / `out_of_control`) based on Nelson rule trips
- [x] **SPC-05**: SPC baseline is an explicit annotated subset of runs, locked at definition with iterative outlier removal applied — never rolling-recomputed
- [x] **SPC-06**: Run identity is `(instrument_id, acquisition_datetime)` captured at annotate-experiment time; not file-import time
- [x] **SPC-07**: Clicking a flagged point on the SPC chart drills into the run's source DataFrame and table view

**P2 (differentiators)**

- [x] **SPC-08**: SPC dashboard displays a Pareto chart showing which metric trips most often across the baseline window

### CAMP — Screening-Campaign Data Model & Cross-Run Comparison

**P1 (table stakes — data model)**

- [ ] **CAMP-01**: User can create a new screening campaign (target name, owner, optional notes) via `Proteomics → Campaign → New Campaign...`; persisted under `System:AppData/Proteomics/campaigns/<id>/`
- [ ] **CAMP-02**: User can save the currently-open analyzed run (DE complete) to an existing campaign with compound identity and run metadata via `Proteomics → Campaign → Save Current Run to Campaign...`
- [ ] **CAMP-03**: Each saved campaign run is tagged with `proteomics.campaign_id`, `campaign_run_id`, `campaign_compound`
- [ ] **CAMP-04**: Compound identity is canonicalized at save time — explicit `compound_id` from a column or Chem SMILES canonicalization (`Chem:convertMolNotation`); fuzzy-name matching is rejected
- [ ] **CAMP-05**: Campaign-index DataFrame opens in its own TableView via `Proteomics → Campaign → Open Campaign...`, listing each run with compound, acquisition date, sample count, protein count, SPC status, and top-hit summary
- [ ] **CAMP-06**: Compound column in the campaign-index uses the Chem cell renderer (`semType = 'Molecule'`, wired via `Chem:detectSmiles`)
- [ ] **CAMP-07**: Two new SEMTYPE values (`Proteomics-CompoundVid`, `Proteomics-SpcStatus`) are registered in `src/utils/proteomics-types.ts` and mirrored in `detectors.js`

**P1 (table stakes — comparison)**

- [ ] **CAMP-08**: User can select two runs from a campaign and open a side-by-side comparison view via `Proteomics → Campaign → Compare Runs...`
- [ ] **CAMP-09**: Comparison view renders two volcano plots side-by-side with shared X-axis (log2FC), independent Y-axis (p-value scale is sample-size-dependent), and shared FC/p threshold lines
- [ ] **CAMP-10**: Comparison view shows a diff table with full outer join on protein ID; columns include log2FC_A, log2FC_B, delta_log2FC, sig_in_A, sig_in_B, sig_in_both, and quantified_in (run list)
- [ ] **CAMP-11**: Comparison view shows a count panel: "Run A: N₁ proteins / Run B: N₂ / Both: N₃ / A only: N₄ / B only: N₅"
- [ ] **CAMP-12**: Selecting a protein in one volcano highlights the same protein in the other volcano and the diff table via a new `CampaignSelectionBus` with per-viewer subscription auto-eviction tied to viewer dispose
- [ ] **CAMP-13**: Comparison view shows a QC-status warning banner when either run's `proteomics.spc_status` is `flagged` or `out_of_control`
- [ ] **CAMP-14**: Comparison view shows a top-impact ranking with configurable scoring rule (significant-hit count / sum-magnitude-of-log2FC / top-N-strength)

**P2 (differentiators)**

- [ ] **CAMP-15**: Enrichment-overlap panel shows GO/KEGG/Reactome terms enriched in A only / B only / both (reuses v1.2 enrichment infrastructure)
- [ ] **CAMP-16**: Chem MCS substructure highlight on the compound pair in the comparison view header (via `DG.Func.find({package: 'Chem', name: 'mcsSearch'})` with arity guard); graceful degrade to plain SMILES if Chem MCS unavailable

---

## Future Requirements (deferred to v1.5+)

- Publish the cross-run comparison view itself as a read-only Story-1 artifact
- Reviewer-side comment threads on published snapshots (notifications, mentions — separate product surface)
- PDF report export from published snapshots (separate engineering project)
- N-way (3+) compound comparison view
- CUSUM / EWMA SPC charts (more sensitive than Shewhart, more tuning parameters)
- Per-protein curated panel SPC (needs FDR correction and multiplicity story)
- First-class target taxonomy (target entity, target search across analyses)
- Compound clustering / proteome-fingerprint SAR
- Optional R `qcc` SPC fallback path (`scripts/spc_qcc.R` mirroring `scripts/limma_de.R` cascade) — blocks on Datagrok R env availability
- Publish-block-on-QC-fail policy (gating publish dialog OK on `proteomics.spc_status`)

---

## Out of Scope

- **Cytokinetics-internal compound-DB connector** — per-customer customization; v1.4 ships per-campaign SDF/CSV upload + Chem canonicalization as the universal path
- **`@datagrok/chem` npm dependency** — cross-package `grok.functions.call('Chem:...')` chosen; avoids RDKit JS bundle bloat
- **Custom Postgres schema for campaigns** — tag-based AppData CSV chosen (HitTriage precedent); avoids coupling milestone ship date to server config
- **`@datagrok-libraries/cruddy` framework** — overkill for a campaign list with a custom comparison viewer; tag-based model preserves v1.3 architecture
- **R-first SPC** — TS-only client path (math is public-domain and deterministic); R `qcc` parked for v1.5
- **Email / Slack / webhook alerting** when a run flags — separate product surface, not v1.4-scoped
- **In-Datagrok comments or annotations** on published snapshots — Datagrok already supports rich comment surfaces at the platform layer; not v1.4-scoped
- **Fuzzy compound name matching** — silent data-quality risk; v1.4 requires canonical `compound_id` or canonical SMILES
- **Compound-specific assay metadata** beyond ID + SMILES (dose, concentration, vehicle batch) — v1.5+ scope

---

## Traceability

| REQ-ID | Phase | Plan | Status |
|--------|-------|------|--------|
| PUB-01 | Phase 15 | — | Pending |
| PUB-02 | Phase 15 | — | Pending |
| PUB-03 | Phase 15 | — | Pending |
| PUB-04 | Phase 15 | — | Pending |
| PUB-05 | Phase 15 | — | Pending |
| PUB-06 | Phase 15 | — | Pending |
| PUB-07 | Phase 15 | — | Pending |
| PUB-08 | Phase 15 | — | Pending |
| PUB-09 | Phase 15 | — | Pending |
| PUB-10 | Phase 15 | — | Pending |
| PUB-11 | Phase 15 | — | Pending |
| PUB-12 | Phase 15 | — | Pending |
| PUB-13 | Phase 15 | — | Pending |
| SPC-01 | Phase 16 | — | Pending |
| SPC-02 | Phase 16 | — | Pending |
| SPC-03 | Phase 16 | — | Pending |
| SPC-04 | Phase 16 | — | Pending |
| SPC-05 | Phase 16 | — | Pending |
| SPC-06 | Phase 16 | — | Pending |
| SPC-07 | Phase 16 | — | Pending |
| SPC-08 | Phase 16 | — | Pending |
| CAMP-01 | Phase 17 | — | Pending |
| CAMP-02 | Phase 17 | — | Pending |
| CAMP-03 | Phase 17 | — | Pending |
| CAMP-04 | Phase 17 | — | Pending |
| CAMP-05 | Phase 17 | — | Pending |
| CAMP-06 | Phase 17 | — | Pending |
| CAMP-07 | Phase 17 | — | Pending |
| CAMP-08 | Phase 18 | — | Pending |
| CAMP-09 | Phase 18 | — | Pending |
| CAMP-10 | Phase 18 | — | Pending |
| CAMP-11 | Phase 18 | — | Pending |
| CAMP-12 | Phase 18 | — | Pending |
| CAMP-13 | Phase 18 | — | Pending |
| CAMP-14 | Phase 18 | — | Pending |
| CAMP-15 | Phase 18 | — | Pending |
| CAMP-16 | Phase 18 | — | Pending |

**Coverage: 37/37 v1.4 requirements mapped (100%).** Phase mappings filled by roadmapper 2026-06-06; plan IDs filled by /gsd:plan-phase as plans are drafted.

---
*Created: 2026-06-06 for milestone v1.4 Cross-Team Review*
*Phase mappings added: 2026-06-06*
