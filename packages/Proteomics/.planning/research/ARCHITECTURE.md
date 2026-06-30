# Architecture Research — v1.4 Cross-Team Review

**Domain:** Subsequent-milestone integration into the v1.3 Proteomics package (mass-spec analysis pipeline)
**Researched:** 2026-06-06
**Confidence:** HIGH (existing package architecture is established; reference packages HitTriage / Bio / Chem inspected directly)

This document answers: *how do read-only publishing, longitudinal SPC tracking, and cross-run comparison integrate with the existing v1.3 Proteomics architecture?* It is **not** a greenfield architecture — it is the contract for the deltas v1.4 adds on top of v1.3.

---

## 1. Architecture Snapshot Going Into v1.4

```
┌──────────────────────────────────────────────────────────────────────┐
│  Datagrok top-menu: Proteomics | Import | Annotate | Analyze |       │
│                     Visualize | Enrichment                           │
│  Entry: src/package.ts  (PackageFunctions class)                     │
└──────────────────┬──────────────────────┬────────────────────────────┘
                   │                      │
                   ▼                      ▼
        ┌──────────────────┐    ┌────────────────────┐
        │   src/parsers/   │    │   src/analysis/    │  mutate df in place
        │   maxquant       │    │  experiment-setup  │  set proteomics.*
        │   spectronaut    │    │  normalization     │  tag on completion
        │   spectronaut-   │    │  imputation        │
        │     candidates   │    │  differential-expr │
        │   fragpipe       │    │  pca               │
        │   generic        │    │  enrichment        │
        └──────┬───────────┘    └────────┬───────────┘
               │                         │
               ▼                         ▼
        ┌──────────────────────────────────────┐
        │       DG.DataFrame (one)             │
        │  rows = proteins                     │
        │  tags = proteomics.source, .groups,  │
        │    .normalized, .imputed,            │
        │    .de_complete, .de_method,         │
        │    .enrichment, .preNormalized       │
        └──────────────┬───────────────────────┘
                       │
                       ▼
        ┌──────────────────┐    ┌────────────────────┐
        │   src/viewers/   │    │   src/panels/      │
        │  volcano         │    │  uniprot-panel     │  semType-triggered
        │  heatmap         │    │                    │  context widget
        │  pca-plot        │    └────────────────────┘
        │  qc-dashboard    │
        │  enrichment-vw   │
        └──────────────────┘
```

**Load-bearing invariants v1.4 must respect:**
1. Pipeline state is **only** on the DataFrame, as `proteomics.*` string tags. No external store.
2. Functions follow `showX` / `openX` / `createX` / `parseX` prefix contracts (see `CLAUDE.md`).
3. New `Proteomics-*` semantic types live in **`src/utils/proteomics-types.ts` AND `detectors.js`** (mirrored — both are sources of truth for different platform load-stages).
4. R-script calls always have a client-side fallback.
5. Cross-package dependencies use `DG.Func.find({package, name})[0]` with arity/existence guards (the Dendrogram pattern in `viewers/heatmap.ts:138-145`).
6. PCA is the only place today that produces a sibling DataFrame opened in its own view. v1.4 will multiply this — sibling DataFrames become a first-class pattern.

---

## 2. New Capability → Existing Architecture Mapping

| New v1.4 capability | Best-fit existing layer | New layer needed? |
|---|---|---|
| Read-only Project publish (target-keyed) | **`src/analysis/`** for trim (mutates df: clones + drops intensity cols) + new files for Project serialize + permissions UI | **`src/publishing/`** — sibling to `analysis/`. Owns trim + Project save + permissions + reviewer-side affordances. Uses the v1.3 `showX/openX/createX` naming. |
| Sample-level SPC tracking | **`src/analysis/`** for math (Shewhart/Nelson rules, intensity median/missingness/control-correlation) + **`src/viewers/`** for the dashboard | **No new layer.** Math goes in `src/analysis/spc.ts`; viewer goes in `src/viewers/spc-dashboard.ts`. |
| Screening-campaign data model + cross-run comparison | Cross-cuts multiple v1.3 layers and owns AppData file I/O — breaks the "one df in memory" assumption | **`src/campaigns/`** — sibling to `analysis/`. Modeled after HitTriage's `src/app/` campaign code. |

**Rationale for the asymmetry:** publishing and SPC operate on a single DataFrame and fit the existing "mutate + tag + viewer" pattern. Campaigns are a step-change — they introduce `campaign.json` + multiple data CSVs in AppData and a new "open a saved campaign" navigation surface. They deserve their own directory. Publishing is in between: small enough to fit in `analysis/`, but reaches into project serialization + permissions UI that don't belong there, so a dedicated directory is cleaner.

---

## 3. Tag Contract Extensions

All eight existing tags keep their current semantics. v1.4 adds these. Tag values follow the existing v1.3 convention (`'true'` for boolean markers; JSON for structured data).

| New tag | Value | When set | Why |
|---|---|---|---|
| `proteomics.published` | `'true'` | At the moment the trimmed clone is built by `publishAnalysisForReview()`. Set on the **trimmed clone**, not the source df. | Lets downstream menu handlers + viewers detect they're inside a published Project and refuse to mutate it. |
| `proteomics.published_at` | ISO-8601 timestamp string | Same time as `proteomics.published` | Audit trail; surfaces in the reviewer-side "Analysis Context" panel. |
| `proteomics.published_by` | `grok.shell.user.id` (UUID string) | Same time as `proteomics.published` | Audit trail; reviewer-side panel resolves to friendly name via `grok.dapi.users.find(id)`. |
| `proteomics.published_target` | Target identifier string (free-text label, gene symbol, or UniProt accession) | Same time as `proteomics.published` | Target-keyed filing — the value the share dialog captures and the project list filters by. |
| `proteomics.published_audit` | JSON `{de_method, fc_threshold, p_threshold, group1_name, group2_name, n_proteins, n_significant, normalize_method, impute_method}` | Same time as `proteomics.published` | Captures the analysis context that produced the trimmed numbers; surfaced in reviewer's "Audit Context" panel. Encoded once at publish time because the source df may be re-analyzed afterward. |
| `proteomics.campaign` | `'true'` | Set on a DataFrame the moment it's claimed by a campaign run (either created fresh by `showNewCampaignRunDialog` or attached via `attachCampaignRun`). | Lets DE / Volcano / Heatmap handlers know this df is **inside a campaign** and offer "Save Run to Campaign" affordances. |
| `proteomics.campaign_id` | Campaign UUID (matches the `campaignId` directory under AppData) | Set together with `proteomics.campaign` | Joins the run df back to its campaign-index df + sibling runs. |
| `proteomics.campaign_run_id` | Run UUID inside the campaign | Set together with `proteomics.campaign` | Used as the AppData filename (`<campaign_id>/runs/<run_id>.csv`) and the row key in the campaign-index df. |
| `proteomics.campaign_compound` | JSON `{vid: string, smiles: string\|null, name: string, week?: string}` | Set when the user picks a compound during `showNewCampaignRunDialog`. Always set when `proteomics.campaign` is set. | Per-run compound identity. SMILES enables Chem-package rendering in the cross-run comparison view. `vid` is the human-readable name. |
| `proteomics.spc_status` | `'pass'` \| `'flag'` \| `'fail'` | Set by `computeRunSpcStatus()` at the end of each campaign run (or on-demand for any df via QC menu). | Per-run pass/flag/fail traffic light. `'flag'` = Nelson rules tripped but not catastrophic; `'fail'` = hard limits exceeded. |
| `proteomics.spc_metrics` | JSON `{median_intensity, missingness_pct, control_correlation, protein_count, sample_count, computed_at}` | Same time as `proteomics.spc_status`. | The four standard SPC metrics surfaced in the SPC dashboard. Persisted so the campaign-index df can read them without re-opening each run. |
| `proteomics.spc_rules_tripped` | JSON array of rule names: `['nelson_1', 'shewhart_3sigma', ...]` | Same time as `proteomics.spc_status` (only when status ≠ `'pass'`). | Auditable trail for why a run was flagged. Empty/missing on pass. |
| `proteomics.campaign_index` | `'true'` | Set on the campaign-index DataFrame produced by `loadCampaign()` (rows = runs, columns = run_id, week, compound_vid, compound_smiles, spc_status, n_significant, top_protein, ...). | Marker that this is the index df, not a run df. Volcano/Heatmap viewers refuse to render against an index df. |

**Naming convention:** `proteomics.<feature>[_<sub-field>]` with snake_case sub-fields. Matches existing `proteomics.de_method`, `proteomics.de_complete`, `proteomics.preNormalized` (the lone camelCase outlier — keep as-is for compatibility).

**No collisions verified:** grepped all `setTag/getTag` calls in `src/`. Existing tags are `de_complete`, `de_method`, `enrichment`, `enrichment_smart_filtered{,_cap,_dropped_parents,_kept,_total}`, `groups`, `imputed`, `location_acc_hash`, `normalized`, `preNormalized`, `source`. All 13 new tags are unique.

**Setter discipline:** v1.4 adds these as helpers, not raw `setTag()` calls scattered through code. Mirror the `getGroups(df)` / `setGroups(df, g)` pattern already in `analysis/experiment-setup.ts`:

```ts
// src/publishing/publish-state.ts
export function markPublished(df: DG.DataFrame, target: string, audit: PublishAudit): void
export function getPublishedAudit(df: DG.DataFrame): PublishAudit | null
export function isPublished(df: DG.DataFrame): boolean

// src/campaigns/campaign-state.ts
export function markCampaignRun(df: DG.DataFrame, campaignId: string, runId: string, compound: Compound): void
export function getCampaignRun(df: DG.DataFrame): { campaignId: string; runId: string; compound: Compound } | null
export function isCampaignRun(df: DG.DataFrame): boolean
export function isCampaignIndex(df: DG.DataFrame): boolean

// src/analysis/spc.ts
export function setSpcStatus(df: DG.DataFrame, status: SpcStatus, metrics: SpcMetrics, rules: string[]): void
export function getSpcStatus(df: DG.DataFrame): { status: SpcStatus; metrics: SpcMetrics; rules: string[] } | null
```

---

## 4. Directory Structure (Delta from v1.3)

```
packages/Proteomics/src/
├── package.ts                  ← extended: ~12 new @grok.decorators.func entries
├── package-test.ts             ← imports 3 new test files
├── analysis/
│   ├── (v1.3 files unchanged)
│   └── spc.ts                  ← NEW: SPC math + setSpcStatus/getSpcStatus
├── parsers/                    ← unchanged
├── viewers/
│   ├── (v1.3 files unchanged)
│   ├── spc-dashboard.ts        ← NEW: openSpcDashboard, createSpcControlChart
│   └── campaign-comparison.ts  ← NEW: createSideBySideVolcano, createCampaignDiffTable
├── panels/
│   ├── uniprot-panel.ts        (v1.3 unchanged)
│   └── published-analysis-panel.ts  ← NEW: reviewer-side audit-context panel
├── publishing/                 ← NEW DIRECTORY: read-only-publish responsibility
│   ├── publish-state.ts        ← markPublished/getPublishedAudit/isPublished
│   ├── trim-dataframe.ts       ← buildTrimmedReviewDf (clone + drop intensity cols)
│   ├── publish-project.ts      ← publishAnalysisForReview (Project.create + permissions)
│   └── share-dialog.ts         ← showShareForReviewDialog (target + reviewer group + note)
├── campaigns/                  ← NEW DIRECTORY: screening-campaign responsibility
│   ├── campaign-state.ts       ← markCampaignRun/getCampaignRun/isCampaignRun
│   ├── campaign-types.ts       ← Campaign, CampaignRun, Compound, ImpactScore interfaces
│   ├── campaign-storage.ts     ← loadCampaign, saveCampaign, listCampaigns (AppData I/O)
│   ├── new-campaign-dialog.ts  ← showNewCampaignDialog (create) + showNewCampaignRunDialog
│   ├── impact-scoring.ts       ← computeImpactScore (count/magnitude/top-N rules)
│   └── compound-integration.ts ← Chem-package bridge (SMILES rendering, MCS via DG.Func.find)
├── utils/
│   ├── proteomics-types.ts     ← extended with new SEMTYPE values (see §5)
│   └── column-detection.ts     ← unchanged
└── tests/
    ├── (v1.3 files unchanged)
    ├── publishing.ts           ← NEW: trim-roundtrip, permissions, audit-tag tests
    ├── spc.ts                  ← NEW: Shewhart/Nelson rule tests with synthetic series
    └── campaigns.ts            ← NEW: campaign save/load, impact scoring, comparison view

packages/Proteomics/detectors.js  ← extended: 2 new detector functions (see §5)
packages/Proteomics/files/        ← unchanged
```

**Why `src/publishing/` and `src/campaigns/` are siblings to `analysis/`, not subfolders:**
- They cross-cut multiple v1.3 concerns. Publishing reads results from `analysis/` + builds viewers from `viewers/` + creates a new Project entity (a platform concern, not an analysis concern).
- Campaigns own *file I/O* (AppData CSV/JSON read/write) which `analysis/` doesn't.
- Keeping them as siblings preserves the v1.3 mental model: `analysis/` mutates the active df; `publishing/` and `campaigns/` operate on a different abstraction layer (a Project / a Campaign).

**Why NOT a new directory for SPC:**
- SPC math operates on the same df via the same "mutate + tag" pattern as normalization. It belongs in `analysis/`.
- The SPC dashboard reads tags and produces viewers — that's the `viewers/` pattern.
- Adding `src/spc/` would split the SPC story across two directories for no architectural gain.

---

## 5. Semantic Type Additions

Two new semantic types are needed; both follow the existing `Proteomics-*` convention. Both must be added in **two places**: `src/utils/proteomics-types.ts` (TypeScript constant) and `detectors.js` (runtime auto-detection). This is the rule documented in `CLAUDE.md`.

| New SEMTYPE key | String value | Used on | Detector heuristic |
|---|---|---|---|
| `SEMTYPE.COMPOUND_VID` | `'Proteomics-CompoundVid'` | Per-run identifier column in campaign-index df, and in compound annotation column of run dfs | Column name contains `vid` or `compound_id` and is unique-per-row; skip detector when ambiguous — let the user tag manually |
| `SEMTYPE.SPC_STATUS` | `'Proteomics-SpcStatus'` | The status column in campaign-index df (one row per run) | Column name matches `spc_status` / `qc_status` AND all values are in `{pass, flag, fail}` |

**SMILES is intentionally NOT given a Proteomics-specific SEMTYPE.** The platform / Chem package already owns `DG.SEMTYPE.MOLECULE` with detectors. v1.4 reads `proteomics.campaign_compound.smiles` from the JSON tag and sets `col.semType = DG.SEMTYPE.MOLECULE` on the SMILES column in the comparison-view df, so the Chem cell renderer takes over automatically.

---

## 6. `package.ts` Menu Structure

The v1.3 menu has five top-level branches: `Import`, `Annotate Experiment`, `Analyze`, `Visualize`, `Enrichment Analysis`. v1.4 adds **two** new top-level branches (`Campaign`, `Share`) and **one** new entry under the existing `Visualize` branch (`SPC Dashboard`). Rationale: the new branches each open a *different kind of artifact* (a campaign index, a Project) that doesn't fit under the existing pipeline-oriented branches.

```
Proteomics
├── Import
│   ├── MaxQuant...                                      (v1.3)
│   ├── Spectronaut Report...                            (v1.3)
│   ├── Spectronaut Candidates...                        (v1.3)
│   ├── FragPipe...                                      (v1.3)
│   └── Generic Matrix...                                (v1.3)
├── Annotate Experiment...                               (v1.3)
├── Analyze
│   ├── Normalize...                                     (v1.3)
│   ├── Impute Missing Values...                         (v1.3)
│   ├── Differential Expression...                       (v1.3)
│   └── Compute SPC Status                               ← NEW v1.4 (immediate; no "..." suffix)
├── Visualize
│   ├── Volcano Plot...                                  (v1.3)
│   ├── Volcano Options...                               (v1.3)
│   ├── Heatmap...                                       (v1.3)
│   ├── PCA...                                           (v1.3)
│   ├── Group-Mean Correlation...                        (v1.3)
│   ├── QC Dashboard...                                  (v1.3)
│   ├── SPC Dashboard...                                 ← NEW v1.4
│   ├── Show All Visualizations...                       (v1.3)
│   └── Enrichment Charts...                             (v1.3)
├── Enrichment Analysis...                               (v1.3)
├── Campaign                                             ← NEW v1.4 BRANCH
│   ├── New Campaign...                                  ← creates campaign + first run
│   ├── Save Current Run to Campaign...                  ← attaches current df to existing campaign
│   ├── Open Campaign...                                 ← navigates to campaign-index df
│   └── Compare Runs...                                  ← side-by-side volcano + diff table
└── Share                                                ← NEW v1.4 BRANCH
    ├── Share Analysis for Review...                     ← target + reviewer-group dialog
    └── Open Published Analysis...                       ← project list filtered to proteomics shares
```

**Naming follows the `...` convention** documented in `CLAUDE.md`: menu items that open a dialog end with `...`; immediate actions do not. `Compute SPC Status` is immediate (uses defaults — current run, all metrics) because it's a single-button QC check; the parameterized variant lives inside `SPC Dashboard...`.

**Precondition tag gating** (in the menu handler — same pattern as v1.3's `proteomics.de_complete` check before Volcano):

| Menu item | Required tag |
|---|---|
| `Analyze → Compute SPC Status` | none — runs on any df with intensity columns |
| `Visualize → SPC Dashboard...` | `proteomics.spc_status` set OR `proteomics.campaign_index = 'true'` (dashboard renders single-run view OR multi-run trend) |
| `Campaign → Save Current Run to Campaign...` | `proteomics.de_complete = 'true'` (a run is only useful after DE) |
| `Campaign → Open Campaign...` | none — picks from AppData |
| `Campaign → Compare Runs...` | none — picks 2 runs from a campaign |
| `Share → Share Analysis for Review...` | `proteomics.de_complete = 'true'` (no point publishing a half-run analysis) |
| `Share → Open Published Analysis...` | none — picks from project list |

---

## 7. DataFrame Proliferation Model

v1.3 had two df shapes (protein-level + sample-level PCA). v1.4 introduces four more. To keep the mental model intact, classify every df by **role** (what it represents) and **lifecycle** (where it lives).

```
                       v1.4 DataFrame Topology

  ┌──────────────────────────────────────────────────────────────┐
  │                       MAIN PROTEIN df                         │
  │  rows = proteins, semType = Proteomics-ProteinId              │
  │  tags = proteomics.source, .groups, .normalized, .imputed,    │
  │         .de_complete, .de_method, [.campaign, .campaign_id,   │
  │         .campaign_run_id, .campaign_compound, .spc_status,    │
  │         .spc_metrics]                                         │
  │  lifecycle: session-scoped (in-memory)                        │
  └──────┬──────────────────┬────────────────────┬───────────────┘
         │                  │                    │
         │ (existing v1.3)  │ (NEW v1.4)         │ (NEW v1.4)
         │ clone+transpose  │ clone+strip        │ snapshot to AppData
         ▼                  ▼                    ▼
  ┌──────────────┐   ┌───────────────────┐   ┌─────────────────────┐
  │ SAMPLE df    │   │ TRIMMED REVIEW df │   │  PERSISTED RUN df   │
  │ (PCA)        │   │ (publish)         │   │ (campaign)          │
  │ rows=samples │   │ rows=proteins     │   │ rows=proteins       │
  │ cols=PC1,PC2,│   │ cols=ID,gene,FC,p,│   │ cols=full DE result │
  │   Group,name │   │   adj.p,sig       │   │  (no intensities)   │
  │              │   │ tags +=published,*│   │ tags = campaign,*   │
  │ lifecycle:   │   │ lifecycle: in     │   │ lifecycle: AppData  │
  │ session      │   │   DG.Project      │   │   CSV file          │
  └──────────────┘   └───────────────────┘   └──────────┬──────────┘
                                                         │
                                                         │ aggregated by
                                                         │ campaign loader
                                                         ▼
                                              ┌─────────────────────┐
                                              │ CAMPAIGN INDEX df   │
                                              │ rows = runs (weeks) │
                                              │ cols = run_id, week,│
                                              │  compound_vid,      │
                                              │  compound_smiles,   │
                                              │  spc_status, top_   │
                                              │  protein, impact    │
                                              │  score, n_signif    │
                                              │ tags = campaign_    │
                                              │   index = 'true'    │
                                              │ lifecycle: rebuilt  │
                                              │  on Open Campaign;  │
                                              │  not persisted      │
                                              └──────────┬──────────┘
                                                         │
                                                         │ when user picks
                                                         │ "Compare Runs"
                                                         ▼
                                              ┌─────────────────────┐
                                              │ COMPARISON DIFF df  │
                                              │ rows = proteins     │
                                              │   (union A ∪ B)     │
                                              │ cols = proteinId,   │
                                              │  gene, log2FC_A,    │
                                              │  log2FC_B, delta_   │
                                              │  log2FC, sig_in_A,  │
                                              │  sig_in_B, sig_both │
                                              │ lifecycle: session  │
                                              └─────────────────────┘
```

**Lifecycle taxonomy:**

| Lifecycle | Examples | Where it lives | When it dies |
|---|---|---|---|
| **Session-scoped** | main protein df, sample PCA df, comparison diff df | In-memory only | View closes / page refresh |
| **In a DG.Project** | trimmed review df | Serialized to platform `Project` entity via `grok.dapi.projects.save` | Project is deleted; otherwise persists across sessions |
| **AppData CSV file** | persisted run df, campaign JSON, campaign-attached layout | `System:AppData/Proteomics/campaigns/<campaign_id>/{campaign.json, runs/<run_id>.csv}` | Manual delete via Files browser |
| **Rebuilt on open** | campaign index df | Constructed from sibling AppData files when user opens a campaign | View closes |

**Sibling-DataFrame rule (extending the v1.3 PCA rule):**
> A DataFrame whose row-axis differs from the main protein df **must** be opened in its own `TableView` (`grok.shell.addTableView(...)`). Viewers built on it must dock into that view, not the main view.

Applies to: PCA (v1.3), Campaign-Index df, Comparison Diff df. The trimmed Review df is also a sibling but it lives inside a Project, so it gets its own TableView automatically when the reviewer opens the Project.

---

## 8. Cross-Package Integration Patterns

### 8.1 The two established patterns in this codebase

Both patterns are already in use in the Proteomics package (the Dendrogram lookup in `heatmap.ts`) and across the monorepo. v1.4 should match these — no new patterns.

**Pattern A — Fire-and-forget side effect** (HitTriage → Chem, line 697 of `hit-design-app.ts`):

```ts
grok.functions.call('Chem:editMoleculeCell', {cell: view.grid.cell(this._molColName, 0)});
```

When to use: when the *only* failure mode you care about is "package not installed" and that's recoverable by skipping. No return value. Fire it off, move on.

**Pattern B — Optional dependency with arity/existence guard** (Bio → Chem, `match-molecules.ts:272`):

```ts
const converterFunc = DG.Func.find({package: 'Chem', name: 'convertMoleculeNotation'})[0];
if (converterFunc && converterFunc.inputs.length === 2) {
  const result = await converterFunc.prepare({mol: smiles, targetNotation: 'molfile'}).call();
  // use result
} else {
  // graceful fallback or skip
}
```

When to use: when you need a return value, or when the call signature might drift across Chem versions. Mirrors the Proteomics → Dendrogram pattern already in `viewers/heatmap.ts:138-145`.

### 8.2 Chem integrations needed for v1.4

| Use case | Pattern | Function | Fallback |
|---|---|---|---|
| Render compound SMILES in campaign-index df + comparison diff df | Set `col.semType = DG.SEMTYPE.MOLECULE` on the SMILES column (no function call needed — Chem's cell renderer auto-fires on this semType) | n/a | If Chem not installed, cell shows SMILES text — degraded but readable |
| Validate SMILES at compound-attach time in `showNewCampaignRunDialog` | Pattern B | `DG.Func.find({package: 'Chem', name: 'validateMolecule'})[0]` | Skip validation; let the user enter any SMILES string and tolerate failures at render time |
| MCS / substructure highlight in `Compare Runs` view (chemist-oriented "what changed structurally") | Pattern B | `DG.Func.find({package: 'Chem', name: 'mcsSearch'})[0]` or similar | Hide the "Highlight Common Scaffold" button if not available |
| Edit/sketch a compound during `showNewCampaignRunDialog` | Pattern A | `grok.functions.call('Chem:editMoleculeCell', {cell})` (already used by HitTriage) | If Chem not installed, surface a text input fallback for SMILES |

**Anti-pattern to avoid:** importing from `@datagrok/chem` directly as an npm dependency. The cross-monorepo packages don't do this (Biologics' `package.json` shows `@datagrok-libraries/chem-meta` for shared *types* only; never `@datagrok/chem` as a value-level dependency). All runtime Chem coupling is via `DG.Func.find` or `grok.functions.call`.

### 8.3 Reference-package summary

| Package | What to study | Files |
|---|---|---|
| **HitTriage** | Campaign persistence model (AppData), per-campaign permissions, target/template entity, cross-package Chem editing | `src/app/hit-triage-app.ts:saveCampaign` (line ~370), `src/app/dialogs/permissions-dialog.ts`, `src/app/dialogs/save-campaign-dialog.ts`, `src/package.ts:saveCampaignJson` (line 38), `src/app/utils.ts:loadCampaigns` |
| **Bio (projects tests)** | The minimal DG.Project save+open roundtrip with tableInfo + layoutInfo + uploadDataFrame | `src/tests/projects-tests.ts:saveAndOpenProject` (line 26) |
| **Bio (Chem lookup)** | The arity-guarded cross-package call pattern | `src/utils/monomer-lib/monomer-manager/match-molecules.ts:272`, `src/utils/seq-helper/seq-helper.ts:93` |
| **Cruddy** | Generic CRUD-app framework — relevant for the "navigate campaigns by target" UX, but probably *not* worth pulling in for v1.4 (overhead > benefit for the small entity set; AppData JSON files like HitTriage are sufficient) | `libraries/cruddy/src/cruddy.ts` |
| **Peptides** | Per-compound analysis viewer composition (mutation cliffs, SAR) — useful pattern for cross-run "compare two compounds' impact" view but different domain math; useful for UX inspiration, not direct reuse | `src/viewers/`, `src/utils/algorithms.ts` |

### 8.4 Patterns v1.4 should NOT invent

There are two things that look like they need new patterns but don't:

1. **Saving a Project for read-only review.** The Bio test `saveAndOpenProject` shows the complete recipe in 15 lines. No new pattern; just compose `DG.Project.create() + project.addChild(tableInfo) + project.addChild(layoutInfo) + grok.dapi.tables.uploadDataFrame + grok.dapi.tables.save + grok.dapi.views.save + grok.dapi.projects.save`. Apply permissions via `project.permissions` (an entity-level permission object — set before saving).

2. **Storing per-target / per-campaign artifacts.** HitTriage uses package AppData (`System:AppData/HitTriage/<appName>/campaigns/<campaignId>/`) — a directory tree of JSON + CSV files. Proteomics should match: `System:AppData/Proteomics/campaigns/<campaignId>/{campaign.json, runs/<runId>.csv, layout.json}`. No databases, no entity tables. Sufficient for the scale (campaigns in the tens, runs per campaign in the tens).

---

## 9. Reviewer-Side Affordances (Detail)

When a reviewer opens a published Project, the package needs to *not* offer pipeline-editing menu items. Implementation:

- All menu handlers for `Analyze | *` items check `isPublished(df)` at the top and short-circuit with `grok.shell.warning('This analysis is published read-only and cannot be modified.')`.
- Volcano works normally (no mutation).
- A new context panel `published-analysis-panel.ts` (in `src/panels/`) fires on any df with `proteomics.published = 'true'` and renders:
  - DE method, FC threshold, p threshold (from `proteomics.published_audit`)
  - Group names (decoded from `proteomics.groups`)
  - Target, share date, sharer's friendly name
  - "Request re-run with different parameters" button (deferred behavior; for v1.4 it just copies a contact link)

The panel is registered with `@grok.decorators.panel` and a semType filter on `Proteomics-ProteinId`. Since "published" is a df-level tag, the panel function's first line checks `isPublished(df)` and returns an empty widget if not — same pattern v1.3 panels use for their preconditions.

---

## 10. Suggested Phase Sequence (Build Order)

**Constraint:** v1.4 phases continue v1.3's numbering, starting at Phase 15.

```
Phase 15 — Publishing Foundation         (lowest risk, highest leverage)
   │
   ├── Phase 16 — SPC Tracking            (independent; can be parallelized)
   │       │
   └────── │
           │
   Phase 17 — Campaign Data Model         (depends on Phase 16 SPC; produces the run df shape)
           │
   Phase 18 — Cross-Run Comparison        (depends on Phase 17; needs ≥1 campaign run schema)
```

### Phase 15 — Publishing Foundation
**Why first:** Compose-existing-artifacts work. Adds no new math, no new data model — just a trim + clone + DG.Project + permissions dialog. Lowest implementation risk. Validates the new `src/publishing/` directory pattern and the `proteomics.published*` tag family before tags get more complex in Phases 16-17. Also: highest demo value for the Cytokinetics check-in, since "I can show my boss the volcano without giving them my Excel" is the most-asked feature in the v1.3 todos.

**Deliverables:** `src/publishing/*` (4 files), `Share` menu branch (2 items), `published-analysis-panel.ts`, tests in `src/tests/publishing.ts`.

**Dependency on platform:** `grok.dapi.projects`, `grok.dapi.tables.uploadDataFrame`, `project.permissions`. All v1.3-stable APIs (used in `Bio/src/tests/projects-tests.ts`).

### Phase 16 — SPC Tracking
**Why second (parallel-with-15 possible):** Adds Shewhart + Nelson rule code and a new viewer, but operates on a *single df* via the same v1.3 pattern (mutate + tag + render). No file I/O, no new entity. The SPC dashboard can render against a single df (no campaign needed) — useful immediately for QC of a single analysis. Phase 17 *consumes* the per-run SPC status, so Phase 16 must land before Phase 17, but Phase 16 can run concurrently with Phase 15.

**Deliverables:** `src/analysis/spc.ts`, `src/viewers/spc-dashboard.ts`, `Analyze → Compute SPC Status` + `Visualize → SPC Dashboard…` menu items, tests in `src/tests/spc.ts`.

**Dependency:** None internal beyond v1.3. SPC math is pure TypeScript — no R fallback needed.

### Phase 17 — Campaign Data Model
**Why third:** Introduces the largest mental-model shift: file I/O, multi-df coordination, new top-level menu branch, new SEMTYPE values, AppData persistence. Depends on Phase 16 because each saved run records its `proteomics.spc_status` — the campaign-index df needs that column to flag bad runs at index time. Phase 17 also depends on Phase 15 if we want "Save Run + Share for Review" composability, but that composability can defer to v1.5 — Phase 17 only needs to *not block* future composability.

**Deliverables:** `src/campaigns/*` (six files), `Campaign` menu branch (3 items — `New`, `Save to`, `Open`), AppData write/read, tests in `src/tests/campaigns.ts`. Defer `Compare Runs` to Phase 18.

**Dependency:** Phase 16 SPC for `proteomics.spc_status`; cross-package Chem for SMILES rendering (Pattern A, fire-and-forget).

### Phase 18 — Cross-Run Comparison
**Why last:** Needs the campaign run schema from Phase 17 to exist with at least 2 runs persisted. Builds the Comparison Diff df + side-by-side volcano viewer + (optional) MCS substructure highlight. Highest UX polish surface — gets the spotlight in the milestone demo.

**Deliverables:** `src/viewers/campaign-comparison.ts`, `Campaign → Compare Runs…` menu item, `src/campaigns/impact-scoring.ts`, cross-DF selection wiring (same pattern as v1.2 enrichment → volcano). Tests for diff df construction + impact scoring; visual tests deferred to UAT.

**Dependency:** Phase 17 campaign storage; Chem for MCS (Pattern B with fallback).

### Why this order minimizes integration risk

| Risk | Mitigation by ordering |
|---|---|
| Tag-collision with v1.3 | Phase 15 introduces only `proteomics.published*` (5 new tags, all under one prefix) — easy to grep + revert if collision found. Phase 16 adds 3 more. Phase 17 adds 4 more on top of a stable base. |
| Cross-package coupling breakage | Phases 15-16 have zero cross-package coupling. Phase 17 adds Chem (low-risk Pattern A only — SMILES rendering via semType). Phase 18 adds Chem MCS (Pattern B with explicit fallback). |
| AppData migration during the milestone | Phase 17 is the first to touch AppData. Schema is locked at Phase 17 entry — Phase 18 is purely additive (read-only sibling files). |
| Reviewer panel breaking when df is not published | The `isPublished(df)` check is added in Phase 15. All later phases' menu handlers gate on it consistently. |
| Build order ↔ test order | Each phase ships its own `src/tests/<phase>.ts` file. Test bundle stays small per phase. |

---

## 11. Anti-Patterns Specific to v1.4

### Anti-Pattern 1: Letting "published" + "campaign" tags coexist on the same df
**What people will do:** Save a Project containing a campaign-run df, accidentally setting both `proteomics.published` and `proteomics.campaign` on the same DataFrame. The reviewer-side panel and the campaign-comparison view both fire.
**Why wrong:** Reviewer panel offers "Request re-run" (which only makes sense for a one-shot analysis); comparison view tries to find sibling runs in a Project that doesn't have any.
**Do instead:** `markPublished(df, ...)` calls `df.deleteTag('proteomics.campaign*')` first. The trim function builds a fresh clone so this is naturally enforced — but document the invariant so future contributors don't break it.

### Anti-Pattern 2: Storing audit context as separate tags instead of one JSON tag
**What people will do:** Add `proteomics.published_fc_threshold`, `proteomics.published_p_threshold`, `proteomics.published_de_method`, etc., one tag per field.
**Why wrong:** Tag namespace pollution; the platform stores tags as a flat string map; audit fields will grow. v1.3 already did this with the smart-filter telemetry (5 tags under `proteomics.enrichment_smart_filtered*`) — don't repeat.
**Do instead:** One `proteomics.published_audit` tag, JSON-encoded. Mirror the v1.3 `proteomics.groups` pattern.

### Anti-Pattern 3: Querying AppData from inside a viewer
**What people will do:** Open the campaign-comparison view → it calls `grok.dapi.files.readCsv` to load the sibling run df → blocks the UI thread.
**Why wrong:** Viewers are render-only. v1.3 invariant: viewers consume an already-loaded df; they don't fetch.
**Do instead:** The `Compare Runs` menu handler in `package.ts` loads both run dfs from AppData *first*, builds the diff df, *then* calls `createSideBySideVolcano(diffDf, runA, runB)`. Viewer stays synchronous.

### Anti-Pattern 4: Tying the SPC rules to the campaign concept
**What people will do:** Put SPC math under `src/campaigns/` because "SPC is for campaigns."
**Why wrong:** SPC is just QC. A scientist running a single one-off analysis benefits from "compute SPC status on this df." Forcing campaign context is hostile.
**Do instead:** Phase 16 SPC math lives in `src/analysis/spc.ts` and works on any df with intensity columns. Phase 17 simply *consumes* the SPC status from the run df's tags when building the campaign index.

### Anti-Pattern 5: Persisting layout JSON inside the campaign-index df
**What people will do:** Add a `layout` column to the campaign-index df.
**Why wrong:** Layouts are per-view artifacts, not per-row data. HitTriage stores layout as a `campaign.layout` field in `campaign.json` (`hit-triage-app.ts:425`).
**Do instead:** Match HitTriage — store run-specific layout in `runs/<runId>.layout.json`; store campaign-level layout in `campaign.json`. Index df stays slim.

---

## 12. Integration Points

### Datagrok platform APIs touched by v1.4 (none new — all v1.3-stable)

| API | Used by | Pattern |
|---|---|---|
| `grok.dapi.projects.save/find/open/delete` | Phase 15 publish/load | Direct call, awaited |
| `grok.dapi.tables.uploadDataFrame` + `grok.dapi.tables.save` | Phase 15 publish | Sequential — upload, then save tableInfo |
| `grok.dapi.views.save` | Phase 15 publish (save volcano layout) | After view exists in shell |
| `grok.dapi.groups.find` (by id) | Phase 15 permissions dialog | Per-group-id lookup; cache in dialog state |
| `grok.dapi.users.find` (by id) | Phase 15 reviewer panel | Resolve author friendly name |
| `_package.files.readAsText` / `writeAsText` | Phase 17 campaign AppData I/O | Async, awaited; mirror HitTriage |
| `grok.dapi.files.readCsv` / `writeAsText` | Phase 17 campaign run data | Platform-level path-based variant for `System:AppData/...` paths |
| `grok.functions.call('Chem:editMoleculeCell', ...)` | Phase 17 compound dialog | Pattern A — fire-and-forget |
| `DG.Func.find({package: 'Chem', name: 'validateMolecule'})` | Phase 17 SMILES validation | Pattern B with `?.[0]` and arity check |
| `DG.Func.find({package: 'Chem', name: 'mcsSearch'})` | Phase 18 substructure highlight | Pattern B with fallback to "skip MCS" |
| `df.setTag` / `df.getTag` / `df.deleteTag` | All phases | Via setter helpers, not raw calls |

### Internal package boundaries

| Boundary | Communication | Notes |
|---|---|---|
| `publishing/` → `analysis/` | Reads tags + `getGroups(df)`; never calls dialogs | One-way: publishing depends on analysis; not the reverse |
| `publishing/` → `viewers/` | Calls `createVolcanoPlot(trimmedDf, opts)` to build the embedded layout | Publishing owns the trimmed df + Project; viewers don't know about publishing |
| `campaigns/` → `analysis/spc.ts` | Reads `getSpcStatus(runDf)` when building the index | One-way; SPC has no campaign awareness |
| `campaigns/` → `viewers/campaign-comparison.ts` | Passes pre-built diff df + run df pair | Comparison viewer is pure render — no AppData I/O |
| `campaigns/` → `publishing/` | None in v1.4 | Composability ("publish a campaign") deferred to v1.5 |
| `panels/published-analysis-panel.ts` → `publishing/publish-state.ts` | Calls `getPublishedAudit(df)` | Panel function reads tags; doesn't import dialog code |

---

## 13. Open Questions for Phase Planning

These are decisions deferred to the per-phase plans, not blockers for the roadmap:

1. **Target identifier shape** (Phase 15): free-text label vs. UniProt accession vs. first-class Datagrok Space? Existing todo `2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md` suggests Spaces; investigation needed on whether `grok.dapi.spaces` exists / is mature. Fallback: free-text + filter on `proteomics.published_target` in the project list.
2. **SPC control-period definition** (Phase 16): is the "in-control" mean/stddev computed from the first N runs of the campaign, a rolling window, or user-specified? Default proposal: rolling 8-run window once campaign has ≥4 runs; full-history before that.
3. **Compound identity uniqueness** (Phase 17): can a single campaign run two analyses of the same compound? Likely yes (technical replicates). Compound `vid` is not a unique key — runId is. Index df has `compound_vid` + `run_id` and may show multiple rows per compound.
4. **Comparison row alignment** (Phase 18): outer join (union of proteins) vs. inner join (intersection)? Proposal: outer join with `null` log2FC entries for the missing side; the viewer's NS-cloud color naturally suppresses these visually.
5. **Permission scope of published projects** (Phase 15): per-project ACL via `project.permissions` (Bio test pattern) vs. Datagrok Space ACL? Bio tests do per-project. Sufficient for v1.4 scope.
6. **Layout snapshot inside the Project** (Phase 15): just the volcano, or the volcano + heatmap + enrichment-table all docked? Smaller surface = lower risk. Start with volcano-only; add others in v1.5 polish.

---

## Sources

- Existing package codebase audit:
  - `.planning/codebase/ARCHITECTURE.md` (2026-05-11)
  - `.planning/codebase/STRUCTURE.md` (2026-05-11)
  - `.planning/codebase/INTEGRATIONS.md` (2026-05-11)
- Backlog todos with detailed scoping notes:
  - `.planning/todos/pending/2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md` (Phase 999.3 backlog)
  - `.planning/todos/pending/2026-05-11-compare-most-impactful-compounds-across-weekly-screening-runs.md` (Phase 999 backlog)
- Reference packages inspected directly in this worktree:
  - `packages/HitTriage/src/app/hit-triage-app.ts` (lines 370-430: `saveCampaign`)
  - `packages/HitTriage/src/package.ts` (line 38: `saveCampaignJson`)
  - `packages/HitTriage/src/app/dialogs/permissions-dialog.ts` (full file)
  - `packages/HitTriage/src/app/dialogs/save-campaign-dialog.ts` (full file)
  - `packages/HitTriage/src/app/types.ts` (Campaign + TriagePermissions interfaces)
  - `packages/Bio/src/tests/projects-tests.ts` (lines 26-47: `saveAndOpenProject` — minimal Project roundtrip)
  - `packages/Bio/src/utils/monomer-lib/monomer-manager/match-molecules.ts:272` (Pattern B Chem lookup)
  - `packages/Proteomics/src/package.ts` (current menu structure)
  - `packages/Proteomics/src/analysis/experiment-setup.ts` (`getGroups`/`setGroups` setter pattern)
  - `packages/Proteomics/src/viewers/heatmap.ts:138-145` (Dendrogram cross-package lookup — Pattern B precedent)
- Tag-collision verification: grep over `setTag/getTag` calls in `packages/Proteomics/src/`.

---

*Architecture research for: v1.4 Cross-Team Review integration into v1.3 Proteomics package*
*Researched: 2026-06-06*
