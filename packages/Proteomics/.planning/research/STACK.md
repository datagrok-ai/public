# Stack Research — v1.4 Cross-Team Review

**Project:** Proteomics Package for Datagrok
**Researched:** 2026-06-06
**Scope:** ADDITIONS for v1.4 only — read-only Project publishing, sample-level SPC, cross-run screening comparison. The TypeScript + webpack + Datagrok-API + R-via-`grok.functions.call` baseline (STACK.md v1.3) is unchanged.
**Confidence:** HIGH

## Headline

**Zero new npm dependencies.** All three v1.4 capabilities ship on platform-native primitives already pulled in by v1.3, plus three new in-package modules (publishing, SPC, campaign model). The single integration point with another Datagrok package is `Chem` — used via the cross-package `grok.functions.call('Chem:...')` and `Chem:chemCellRenderer` registration mechanism, not as an npm dep.

| Capability | Headline tech choice | Why |
|---|---|---|
| Read-only publishing | `DG.Project` + `addChild(tableInfo)` + `grok.dapi.projects.save` + `grok.dapi.permissions.grant` | Native, durable, opens via standard UI, permission API already supports group-scoped view/edit, no schema needed. |
| Cross-DataFrame "campaign" | DataFrame tags (`proteomics.campaign_id`) + a campaign-index DataFrame stored in `_package.files` | Matches the existing pipeline-state model (everything is a `proteomics.*` tag); avoids over-engineering. HitTriage precedent uses the same shape. |
| SPC rules engine | **Build in-package** (~390 LOC TS for Shewhart I-MR + Nelson rules 1–8 + run-QC extractor + viewer + store) | No usable npm; `@datagrok-libraries/statistics` doesn't ship SPC; math is small, well-published, and trivially unit-testable. |
| Compound structure rendering | `grok.functions.call('Chem:detectSmiles', ...)` + `col.meta.cellRenderer = 'Molecule'` + `SEMTYPE.MOLECULE` | Chem auto-registers a `'Molecule'` cell renderer; setting the semType on the column wires it up. Zero coupling. |

---

## Recommended Stack — Additions Only

### Datagrok JS-API surface (already pulled in via `datagrok-api ^1.25.0`)

| API | Method/property | Use in v1.4 | Confirmed in current public js-api |
|---|---|---|---|
| `DG.Project` | `Project.create()` (static) | Construct a new project | `js-api/src/entities/project.ts:31` |
| `DG.Project` | `addChild(tableInfo \| dataframe)` | Add a TableInfo / DataFrame to the project | `entities/project.ts:101–105` (overload accepts `DataFrame`; auto-resolves to `getTableInfo()`) |
| `DG.Project` | `description`, `options` (MapProxy) | Capture publish-time metadata (target, run date, source DF id, group filter) | `entities/project.ts:26, 50` |
| `grok.dapi.tables` | `uploadDataFrame(df)` then `save(tableInfo)` | Required before `projects.save` to persist trimmed DataFrame | `dapi.ts:990` |
| `grok.dapi.projects` | `save(project, {saveRelations: true})` | Persist project + child tables in one round-trip; `saveRelations` defaults true | `dapi.ts:794–799` |
| `grok.dapi.permissions` | `grant(entity, group, edit:boolean)` | Reviewer-group access. `edit=false` → view-only | `dapi.ts:700` |
| `grok.dapi.permissions` | `get(entity)` / `check(entity, 'View'\|'Edit'\|'Share'\|'Delete')` | Read-back of grants; UAT verification | `dapi.ts:682, 693` |
| `grok.dapi.groups` | `filter('name = "X"').first()` / `createNew(name)` | Resolve or create a reviewer group | `dapi.ts:450` + `HttpDataSource.filter` |
| `grok.dapi.groups` | `addMember(g, m)`, `addAdminMember(g, m)` | Populate reviewer group from a curator dialog | `dapi.ts:462, 469` |
| `DG.DataFrame` | `clone(bitset, columns)` + `columns.remove(name)` | "Trimmed" projection: keep protein-id / gene / log2FC / p / adj.p / sig (+ enrichment if present); drop intensities and per-sample columns | DataFrame API (used in v1.3 heatmap clone pattern) |
| `DG.DataFrame.tags` | `df.setTag('proteomics.X', v)` | Mark the published table as `proteomics.published=true`, carry `proteomics.target`, `proteomics.de_method`, etc. into the project copy | Used pervasively in v1.3 |
| `grok.dapi.fetchProxy` | already used for UniProt/g:Profiler | Not new, but the only sanctioned external-URL primitive — keep using it if SPC adds any external lookups (none anticipated) | `dapi.ts:193` |

**Canonical publish flow (matches `packages/ApiSamples/scripts/dapi/projects.js`):**

```ts
const trimmed = df.clone(df.filter, ['Protein.Id', 'Gene', 'log2FC', 'p-value', 'adj.p-value', 'significant']);
trimmed.name = `${target} — ${runDate} — ${df.name}`;
trimmed.setTag('proteomics.published', 'true');
trimmed.setTag('proteomics.target', target);

const project = DG.Project.create();
project.description = `Published proteomics analysis for ${target}. Source DF: ${df.name}. Method: ${df.getTag('proteomics.de_method')}.`;
project.addChild(trimmed.getTableInfo());

await grok.dapi.tables.uploadDataFrame(trimmed);
await grok.dapi.tables.save(trimmed.getTableInfo());
await grok.dapi.projects.save(project);

const reviewers = await grok.dapi.groups.filter(`name = "${reviewerGroupName}"`).first()
  ?? await grok.dapi.groups.createNew(reviewerGroupName);
await grok.dapi.permissions.grant(project, reviewers, /* edit */ false);
```

Verified against `packages/ApiSamples/scripts/dapi/projects.js` (3-line canonical form) and `.../layouts-and-permissions.js` (permissions.grant). This is the Datagrok-native shape — not a new pattern.

**Project vs Space:** the milestone context already declared **Project** as the chosen primitive. Confirmed correct: `Spaces` (`grok.dapi.spaces`, `dapi.ts:805`) are hierarchical containers with their own storage, designed for organizing collections of entities — overkill for "one frozen published analysis = one project". Projects open as table views via standard UI flow (`project.open()`), Spaces don't.

### Cross-DataFrame "campaign" data model

**Pattern (confirmed Datagrok-native):** DataFrame tags + a package-owned **campaign-index DataFrame** persisted under `_package.files`. No Postgres schema, no `@datagrok-libraries/cruddy`, no custom entity type.

| Layer | Shape | Where |
|---|---|---|
| Per-run tag | `df.setTag('proteomics.campaign_id', <uuid>)` + `proteomics.run_id` + `proteomics.run_date` + `proteomics.compound` (SMILES or registry id) | On every analysis DataFrame |
| Per-run index row | A row in `campaign-<id>.csv` stored in `System:AppData/Proteomics/campaigns/<id>/index.csv`: campaign_id, run_id, run_date, compound, source_df_id, de_method, published_project_id, qc_status | `_package.files.writeAsText(...)` (`dapi.ts:1465`) |
| Campaign loader | `loadCampaign(id) → DG.DataFrame[]` — reads index, fetches each saved TableInfo via `grok.dapi.tables.getTable(id)` | New `src/campaigns/campaign-store.ts` |

**Why this and not cruddy / a new schema:**
- The existing v1.3 pipeline already stores **all workflow state in `proteomics.*` tags** (see `packages/Proteomics/CLAUDE.md` → "Tags"). Introducing a parallel state store would break the single-source-of-truth contract.
- `@datagrok-libraries/cruddy` is for "CRUD applications on top of Datagrok" — full forms-over-data UIs (HitDesign, MolTrack style). v1.4's campaign object has ≤6 fields and ≤O(100) rows per project; cruddy is several orders of magnitude over-spec.
- HitTriage uses the same files-in-AppData pattern: `System:AppData/HitTriage/{appName}/campaigns/{campaignId}/` (see `packages/HitTriage/src/app/hit-triage-app.ts:375`). Different domain, same shape — a precedent.

**When to escalate to a Postgres schema (NOT v1.4):** if campaigns exceed O(1000) rows per project, need server-side filtering, or need to be shared across multiple deployments. At that point the model becomes the same as `packages/Plates/databases/plts/0000_init.sql` — package-owned schema in `databases/proteomics/`. Park it for v1.5+.

### Chem package integration (compound structure rendering)

Chem exposes its full surface via `grok.functions.call('Chem:<funcName>', {...})` (this is the universal cross-package call shape — see `tools/GROK_S.md`). v1.4 needs only the renderer registration + (optionally) SMILES normalization. The functions below were grepped directly from `packages/Chem/src/package.ts`.

| Chem function | Signature | v1.4 use |
|---|---|---|
| `Chem:detectSmiles(col, min)` | `(col: DG.Column, min: int) → void`; if detected, sets `col.semType = 'Molecule'`, `col.meta.units = 'smiles'`, `col.meta.cellRenderer = 'Molecule'` | Call after parsing each campaign run's `compound` column — auto-wires rendering platform-wide (`packages/Chem/src/package.ts:2106`) |
| `Chem:convertMolNotation(molecules, targetNotation, kekulize?)` | `(DG.Column, 'smiles'\|'molblock'\|'inchi', boolean?) → Promise<DG.Column>` | Normalize compound IDs across runs (different vendors may export different notations) (`packages/Chem/src/package.ts:1690, 1707`) |
| `Chem:isSmiles(s)` | `(string) → boolean` | Defensive check before setting semType when compound column is heterogeneous (`packages/Chem/src/package.ts:2090`) |
| `Chem:getProperties(molecules, selected?)` | `(DG.Column, string[]?) → Promise<DG.DataFrame>` | Optional — adds MW/cLogP/PSA columns if scientists want structure-property side-by-side in cross-run comparison (`packages/Chem/src/package.ts:2180`) |
| `Chem:chemCellRenderer` (registered, not called) | meta `{cellType: 'Molecule', role: 'cellRenderer'}` | Auto-resolves whenever any column has `semType=Molecule`. We don't call it — we just set the semType (`packages/Chem/src/package.ts:411–429`) |

**Minimum integration shape (the one-liner that does the job):**

```ts
// In src/campaigns/parse-compound.ts after a per-run parser identifies the compound column
import * as grok from 'datagrok-api/grok';
await grok.functions.call('Chem:detectSmiles', {col: compoundCol, min: 5});
// That's it. Datagrok now renders structures everywhere this column appears.
```

**Side-by-side volcano comparison:** consume the existing v1.3 `createVolcanoPlot` factory unchanged. The comparison viewer is a thin orchestrator that opens two `createVolcanoPlot` calls in a 1-row docked layout and links their selections via the **same `onCurrentRowChanged` / `df.selection` pattern** already used in `viewers/enrichment-viewers.ts` (the v1.2 cross-DF link). No new viewer infrastructure.

### SPC (Shewhart I-MR + Nelson rules) — Build, don't buy

**Build verdict: HIGH confidence.**

**What's in `@datagrok-libraries/statistics` today** (grep `libraries/statistics/src/*.ts`):

| File | Public surface | SPC-relevant? |
|---|---|---|
| `tests.ts` | `tTest`, `uTest` | No |
| `multiple-tests.ts` | `fdrcorrection` | No |
| `box-plot-statistics.ts` | `calculateBoxPlotStatistics` (min/Q1/median/Q3/max) | Marginally — gives quartiles, but not σ̂/d2/d4 constants or rules engine |
| `confidence-intervals.ts` | `getConfidence` | No |
| `correlation-coefficient.ts` | `kendallsTau` | No |
| `fit/` | dose-response fits | No |
| `mpo/` | multi-property optimization | No |
| `compute-functions/` | misc compute | No |

**No Shewhart, no Nelson, no Western Electric, no moving range, no I-MR constants, no rules engine.** Grep across all `libraries/` and `packages/` confirms.

**npm landscape (verified 2026-06-06):**

| Package | Latest | License | Types | Implements | Verdict |
|---|---|---|---|---|---|
| `process-control-charts` | 1.0.6 | MIT | None | Attribute charts only (p/np/c/u). Wrong type — we need continuous (intensity median, missingness %). | **Reject** |
| `spc-limits-js` (GitHub: AUS-DOH-Safety-and-Quality) | repo not resolvable via GitHub API (likely renamed/deleted) | unknown | unknown | Limits only (not rules engine) per advertised description | **Reject — non-existent / untrusted** |
| `qcc` (R) | already noted in `packages/Minitab/docs/MINITAB_STATISTICAL_METHODS.md:17, 537` | GPL | n/a | R-only; can be called via existing R-script path | **Optional fallback** — not the primary path |

**LOC estimate for in-package build:**

| Module | Function | LOC est |
|---|---|---|
| `src/spc/shewhart.ts` | `computeIndividualsAndMovingRange(values: number[]) → {center, ucl, lcl, sigmaHat, mr[]}` (constants d2=1.128 for n=2) | ~40 |
| `src/spc/nelson-rules.ts` | `applyNelsonRules(values, center, sigma) → Array<{ruleId, idxs}>` for rules 1–8 | ~150 |
| `src/spc/run-qc-metrics.ts` | `extractRunQcMetrics(df, groups) → {runId, intensityMedian, missingPct, ctrlCorr, proteinCount}` | ~80 |
| `src/spc/spc-store.ts` | Index-DataFrame persistence (read/write to `AppData/Proteomics/spc/spc-history.csv`) | ~50 |
| `src/viewers/spc-chart.ts` | I-MR viewer factory (delegates to `DG.Viewer.lineChart` with formula lines for UCL/CL/LCL + flagged-point coloring) | ~70 |
| **Total** | | **~390 LOC** |

For comparison, v1.3 shipped 7,607 src LOC. SPC is ~5% of that; well within the "small, well-known math" zone.

**References (math is in the public domain):**
- Nelson Rules 1–8 — Lloyd S. Nelson, J. Quality Technology, 1984. Reference summary at [Quality Gurus — Nelson Rules](https://www.qualitygurus.com/nelson-rules-and-western-electric-rules-for-control-charts/) and [QI Macros — Nelson Rules](https://www.qimacros.com/control-chart/nelson-rules/).
- Individuals & Moving Range constants (d2=1.128, D3=0, D4=3.267 for n=2) — Wheeler, *Understanding Variation*. See [SPC Press — Using Extra Detection Rules](https://www.spcpress.com/pdf/DJW322.Oct.17.Using%20Extra%20Detection%20Rules.pdf).
- Nelson rule patterns: 1 (point > 3σ), 2 (9-in-a-row same side), 3 (6-in-a-row trending), 4 (14-in-a-row alternating), 5 (2-of-3 > 2σ same side), 6 (4-of-5 > 1σ same side), 7 (15-in-a-row within 1σ), 8 (8-in-a-row > 1σ either side).

**Optional R fallback (not blocking):** the existing `scripts/limma_de.R` pattern can be reused for an `scripts/qcc_spc.R` that calls R's `qcc` package, mirroring the v1.3 DE three-level fallback. Recommend deferring this to v1.5 — the rules math is deterministic enough that the client-side path will give bit-identical results across runs.

### What's NOT being added

| Avoided | Why | Replacement |
|---|---|---|
| `@datagrok-libraries/cruddy` | Full CRUD-app framework. Campaign object is 6 fields, ≤O(100) rows — wrong scale. | Tag-based + AppData CSV (HitTriage precedent) |
| `@datagrok-libraries/db-explorer` | Database schema browsing UI. We have zero new DB schemas. | n/a — no DB |
| Custom Postgres schema under `databases/proteomics/` | Adds ops burden (migrations, GRANT scripts), couples v1.4 ship to a server config. Cytokinetics check-in is the milestone close — tight scope. | AppData CSV with documented "promote to DB at v1.5+ if O(1000)+ campaigns" trigger |
| `process-control-charts` npm | Attribute charts only; no rules engine; no TS types; 0 deps but also 0 fit | In-package `src/spc/` |
| `qcc` R package on the critical path | Adds R-environment dependency to QC flow. v1.3 R is for DE only, and even there has a TS fallback. | TS-first, R as optional v1.5+ fallback |
| Custom JsViewer for I-MR / volcano comparison | Datagrok's `lineChart` + formula lines (`DG.Viewer.lineChart`) handles UCL/CL/LCL natively. Volcano comparison reuses `createVolcanoPlot`. | Factory functions returning configured `DG.Viewer` |
| New top-level dependency on `@datagrok-libraries/chem-meta` | Chem package is called via `grok.functions.call('Chem:...')`. Pulling chem-meta would bundle RDKit JS into Proteomics — bloats the bundle for a feature 99% of users won't trigger. | Cross-package function calls |
| Custom `DG.Entity` subclass for Campaign | Datagrok platform doesn't let third-party packages register new entity types. The "Project containment" + tag pattern is the platform-blessed substitute. | Tag-based |

### Where new code lands (extends existing convention)

This slots into the existing v1.3 layout (see `packages/Proteomics/CLAUDE.md` → "Architecture"):

```
src/
  package.ts              # ADD: top-menu items "Proteomics | Publish for Review...",
                          #      "Proteomics | SPC Drift...", "Proteomics | Compare Runs..."
  analysis/               # unchanged
  parsers/                # unchanged
  viewers/
    volcano.ts            # reused as-is for cross-run comparison
    spc-chart.ts          # NEW — I-MR viewer
    comparison.ts         # NEW — side-by-side orchestrator (calls createVolcanoPlot twice)
  panels/                 # unchanged
  publishing/             # NEW — read-only Project export
    publish-dialog.ts     # showPublishDialog(df)
    trim-dataframe.ts     # trimForPublish(df, opts) → DG.DataFrame
    publish.ts            # publishToProject(trimmed, target, reviewers) → Promise<Project>
  campaigns/              # NEW — cross-run model
    campaign-store.ts     # load/save campaign index from AppData
    campaign-tag.ts       # tag setter/reader (mirrors getGroups/setGroups in analysis/experiment-setup)
    compound-utils.ts     # detectSmiles wrapper, optional normalization
  spc/                    # NEW — Shewhart + Nelson + per-run QC extractor
    shewhart.ts
    nelson-rules.ts
    run-qc-metrics.ts
    spc-store.ts
  utils/proteomics-types.ts  # ADD: SEMTYPE.CAMPAIGN_ID, SEMTYPE.RUN_ID, SEMTYPE.COMPOUND
                            #   (and mirror in detectors.js — per the package contract)
  tests/                  # one test file per new module (existing pattern)
```

This matches the "Function-naming convention (load-bearing)" contract documented in `packages/Proteomics/CLAUDE.md`: `showPublishDialog` (dialog), `publishToProject` (mutating action), `createSpcChart` (pure factory), `parseCampaignIndex` (parser).

### New `proteomics.*` tags (extends v1.3 tag contract)

Following the v1.3 convention that **all workflow state lives in `proteomics.*` tags**:

| Tag | Value | Set by | Read by |
|---|---|---|---|
| `proteomics.published` | `'true'` | publish flow | reviewers' view; prevents accidental re-publishing of a frozen copy |
| `proteomics.target` | string (target name, e.g. `HER2`) | publish flow + campaign-aware analysis | publish dialog default, project description, campaign index key |
| `proteomics.campaign_id` | UUID string | campaign assignment | campaign loader, comparison viewer, SPC aggregator |
| `proteomics.run_id` | string | campaign assignment | SPC x-axis, comparison labels |
| `proteomics.run_date` | ISO date | campaign assignment | SPC ordering, comparison title |
| `proteomics.compound` | SMILES or registry id | campaign assignment | Chem renderer trigger, comparison hover |
| `proteomics.qc_pass` | `'true'\|'flag'\|'fail'` | SPC engine after `applyNelsonRules` | grid coloring, publish-blocker policy |

Add these to `packages/Proteomics/CLAUDE.md` `DataFrame tags` table at v1.4 close.

---

## Installation

**None required.** All capabilities ship on the existing `package.json` dependency set:

```jsonc
// packages/Proteomics/package.json  — UNCHANGED for v1.4
{
  "dependencies": {
    "@datagrok-libraries/bio": "^5.61.3",
    "@datagrok-libraries/math": "^1.2.6",
    "@datagrok-libraries/ml": "^6.10.7",
    "@datagrok-libraries/statistics": "^1.2.12",
    "@datagrok-libraries/utils": "^4.6.9",
    "@datagrok-libraries/test": "^1.1.0",
    "cash-dom": "^8.1.5",
    "datagrok-api": "^1.25.0",
    "rxjs": "^6.5.5",
    "wu": "^2.1.0"
  }
}
```

This continues the v1.3 "zero new npm dependencies" decision (`PROJECT.md`, Key Decisions table). Bundle size: unchanged.

---

## Alternatives Considered

| Recommended | Alternative | When alternative would be better |
|---|---|---|
| `DG.Project` for publishing | Datagrok **Space** (`grok.dapi.spaces`) | If we needed multiple analyses + supplementary files + nested subspaces per published unit. v1.4 ships one analysis per publish, no extra files, no hierarchy. Space adds ceremony with zero payoff. |
| Tags + AppData CSV for campaigns | Package-owned Postgres schema (`databases/proteomics/0000_init.sql`) | When campaign rows exceed O(1000) per project, or when we need server-side filtering, or when campaigns must be shared across deployments. Document the trigger; defer to v1.5+. |
| Tags + AppData CSV for campaigns | `@datagrok-libraries/cruddy` | If campaigns grew into a full editable forms-over-data UI (HitDesign-class). v1.4 campaigns are read-mostly with curator-only edit. |
| In-package SPC math | `process-control-charts` npm | Never — wrong chart type (attribute only, we need continuous), no TS types, abandoned-looking publish history. |
| In-package SPC math | R `qcc` via `scripts/spc_qcc.R` | Acceptable as a v1.5+ optional fallback (matches v1.3 DE three-level cascade pattern). Not the primary path because (a) adds server R dep for a feature most users will not trigger; (b) deterministic math means TS and R will agree to floating-point precision. |
| `grok.functions.call('Chem:...')` for structures | Direct dep on `@datagrok-libraries/chem-meta` | Direct dep bundles RDKit JS (~MB) into Proteomics. The `grok.functions.call` route runs in the Chem package's already-loaded RDKit instance. |
| Reuse `createVolcanoPlot` for comparison | Build a new comparison-specific JsViewer | Custom JsViewer doubles maintenance for ~0 functional gain. v1.3 already shipped 7 reusable viewers; this is the established pattern. |

---

## What NOT to Use

| Avoid | Why | Use Instead |
|---|---|---|
| `process-control-charts` (npm) | Attribute charts only (p/np/c/u — for defect-rate data). We need continuous (I-MR). 0 deps but also 0 types, 0 rules engine. | In-package `src/spc/` (~390 LOC TS) |
| `spc-limits-js` (npm/GitHub) | Repo not resolvable via GitHub API; cannot verify license, maintenance, or scope. Risk of importing abandonware. | In-package `src/spc/` |
| `@datagrok-libraries/cruddy` for campaigns | Designed for full CRUD-app UIs. v1.4 campaign object has 6 fields and ≤O(100) rows. | Tags + AppData CSV |
| Direct `fetch()` for cross-DF wiring | Datagrok anti-pattern — bypasses auth/routing. Already documented in `packages/CLAUDE.md` and the project rules. | `grok.dapi.*` for entities, `grok.dapi.fetchProxy()` for external URLs |
| Custom `DG.Entity` subclass for Campaign | Platform doesn't expose entity-type registration to packages | `DG.Project.options` (MapProxy) + `proteomics.*` tags |
| Hand-rolled SMILES parser / structure renderer | Chem already ships RDKit + cell renderer | `grok.functions.call('Chem:detectSmiles', ...)` |
| `webpack.config.js` change for new modules | New `src/{publishing,campaigns,spc}/` import normally — no externals needed (no new platform-provided libs) | Leave webpack config alone |

---

## Stack Patterns by Variant

**If reviewers belong to multiple existing Datagrok groups (e.g., Cytokinetics has `Biology`, `Computational`, `Mgmt`):**
- Loop `grok.dapi.permissions.grant(project, group, false)` once per group
- This is documented and supported — `permissions.grant` is per-entity-per-group, no upper bound on calls

**If a single campaign spans >100 runs (rare for proteomics; SPC literature points are usually weekly):**
- Switch the campaign-index store from CSV-in-AppData to `databases/proteomics/0000_init.sql` (Postgres schema, following `packages/Plates/databases/plts/` precedent)
- Migration is mechanical: add a SQL schema, queries under `queries/`, swap the campaign-store load/save implementation

**If SPC needs to flag a run for "no publish" (block-on-fail policy):**
- Publish dialog reads `df.getTag('proteomics.qc_pass')` and disables the OK button + shows the flag message
- No new infrastructure needed — same gating pattern as `proteomics.de_complete` already used by volcano/heatmap

**If a deployment lacks R (which is the default assumption per v1.3):**
- Everything in v1.4 still works — none of the three capabilities depend on R
- The optional `scripts/spc_qcc.R` fallback (v1.5+) is the only R touchpoint and is opt-in

---

## Version Compatibility

| Package | Compatible With | Notes |
|---|---|---|
| `datagrok-api ^1.25.0` (current) | `Project.addChild`, `permissions.grant`, `groups.createNew`, `tables.uploadDataFrame`, `tables.save`, `fetchProxy` | All confirmed in `js-api/src/dapi.ts` and `js-api/src/entities/project.ts` at HEAD of this worktree. No version bump needed. |
| `@datagrok-libraries/statistics ^1.2.12` | `tTest`, `fdrcorrection`, `calculateBoxPlotStatistics` | Already used; no SPC primitives — verified by grep across all `*.ts` in `libraries/statistics/src/`. We don't bump this. |
| Chem package (cross-package call) | `Chem:detectSmiles`, `Chem:convertMolNotation`, `Chem:isSmiles`, `Chem:getProperties`, `Chem:chemCellRenderer` (registered) | All present in `packages/Chem/src/package.ts` at HEAD. Cross-package calls are version-decoupled via the `grok.functions.call` dispatch — Chem can upgrade independently. |

---

## Sources

**Authoritative — repo-local (HIGH confidence):**
- `js-api/src/dapi.ts` — Dapi surface (Projects, Permissions, Groups, Tables, Files, Spaces). Read at HEAD `21d37776ee`.
- `js-api/src/entities/project.ts` — `Project.create()`, `addChild`, `addLink`, `description`, `options`. Read at HEAD.
- `js-api/src/const.ts` (line 227) — `SEMTYPE.MOLECULE = 'Molecule'`.
- `packages/Chem/src/package.ts` (lines 411, 1690, 1707, 2090, 2106, 2180) — Chem function entrypoints. Read at HEAD.
- `packages/ApiSamples/scripts/dapi/projects.js` — canonical publish + addChild + save pattern (3 lines).
- `packages/ApiSamples/scripts/dapi/layouts-and-permissions.js` — canonical permissions.grant pattern.
- `packages/ApiSamples/scripts/dapi/groups.js` — group create + member management pattern.
- `packages/HitTriage/src/app/hit-triage-app.ts` (line 375) — AppData/campaigns pattern precedent.
- `packages/Plates/databases/plts/0000_init.sql` — package-owned schema pattern (for the "escalate-to-DB" trigger).
- `libraries/statistics/src/{tests,multiple-tests,box-plot-statistics,confidence-intervals,correlation-coefficient}.ts` — confirms no SPC primitives.
- `packages/Minitab/docs/MINITAB_STATISTICAL_METHODS.md` — internal Minitab plugin's own notes corroborate the "implement Nelson rules in JS" verdict.
- `packages/Proteomics/CLAUDE.md` — current pipeline tag contract, function-naming convention.

**External (MEDIUM confidence, math is in the public domain so OK):**
- [Quality Gurus — Nelson Rules and Western Electric Rules](https://www.qualitygurus.com/nelson-rules-and-western-electric-rules-for-control-charts/) — Nelson rule 1–8 reference
- [QI Macros — Nelson Rules](https://www.qimacros.com/control-chart/nelson-rules/) — corroborating reference
- [SPC Press — Using Extra Detection Rules (PDF)](https://www.spcpress.com/pdf/DJW322.Oct.17.Using%20Extra%20Detection%20Rules.pdf) — Wheeler's I-MR constants
- [Parsec — Nelson vs Western Electric](https://www.parsec-corp.com/blog/nelson-vs-western-electric) — confirms 7 of 8 rules are shared; rule 4 differs

**External (LOW confidence — used only to verify the "no usable npm" finding):**
- [npm — process-control-charts](https://www.npmjs.com/package/process-control-charts) — confirmed attribute-chart only, no types (registry direct: latest 1.0.6, MIT, description "methods to solve p-chart, np-chart, c-chart, u-chart", no dependencies, no `types` field)
- GitHub API call for `AUS-DOH-Safety-and-Quality/spc-limits-js` — returned empty; repo not currently accessible

**Gaps to surface to roadmapper:**
- The optional R `qcc` SPC fallback path is deferred (v1.5+). If milestone scope ever expands to require R, the work is "new `scripts/spc_qcc.R` mirroring `scripts/limma_de.R` pattern" — straightforward but blocks on Datagrok R env availability.
- The "publish-block on QC fail" policy is mentioned as an option but not a confirmed requirement — flag for the planner to confirm with the product owner.

---
*Stack research for: Proteomics v1.4 — cross-team review, SPC drift tracking, cross-run screening comparison*
*Researched: 2026-06-06*
*Confidence: HIGH (every Datagrok API and Chem entrypoint verified against repo HEAD; SPC math is public domain; no-suitable-npm verdict verified against the two candidate libraries surfaced by search)*
