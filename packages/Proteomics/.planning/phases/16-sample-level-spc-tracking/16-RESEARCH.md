# Phase 16: Sample-Level SPC Tracking - Research

**Researched:** 2026-06-08
**Domain:** Per-run Shewhart I-MR + Nelson rules engine, AppData-backed multi-run log, in-dashboard baseline operator workflow — all platform-native (Datagrok JS-API + `_package.files` AppData), zero new npm dependencies.
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**D-01 — Run identity capture (SPC-06, Pitfall 7):** Inline in `showAnnotationDialog` (Annotate Experiment), optional. Two new inputs (`instrument_id` free-text string + `acquisition_datetime` date input), default empty. Spectronaut PG-report parser augmentation: extract `R.RunDate` (earliest observed) and `R.InstrumentMethod` (first observed); stash on the parsed df as `proteomics.spc_run_meta_seed` (JSON `{instrument_id, acquisition_datetime}`); annotation dialog reads this seed and prefills the two inputs. Other parsers leave the seed unset; the operator types both fields manually. Persisted on the df as `proteomics.spc_run_meta` (JSON `{instrument_id: string, acquisition_datetime: ISO-8601 string}`). Empty values tolerated at annotation time. `Analyze → Compute SPC Status` REFUSES on empty `spc_run_meta` with `"Open Annotate Experiment to set instrument + acquisition datetime."`. One df = one run = one `(instrument_id, acquisition_datetime)` pair. Spectronaut Candidates source is HARD-REFUSED for Compute SPC Status: `"SPC requires per-sample intensities. Re-import this analysis from the Spectronaut PG report (not the Candidates report) to compute SPC."`

**D-02 — Multi-run aggregation surface (SPC-02, SPC-05, SPC-07, Pitfalls 6/8):** AppData append-only SPC run log at `System:AppData/Proteomics/spc/runs.csv`. Compute SPC Status appends one row per (instrument_id, acquisition_datetime) primary key — second invocation OVERWRITES the prior row (idempotent). Columns: `run_id` (UUID), `instrument_id`, `acquisition_datetime` (ISO-8601), `run_label` (df.name), `median_intensity`, `missing_pct`, `control_corr`, `protein_count`, `status`, `rules_tripped` (JSON array), `source_project_id` (null when source df unpublished), `source_df_name`, `computed_at` (ISO-8601). Phase 17 additivity contract: Phase 17 reads the SAME shape and adds a `campaign_id` column when saving a run to a campaign (back-filling existing rows with NULL is fine). Phase 16 dashboard ignores `campaign_id` if present. No migration script. Schema version recorded in sibling `runs-meta.json` (`{spc.runs.schema_version: '1'}`). Locked baseline files (D-04) live next to runs.csv in the same `System:AppData/Proteomics/spc/` directory.

**D-03 — Control-replicate Pearson correlation definition (SPC-01):** Mean of pairwise Pearson correlation across Group 1 (Control) intensity columns. With Group 1 columns `[c_1, c_2, ..., c_k]`, compute `mean(Pearson(c_i, c_j))` for all unordered pairs `(i, j)` where `i < j`. NaN-tolerant pairwise: pairs with fewer than 3 jointly-observed proteins fall back to NaN and are excluded from the mean. If fewer than 2 control columns exist, `control_corr = NaN` and the rules engine treats NaN as "metric unavailable" (skip, not flag). Documented limitation: measures within-run reproducibility of biological controls; does NOT detect inter-run batch drift (Pitfall 12). Pooled-QC option B deferred to v1.5.

**D-04 — Baseline definition + Nelson rule toggling + persistence (SPC-03, SPC-05, SPC-08, Pitfall 6):** SPC Dashboard is the single surface. First open of `Visualize → SPC Dashboard...` with no baseline locked: empty-state banner `"No baseline locked for instrument <N>. [Define baseline...]"`. Click → modal listing all runs for that instrument from runs.csv with checkboxes (default: first 7 checked). Below: "Iterate 3σ outlier removal" toggle (default ON; capped at 2 iterations per Pitfall 6) and per-metric Nelson rule toggle grid (rules 1+5 default ON; 2-4 / 6-8 default OFF). On OK: compute mean+SD per metric over checked-and-non-outlier runs, write `System:AppData/Proteomics/spc/baseline-<instrument_id>.json` with schema (instrument_id, locked_at, included_run_ids, excluded_run_ids, iteration_trace, metrics{4 × {mean, sd}}, rules_enabled{4 × {nelson_1..nelson_8}}). One baseline.json per instrument_id. Rebuild OVERWRITES (no version history — v1.5+). Rule toggle = persisted (no separate save step). Pass/flagged/out_of_control classification: `pass` = no enabled rule trips; `out_of_control` = Nelson rule 1 (single point beyond 3σ) trips; `flagged` = any other enabled rule trips. Drill-down on flagged point: lookup runs.csv row's `source_project_id`; if non-null AND project still exists, opens via `grok.dapi.projects.find(id).open()`. Otherwise toast: `"Source for '<run_label>' is not currently available. Open the run's analyzed file and re-run Compute SPC Status."`. Pareto chart (SPC-08, P2): second tab/panel inside the dashboard; reads each row's `rules_tripped` JSON, counts per (rule, metric), bar-charts descending. Rules-tripped tag set unconditionally so deferred Pareto can land in v1.4.1.

### Claude's Discretion

- `proteomics.spc_metrics` tag shape: JSON `{median_intensity, missing_pct, control_corr, protein_count, sample_count, computed_at}` — matches D-02 columns minus run-identity fields (which live in `spc_run_meta`). Belt-and-braces follows Phase 15 D-05: write the same six values into a single-row `spc_metrics_meta` column so a published df's metrics survive Project serialization.
- `proteomics.spc_status` tag values: lowercase `'pass'` / `'flagged'` / `'out_of_control'`. New SEMTYPE deferred to Phase 17 (CAMP-07 plans `Proteomics-SpcStatus`); Phase 16 sets the tag and the dashboard reads it directly without semType-based detection.
- `proteomics.spc_rules_tripped` tag shape: JSON array of canonical rule names `['nelson_1@median_intensity', 'nelson_5@control_corr', ...]`. Empty array for `pass`.
- "Protein count above threshold" (SPC-01): `# proteins with non-null intensity in at least one sample column` — row count minus all-missing-rows. No user-tunable threshold for v1.4.
- Spectronaut PG-report run-identity coverage: `R.RunDate` parsed as ISO-8601, earliest observed; `R.InstrumentMethod` verbatim, first observed.
- Menu positions: `Analyze → Compute SPC Status` (immediate, no `...`); `Visualize → SPC Dashboard...` (`...` suffix). Last-picked instrument cached via Phase 13 D-02 precedent (see Pitfall 1 below — precedent code uses `grok.dapi.userDataStorage` despite a soft-deprecation; Phase 16 follows the same precedent).
- Dashboard implementation primitive: `DG.Viewer.lineChart` for I-chart + MR-chart with formula lines for UCL/CL/LCL. MR-chart against `abs(diff(metric))` columns computed at dashboard-open time. Pareto = `DG.Viewer.barChart`.
- runs.csv read concurrency: single-writer assumption. File-lock sentinel deferred unless concurrent corruption surfaces.

### Deferred Ideas (OUT OF SCOPE)

- Pooled-QC sample column as a first-class concept (revisit if Cytokinetics surfaces inter-run batch drift).
- Optional R `qcc` SPC fallback (REQUIREMENTS.md "Future Requirements", v1.5+).
- CUSUM / EWMA SPC charts.
- Per-protein curated panel SPC (Pitfall 5 explicit ban at this scope).
- `Proteomics-SpcStatus` SEMTYPE registration (Phase 17 CAMP-07 owns this).
- Cross-instrument comparison view ("compare instrument-1 to instrument-2 trend").
- Publish-block-on-QC-fail policy.
- Email / Slack / webhook alerting.
- Baseline version history / undo "Rebuild baseline".
- runs.csv concurrent-writer file lock.
- Per-protein drill-down from a flagged SPC point (drill-down opens source PROJECT only).
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| SPC-01 | Per-run sample-level QC metrics (median log2 intensity, percent missing, control-replicate Pearson correlation, protein count above threshold) computed when a run completes | `src/analysis/spc.ts` math (~80 LOC for the four metrics) on the source df; NaN-tolerant aggregation patterns already in `qc-computations.ts`; column lookup via `findColumn` + `getGroups(df)`; "protein count above threshold" definition lifted from Claude's discretion in CONTEXT.md |
| SPC-02 | SPC dashboard renders Shewhart I-chart and MR-chart per metric across a campaign's runs (via `DG.Viewer.lineChart` + formula lines for UCL/CL/LCL) | `DG.Viewer.lineChart` factory + `FormulaLinesHelper.addLine` API verified in js-api/src/{viewer.ts, helpers.ts}; MR derived from `abs(diff(metric))` over an instrument-filtered slice of runs.csv; existing `qc-dashboard.ts` reference-line pattern (M=0 line) is the direct precedent |
| SPC-03 | Nelson rules 1 (single point beyond 3σ) and 5 (2-of-3 beyond 2σ) enabled by default; rules 2-4 and 6-8 user-toggleable | `src/analysis/spc.ts` Nelson rules engine (~150 LOC TS, pure functions); per-metric toggle grid persisted into baseline-<instrument_id>.json's `rules_enabled` map; false-alarm-rate disclosure tooltip in dashboard sidebar (Pitfall 5 hygiene from research/PITFALLS.md) |
| SPC-04 | Each run tagged with `proteomics.spc_status` (`pass` / `flagged` / `out_of_control`) based on Nelson rule trips | Classification rule from D-04: rule-1 trip → `out_of_control`; any other enabled rule trip → `flagged`; none → `pass`. Set as a plain DataFrame tag on the analyzed df (no SEMTYPE; Phase 17 CAMP-07 owns the SEMTYPE addition) |
| SPC-05 | SPC baseline is an explicit annotated subset of runs, locked at definition with iterative outlier removal applied — never rolling-recomputed | D-04 modal flow + iterative 3σ removal capped at 2 iterations per Pitfall 6; baseline persisted as `baseline-<instrument_id>.json` next to runs.csv; rebuild OVERWRITES (no version history) |
| SPC-06 | Run identity is `(instrument_id, acquisition_datetime)` captured at annotate-experiment time; NOT file-import time | D-01 Annotate Experiment dialog extension (two new optional inputs); Spectronaut PG parser augmentation to populate `proteomics.spc_run_meta_seed` from `R.RunDate` + `R.InstrumentMethod`. Existing Spectronaut streaming parser does NOT capture these columns today — the column-index resolution at parsers/spectronaut-parser.ts:398-406 drops them; the augmentation MUST add them to that index map |
| SPC-07 | Clicking a flagged point on the SPC chart drills into the run's source DataFrame and table view | D-04 click handler reads runs.csv row's `source_project_id`; non-null → `grok.dapi.projects.find(id).open()` (pattern verified in `src/publishing/post-open-recovery.ts` and `src/publishing/assert-published-shape.ts` from Phase 15); null → toast (specific copy in UI-SPEC) |
| SPC-08 | SPC dashboard displays a Pareto chart showing which metric trips most often across the baseline window (P2) | `DG.Viewer.barChart` against an aggregated `(rule, metric, trip_count)` DataFrame computed at dashboard-open time from runs.csv `rules_tripped` JSON; rules-tripped tag set UNCONDITIONALLY so a deferred Pareto can land in a v1.4.1 hotfix |
</phase_requirements>

## Summary

Phase 16 ships a complete SPC story on platform-native primitives: `DG.Viewer.lineChart` + `FormulaLinesHelper` for I-MR charts (the math is locked in CONTEXT.md — Shewhart I-MR with d2=1.128, D4=3.267 for n=2 + Nelson rules 1-8), an append-only CSV in `System:AppData/Proteomics/spc/runs.csv` for the multi-run log (HitTriage's `_package.files.readAsText`/`writeAsText` precedent at `hit-design-app.ts:1056-1074` is the closest analog and confirms there is no Datagrok "append" primitive — the package reads-modifies-writes the whole CSV per Compute SPC Status invocation), a JSON baseline file per instrument, and a sidebar-driven dashboard pattern that mirrors `qc-dashboard.ts:25` and `enrichment-viewers.ts:191` (multi-viewer TableView built with `tv.dockManager.dock(...)`).

Three platform realities require planner attention (none are blockers, but they all need to be locked in plan tasks):

1. **`ui.input.dateTime` does NOT exist.** The UI-SPEC references it but the js-api ships only `ui.input.date(name, {value: dayjs.Dayjs})` (verified at `js-api/ui.ts:1080`). The plan must use `ui.input.date` and capture time-of-day either via a sibling `ui.input.string` for HH:MM or by accepting date-only precision (acquisition_datetime stored at midnight UTC of that calendar day). **Recommend: date-only precision.** The downstream SPC chart's X-axis is `acquisition_datetime` ordering — date-only resolution is sufficient for weekly assays and avoids two-input complexity. CONTEXT.md D-01 uses ISO-8601 strings, which accommodates both shapes.

2. **`grok.dapi.userDataStorage` is soft-deprecated** in favor of `grok.userSettings: UserSettingsStorage` (deprecation note at `js-api/src/dapi.ts:168`). Phase 13 D-02 precedent uses `grok.dapi.userDataStorage` (see `src/analysis/subcellular-location.ts:246` and `src/utils/gene-label-resolver.ts:139`). The dashboard's last-picked-instrument cache should follow the same precedent for consistency with the existing v1.3 codebase, with a note that v1.5 may migrate to `grok.userSettings` as part of a broader cleanup. Either API will work; the precedent wins.

3. **The streaming Spectronaut PG parser does NOT capture `R.RunDate` / `R.InstrumentMethod` today.** The column-index resolution at `parsers/spectronaut-parser.ts:398-406` only resolves `PG.ProteinGroups`, `R.Condition`, `R.Replicate`, the quantity column, `EG.Qvalue`, `R.FileName`, and `PG.Organisms`. The augmentation must (a) add these two columns to the index map IF the header includes them (optional — not in REQUIRED_COLUMNS so backward compatibility holds), (b) extract earliest-`R.RunDate` and first-`R.InstrumentMethod` during streaming aggregation, (c) store both in `proteomics.spc_run_meta_seed` on the parsed df. The same logic applies to the non-streaming text path. MaxQuant / FragPipe / Generic / Spectronaut Candidates parsers leave the seed unset.

**Primary recommendation:** Build `src/analysis/spc.ts` (math + tag helpers), `src/analysis/spc-storage.ts` (AppData I/O + sibling JSON baseline), `src/viewers/spc-dashboard.ts` (sidebar + four I-charts + four MR-charts + Pareto tab), plus parser augmentation in `src/parsers/spectronaut-parser.ts` (seed population) and dialog extension in `src/analysis/experiment-setup.ts` (two new optional inputs). Five new file additions, two existing files extended, one menu registration update. No new npm deps. Tests under `src/tests/spc.ts`.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Per-run SPC metric math (median, missing %, control-corr, protein count) | `src/analysis/spc.ts` (pure TS, deterministic) | — | Math layer — runs in browser, no R, no external API. Matches v1.3 `analysis/` pattern. |
| Nelson rules engine (rules 1-8) | `src/analysis/spc.ts` | — | Pure TS; ~150 LOC; public-domain math (Nelson 1984). |
| Run identity capture (instrument + acquisition_datetime) | `src/analysis/experiment-setup.ts` (existing dialog extension) | `src/parsers/spectronaut-parser.ts` (seed prefill from R.RunDate / R.InstrumentMethod) | D-01 inline extension. Other parsers leave seed unset. |
| AppData CSV runs.csv read/write + baseline JSON read/write | `src/analysis/spc-storage.ts` (new sibling to spc.ts) | — | File I/O kept separate from math so `computeSpcMetrics` stays pure. HitTriage `_package.files` precedent. |
| SPC Dashboard (4 I-charts + 4 MR-charts + Pareto + sidebar) | `src/viewers/spc-dashboard.ts` | `DG.Viewer.lineChart` + `formulaLines.addLine` (UCL/CL/LCL) + `DG.Viewer.barChart` (Pareto) | Mirrors `qc-dashboard.ts` shape: multi-viewer TableView with `tv.dockManager.dock(...)`. |
| Sidebar (instrument picker + baseline status + per-metric rule grid) | `src/viewers/spc-dashboard.ts` (composes `ui.divV` + `ui.input.choice` + `ui.input.bool` grid) + `tv.dockManager.dock(sidebar, DG.DOCK_TYPE.LEFT, null, ..., 0.25)` | — | Platform-native primitives; no shadcn, no third-party blocks. |
| Drill-down click-through (flagged point → source project) | `src/viewers/spc-dashboard.ts` click handler | `grok.dapi.projects.find(id).open()` (Phase 15 precedent) | Pure platform call. Toast fallback when project missing. |
| Last-picked-instrument cache across sessions | `src/viewers/spc-dashboard.ts` open path | `grok.dapi.userDataStorage` (Phase 13 D-02 precedent) | userDataStorage soft-deprecated but matches existing v1.3 codebase usage. |
| Test surface (SPC math + AppData round-trip + dashboard dock) | `src/tests/spc.ts` | `@datagrok-libraries/test` (existing test framework) | One file per Pitfall 5/6/7/8 lever + idempotency + drill-down. |

**Why no new directory:** SPC math is small (~390 LOC total) and fits cleanly under the existing `src/analysis/` / `src/viewers/` split (per ARCHITECTURE.md §4 — resolved divergence). Creating `src/spc/` would split the SPC story across two directories for no architectural gain.

## Standard Stack

### Core (platform-native — no installation)

| API | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `DG.Viewer.lineChart(df, options?)` | datagrok-api ^1.25.0 | I-chart + MR-chart per metric | Native viewer factory; supports the formula-lines API for UCL/CL/LCL [CITED: js-api/src/viewer.ts:243] |
| `Viewer.meta.formulaLines.addLine({formula, color, width, visible})` | datagrok-api ^1.25.0 | UCL = mean+3σ, CL = mean, LCL = mean-3σ horizontal lines on I-chart; UCL = D4×mean(MR), CL = mean(MR) on MR-chart | Verified in `js-api/src/helpers.ts:107` + existing `qc-dashboard.ts:88-94` precedent (M=0 reference line on MA scatter) [VERIFIED: js-api codebase] |
| `DG.Viewer.barChart(df, options?)` | datagrok-api ^1.25.0 | Pareto chart (P2) | Native viewer factory; standard for descending bar charts [CITED: js-api/src/viewer.ts] |
| `grok.dapi.files.exists(path)` / `.readAsText(path)` / `.readCsv(path)` / `.writeAsText(path, csvString)` | datagrok-api ^1.25.0 | runs.csv + baseline-<id>.json + runs-meta.json read/write under `System:AppData/Proteomics/spc/` | Verified in `js-api/src/dapi.ts:1298,1418,1426,1465`; HitTriage precedent at `hit-design-app.ts:1056-1074` and `hit-triage-app.ts:420` [VERIFIED: js-api codebase + HitTriage codebase] |
| `_package.files.readAsText('path')` / `_package.files.writeAsText('path', 'content')` | datagrok-api ^1.25.0 | Sibling baseline JSON read/write (relative paths within package's AppData) | HitTriage uses both `grok.dapi.files.*` (absolute `System:AppData/...` paths) and `_package.files.*` (relative within the package). Either works; planner picks the form consistent with v1.3 publishing code [VERIFIED: js-api codebase + HitTriage codebase] |
| `grok.dapi.projects.find(id)` + `project.open()` | datagrok-api ^1.25.0 | D-04 drill-down click → opens source project of flagged run | Verified in Phase 15 code at `src/publishing/post-open-recovery.ts` and `src/publishing/assert-published-shape.ts:64`; canonical pattern from `packages/Bio/src/tests/projects-tests.ts:26` [VERIFIED: in-repo Phase 15 code] |
| `grok.dapi.userDataStorage.get(name)` / `.put(name, data)` | datagrok-api ^1.25.0 | Last-picked-instrument cache (Phase 13 D-02 precedent) | Verified in `js-api/src/dapi.ts:716-751`; existing usage at `subcellular-location.ts:246` and `gene-label-resolver.ts:139`. Soft-deprecated in favor of `grok.userSettings` but matches v1.3 codebase precedent [VERIFIED: js-api codebase + Proteomics codebase] |
| `ui.input.string(name, {value, tooltipText})` | datagrok-api ^1.25.0 | `instrument_id` input on Annotate Experiment dialog | Standard pattern across the codebase [VERIFIED: package usage] |
| `ui.input.date(name, {value: dayjs.Dayjs, tooltipText})` | datagrok-api ^1.25.0 | `acquisition_datetime` input on Annotate Experiment dialog (date-only precision — see Summary §1 reality check) | Verified at `js-api/ui.ts:1080`; ApiSamples precedent at `packages/ApiSamples/scripts/ui/inputs/inputs.js:9` [VERIFIED: js-api codebase] |
| `ui.input.bool(name, {value, tooltipText})` | datagrok-api ^1.25.0 | "Iterate 3σ outlier removal" toggle + per-metric Nelson rule checkboxes in baseline modal | Existing usage in `enrichment.ts:500-505` [VERIFIED: package usage] |
| `ui.input.choice(name, {value, items, tooltipText})` | datagrok-api ^1.25.0 | Sidebar instrument picker | Existing usage at `share-dialog.ts:59` and elsewhere [VERIFIED: package usage] |
| `DG.TaskBarProgressIndicator.create(msg)` | datagrok-api ^1.25.0 | Optional progress indicator on Compute SPC Status (small compute — optional but consistent) | Standard pattern across `analysis/` files [VERIFIED: package usage] |
| `grok.shell.tv.dockManager.dock(element, DG.DOCK_TYPE.*, null, title, ratio)` | datagrok-api ^1.25.0 | Sidebar + chart panels in SPC Dashboard | Existing usage at `qc-dashboard.ts:142-149` and `enrichment-viewers.ts:222-240` [VERIFIED: package usage] |
| `grok.shell.addTableView(df)` | datagrok-api ^1.25.0 | Opens runs.csv-derived df in its own TableView (sibling-DF rule from ARCHITECTURE.md §7) | Existing usage at `enrichment.ts:574` and `pca-plot.ts` [VERIFIED: package usage] |

### Supporting (in-package, no external dependency)

| Library / Helper | Purpose | When to Use |
|---------|---------|-------------|
| `getGroups(df)` from `src/analysis/experiment-setup.ts` | Resolve Group 1 column names for control-corr computation | Every Compute SPC Status invocation [VERIFIED: existing helper] |
| `findColumn(df, semType, nameHints)` from `src/utils/column-detection.ts` | Locate intensity columns / gene columns when Group annotation is incomplete | Defensive lookups inside `computeSpcMetrics` [VERIFIED: existing helper] |
| `ensureFreshFloat(df, name)` pattern from `src/viewers/qc-computations.ts:27` | Idempotent re-write of `spc_metrics_meta` belt-and-braces column | Compute SPC Status second invocation must not duplicate columns [VERIFIED: existing pattern] |
| `dayjs` (already a webpack external) | Date/time math for acquisition_datetime + run sorting | Per `webpack.config.js` externals — already in `datagrok-api` peer set [VERIFIED: webpack config] |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| AppData CSV runs.csv | `grok.dapi.userDataStorage` (per-user JSON blob) | userDataStorage is per-user; runs.csv on AppData is per-machine. Cytokinetics-team shareability requires AppData. (CONTEXT.md D-02 locked the AppData choice.) |
| AppData CSV runs.csv | Custom Postgres schema under `databases/proteomics/` | Couples Phase 16 ship to a server config; over-engineered for ~104 rows/year. STACK.md "promote to DB at v1.5+ if O(1000)+" trigger applies. (Phase 17 CAMP scope, not Phase 16.) |
| `DG.Viewer.lineChart` formula lines | Hand-rolled overlay viewer | Loses tooltip integration and pan/zoom; reinventing UCL/CL/LCL rendering. (Phase 14 D-04 magenta/cyan/gray color lock does NOT apply — SPC is univariate.) |
| `grok.dapi.userDataStorage` (last-picked instrument) | `grok.userSettings` (newer API) | userSettings is the modern API; userDataStorage is soft-deprecated. v1.3 codebase uses userDataStorage exclusively; consistency wins. Track migration as future cleanup per memory `feedback_keep_workaround_capture_future.md`. |
| `ui.input.date` + sibling time input for acquisition_datetime | Use date-only precision (store at midnight UTC) | Two inputs add UX complexity; date-only is sufficient for weekly assays (CONTEXT.md describes weekly cadence). Recommend date-only. |

**Installation:**

No new packages. v1.3 `package.json` dependency set is unchanged.

**Version verification:** Not applicable — zero new npm dependencies. All capabilities are platform-native (datagrok-api ^1.25.0 already in v1.3 baseline).

## Package Legitimacy Audit

> **Not applicable.** Phase 16 installs zero external packages. The Package Legitimacy Gate protocol returns trivially: no slopcheck input, no registry verification, no postinstall script audit needed.

| Package | Registry | Disposition |
|---------|----------|-------------|
| (none) | — | — |

## Architecture Patterns

### System Architecture Diagram

```
                          PHASE 16 SPC DATA FLOW
                          ─────────────────────

   ┌──────────────────────────────────────────────────────┐
   │  Spectronaut PG parser (parsers/spectronaut-parser)  │
   │  EXTENSION: also extract earliest R.RunDate +        │
   │  first R.InstrumentMethod from header → set          │
   │  proteomics.spc_run_meta_seed (JSON) on parsed df    │
   └────────────────────┬─────────────────────────────────┘
                        │
                        ▼
   ┌──────────────────────────────────────────────────────┐
   │  Annotate Experiment dialog (analysis/experiment-    │
   │  setup.ts:showAnnotationDialog)                      │
   │  EXTENSION: 2 new optional inputs (instrument_id +   │
   │  acquisition_datetime), prefilled from seed.         │
   │  On OK → proteomics.spc_run_meta = JSON              │
   └────────────────────┬─────────────────────────────────┘
                        │
                        ▼
   ┌──────────────────────────────────────────────────────┐
   │  Analyze → Compute SPC Status (NEW menu)             │
   │  Gate checks:                                        │
   │    - df present                                      │
   │    - proteomics.source !== 'spectronaut-candidates'  │
   │    - proteomics.spc_run_meta non-empty               │
   │    - getGroups(df) !== null                          │
   │  Pipeline:                                           │
   │    1. computeSpcMetrics(df, groups, runMeta)         │
   │       → {median_intensity, missing_pct,              │
   │          control_corr, protein_count, sample_count}  │
   │    2. loadBaseline(instrument_id)                    │
   │    3. evaluateNelsonRules(metrics, baseline,         │
   │       rules_enabled) → {status, rulesTripped}        │
   │    4. setSpcStatus(df, metrics, ruleResult)          │
   │       → tags + spc_metrics_meta column               │
   │    5. appendRun(row)                                 │
   │       → AppData runs.csv (read whole CSV,            │
   │          upsert by (instrument_id, datetime),        │
   │          write whole CSV — idempotent)               │
   │    6. shell.info success toast                       │
   └────────────────────┬─────────────────────────────────┘
                        │
                        ▼
   ┌──────────────────────────────────────────────────────┐
   │  System:AppData/Proteomics/spc/                      │
   │    runs.csv                ← append-only (one row    │
   │                              per (instrument_id,     │
   │                              acquisition_datetime))  │
   │    baseline-<instrument>.  ← one per instrument;     │
   │      json                    overwritten on rebuild  │
   │    runs-meta.json          ← {schema_version: '1'}   │
   └────────────────────┬─────────────────────────────────┘
                        │
                        ▼
   ┌──────────────────────────────────────────────────────┐
   │  Visualize → SPC Dashboard... (NEW menu)             │
   │  Open path:                                          │
   │    - Check runs.csv exists with ≥1 row               │
   │    - Read last-picked-instrument from                │
   │      userDataStorage                                 │
   │    - Build instrument list (deduped from runs.csv)   │
   │    - Open instrument picker dialog                   │
   │    - On OK:                                          │
   │      → grok.shell.addTableView(loadRuns(instrument)) │
   │      → tv.dockManager.dock(sidebar, LEFT, 0.25)      │
   │      → tv.dockManager.dock(4 I-charts + 4 MR-charts) │
   │      → Pareto tab (P2)                               │
   │      → Wire I-chart point-click → drill-down handler │
   └────────────────────┬─────────────────────────────────┘
                        │
            ┌───────────┴────────────┐
            ▼                        ▼
   ┌────────────────┐      ┌─────────────────────┐
   │  Drill-down:   │      │  Rule toggle (live  │
   │  click flagged │      │  update from        │
   │  point         │      │  sidebar):          │
   │  → lookup      │      │  → update           │
   │    source_proj │      │    rules_enabled    │
   │  → projects.   │      │  → writeBaseline    │
   │    find(id).   │      │  → re-evaluate over │
   │    open()      │      │    all runs in slice│
   │  OR toast if   │      │  → update point     │
   │  null          │      │    colors in place  │
   └────────────────┘      └─────────────────────┘
```

### Recommended Project Structure

```
src/
├── analysis/
│   ├── (existing files unchanged)
│   ├── experiment-setup.ts          # EXTENDED: 2 new inputs in showAnnotationDialog
│   │                                #           + setRunMeta/getRunMeta helpers
│   ├── spc.ts                       # NEW: math layer
│   │                                #   - computeSpcMetrics(df, groups, runMeta)
│   │                                #   - evaluateNelsonRules(values, baseline, enabledRules)
│   │                                #   - setSpcStatus/getSpcStatus
│   │                                #   - I-chart constants (d2=1.128, D4=3.267)
│   └── spc-storage.ts               # NEW: AppData I/O layer
│                                    #   - appendRun(row)
│                                    #   - loadRuns(instrumentId?)
│                                    #   - loadBaseline(instrumentId)
│                                    #   - saveBaseline(instrumentId, baseline)
│                                    #   - readSchemaVersion / writeSchemaVersion
├── parsers/
│   └── spectronaut-parser.ts        # EXTENDED: streaming + text paths capture
│                                    #           R.RunDate (earliest) + R.InstrumentMethod
│                                    #           (first); set spc_run_meta_seed
├── viewers/
│   ├── (existing files unchanged)
│   └── spc-dashboard.ts             # NEW: openSpcDashboard,
│                                    #      createSpcChartPanel,
│                                    #      showDefineBaselineDialog,
│                                    #      createParetoChart (P2)
├── tests/
│   └── spc.ts                       # NEW: unit + integration tests
└── package.ts                       # EXTENDED: register 2 new menu items;
                                     #           gate Compute SPC Status by
                                     #           source / spc_run_meta / groups
```

### Pattern 1: Tag-helper getters/setters mirroring `setGroups` / `getGroups`

**What:** Wrap every `proteomics.spc_*` tag access through typed helpers in `src/analysis/spc.ts`. Never raw `df.setTag('proteomics.spc_*', ...)` in handlers.

**When to use:** All five new tags (`spc_status`, `spc_metrics`, `spc_rules_tripped`, `spc_run_meta`, `spc_run_meta_seed`).

**Example:**
```typescript
// Source: src/analysis/experiment-setup.ts:13-26 (existing precedent)
// Mirror in src/analysis/spc.ts:
export interface RunMeta {
  instrument_id: string;
  acquisition_datetime: string; // ISO-8601
}
export function setRunMeta(df: DG.DataFrame, meta: RunMeta): void {
  df.setTag('proteomics.spc_run_meta', JSON.stringify(meta));
}
export function getRunMeta(df: DG.DataFrame): RunMeta | null {
  const raw = df.getTag('proteomics.spc_run_meta');
  if (!raw) return null;
  try { return JSON.parse(raw) as RunMeta; } catch { return null; }
}
```

### Pattern 2: AppData CSV read-modify-write (no append primitive)

**What:** `grok.dapi.files.writeAsText(path, csvString)` overwrites the whole file every time. There is no Datagrok "append" primitive [VERIFIED: js-api/src/dapi.ts:1465 + HitTriage pattern `hit-design-app.ts:1056-1074` reads whole CSV via `grok.dapi.files.readCsv`, computes new content, writes whole CSV via `grok.dapi.files.writeAsText`]. The Phase 16 "append-only log" is therefore implemented as: read whole runs.csv → in-memory upsert by `(instrument_id, acquisition_datetime)` primary key → write whole runs.csv.

**When to use:** Every `appendRun(row)` invocation.

**Example:**
```typescript
// Pattern derived from packages/HitTriage/src/app/hit-design-app.ts:1056-1074
const path = 'System:AppData/Proteomics/spc/runs.csv';
let runsDf: DG.DataFrame;
if (await grok.dapi.files.exists(path)) {
  runsDf = await grok.dapi.files.readCsv(path);
} else {
  runsDf = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('run_id', []),
    DG.Column.fromStrings('instrument_id', []),
    // ... 13 columns total per CONTEXT.md D-02
  ]);
}
// upsert by composite key (instrument_id, acquisition_datetime)
upsertRow(runsDf, row);
await grok.dapi.files.writeAsText(path, runsDf.toCsv());
```

**Concurrency caveat:** CONTEXT.md D-02 explicitly defers file-locking ("single-writer assumption … `runs.csv.lock` sentinel only if research shows necessary"). The HitTriage pattern is also single-writer. No file lock for Phase 16.

### Pattern 3: Sidebar dock + multi-viewer dock with `tv.dockManager.dock`

**What:** Mirror `qc-dashboard.ts:142-149` and `enrichment-viewers.ts:222-240`: get `const dm = tv.dockManager`, dock the sidebar first (LEFT, 0.25), then dock each I-chart RIGHT/DOWN relative to the previous node.

**Example:**
```typescript
// Source: src/viewers/qc-dashboard.ts:142-149 (existing precedent)
const dm = tv.dockManager;
const sidebarNode = dm.dock(sidebarElement, DG.DOCK_TYPE.LEFT, null, 'SPC Controls', 0.25);
const iMedianNode = dm.dock(iChartMedian, DG.DOCK_TYPE.RIGHT, sidebarNode, 'Median Intensity (I)', 0.5);
dm.dock(mrChartMedian, DG.DOCK_TYPE.DOWN, iMedianNode, 'Median Intensity (MR)', 0.3);
// ... and so on for missing_pct, control_corr, protein_count
```

### Pattern 4: Reference lines on `DG.Viewer.lineChart` via `viewer.meta.formulaLines.addLine`

**What:** Render UCL = mean+3σ, CL = mean, LCL = mean-3σ as horizontal formula lines on each I-chart. For the MR-chart, only UCL = D4 × mean(MR) and CL = mean(MR) (LCL = 0 is the X-axis itself for n=2).

**Example:**
```typescript
// Pattern derived from src/viewers/qc-dashboard.ts:88-94 (existing precedent for M=0 reference line)
const iViewer = DG.Viewer.lineChart(slice, {
  xColumnName: 'acquisition_datetime',
  yColumnName: 'median_intensity',
});
iViewer.meta.formulaLines.addLine({
  formula: `\${median_intensity} = ${baseline.median_intensity.mean + 3 * baseline.median_intensity.sd}`,
  color: '#666666', width: 1, style: 'dashed', visible: true,
  title: 'UCL',
});
iViewer.meta.formulaLines.addLine({
  formula: `\${median_intensity} = ${baseline.median_intensity.mean}`,
  color: '#888888', width: 1, visible: true,
  title: 'CL',
});
iViewer.meta.formulaLines.addLine({
  formula: `\${median_intensity} = ${baseline.median_intensity.mean - 3 * baseline.median_intensity.sd}`,
  color: '#666666', width: 1, style: 'dashed', visible: true,
  title: 'LCL',
});
```

The `formulaLines.addLine` API is verified in `js-api/src/helpers.ts:107` and `FormulaLine` accepts arbitrary formula strings (interpreted by Datagrok's formula engine — `${col} = constant` renders a horizontal line at `y = constant` against any xColumnName).

### Pattern 5: Belt-and-braces tag + column encoding (Phase 15 D-05)

**What:** Critical metadata lives in BOTH a `proteomics.spc_*` tag AND a single-row `spc_metrics_meta` column on the analyzed df. Reason: serializer-strip behavior in `DG.Project.save` may drop tags but not column data (Phase 13 round-3 commit `e527d07ba1` evidence; Phase 15 Pitfall 3 mitigation).

**When to use:** Whenever SPC metric values must survive a Phase 15 publish-roundtrip.

**Example:**
```typescript
// Pattern from src/publishing/trim-dataframe.ts (Phase 15)
// In src/analysis/spc.ts:setSpcStatus(...)
df.setTag('proteomics.spc_metrics', JSON.stringify(metrics));
// AND
ensureFreshFloat(df, '~spc_median_intensity').init(() => metrics.median_intensity);
ensureFreshFloat(df, '~spc_missing_pct').init(() => metrics.missing_pct);
// ... per-field columns prefixed `~` so they don't render in default grids
// (HitTriage convention — `hit-design-app.ts:1072` filters out `~`-prefixed cols on save)
```

(Phase 15 publishing helper does NOT exclude SPC tags by default per CONTEXT.md "Claude's Discretion" — published projects carrying an SPC-annotated run can surface "this run was QC'd `pass` on 2026-07-10" via the reviewer panel.)

### Anti-Patterns to Avoid

- **Raw `df.setTag('proteomics.spc_*', ...)` in handlers:** breaks the v1.3 helper-discipline pattern (`getGroups` / `setGroups` / `markPublished` / `isPublished` in `publishing/publish-state.ts`). Always go through typed helpers.
- **Computing SPC math inside a viewer factory:** Pitfall 4-of-this-doc anti-pattern. Compute SPC Status mutates the source df + appends to runs.csv FIRST; dashboard reads pre-computed values from the AppData log.
- **Reusing v1.2 module-level `activeSubscriptions[]` for the I-chart click-handler:** Pitfall 13 from research/PITFALLS.md applies to N-producer/N-consumer wiring. Phase 16 is single-producer (one I-chart at a time), single-consumer (one drill-down handler). A plain `viewer.onCurrentRowChanged.subscribe(...)` + explicit unsubscribe-on-close is enough; **do NOT** build a CampaignSelectionBus here (Phase 18's territory).
- **Adding intensity-derived columns to runs.csv:** Pitfall 8 ban. runs.csv stores ONLY the four metric values per run × 13 cols. Per-protein drill-down opens the source project; NOT a row-expanded view.
- **Re-running SPC math from the viewer on every dashboard open:** Pitfall 8 performance trap. Read pre-computed metrics from runs.csv; only re-evaluate Nelson rules (cheap) when the user toggles a rule checkbox.
- **Order runs.csv by `computed_at` or insertion order:** Pitfall 7 ban. X-axis is always `acquisition_datetime` ASCENDING. Sort at dashboard-open time.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Horizontal reference lines on I-chart / MR-chart | Custom canvas overlay | `viewer.meta.formulaLines.addLine({formula: '${col} = constant', ...})` | Native platform feature; tooltip integration; pan/zoom-aware [VERIFIED: js-api/src/helpers.ts:107 + qc-dashboard.ts:88-94 precedent] |
| Pareto bar chart | Custom bar renderer | `DG.Viewer.barChart(df, {splitColumnName, valueColumnName, ...})` | Native viewer factory [VERIFIED: existing usage in qc-dashboard.ts:127] |
| AppData CSV write | `Blob` + `URL.createObjectURL` + hand-rolled download | `grok.dapi.files.writeAsText(path, df.toCsv())` | Native files API; supports `System:AppData/...` namespace; HitTriage precedent at `hit-design-app.ts:1074` [VERIFIED: js-api/src/dapi.ts:1465] |
| AppData CSV read | `fetch('/files/...')` + manual CSV parse | `grok.dapi.files.readCsv(path)` → `DG.DataFrame` | Native files API with typed columns; HitTriage precedent at `hit-triage-app.ts:90` [VERIFIED: js-api/src/dapi.ts:1426] |
| Cross-session "last picked instrument" cache | `localStorage` | `grok.dapi.userDataStorage.put/.get` | Per-user, cross-machine sync; existing Proteomics precedent in `subcellular-location.ts:246` [VERIFIED: existing pattern] |
| Pearson correlation between two columns | Custom implementation | `@datagrok-libraries/statistics`'s correlation primitives — BUT it ships `kendallsTau` only (verified at `libraries/statistics/src/correlation-coefficient.ts`), NOT Pearson. Build inline in `src/analysis/spc.ts` (~15 LOC NaN-tolerant pairwise; CONTEXT.md D-03 lock) | The library doesn't ship Pearson; inline build is small, public-domain, and deterministic [VERIFIED: grep across libraries/statistics] |
| Idempotent column re-creation | `df.columns.remove + addNewFloat` boilerplate at every call site | `ensureFreshFloat(df, name)` helper from `qc-computations.ts:27` | Existing precedent — exported helper, no need to duplicate [VERIFIED: existing helper] |
| Datetime parse / format for `acquisition_datetime` | `Date.parse` / `toISOString` boilerplate | `dayjs(value).toISOString()` — `dayjs` is already a webpack external for the package | Existing dependency; consistent across the codebase [VERIFIED: webpack config + ApiSamples usage] |
| Tag value serialization | Repeated `JSON.stringify` / `JSON.parse` with try/catch | Typed `setX(df, value)` / `getX(df): X \| null` helpers (Pattern 1 above) | v1.3 codebase convention from `experiment-setup.ts:13-26` [VERIFIED: existing pattern] |
| Dialog construction | Hand-rolled HTML | `ui.dialog(title).add(input).onOK(...).show()` | Existing precedent across every dialog in `analysis/` and `publishing/` [VERIFIED: existing pattern] |

**Key insight:** Phase 16 has *no math research to do* (Nelson 1984 + Wheeler constants locked in CONTEXT.md), and *no new external code* to write. Every primitive needed is already in v1.3's dependency set. The work is composition, not invention.

## Common Pitfalls

> Phase 16 inherits four named pitfalls from `.planning/research/PITFALLS.md` (5/6/7/8) and one platform-API pitfall surfaced during this research.

### Pitfall 1: userDataStorage API is soft-deprecated but v1.3 codebase uses it exclusively

**What goes wrong:** Planner reads the modern `grok.userSettings` docs, picks that API for the last-picked-instrument cache, and ends up with inconsistent storage layers across the v1.3 codebase (subcellular-location.ts and gene-label-resolver.ts use userDataStorage; the new SPC code would use userSettings).

**Why it happens:** `js-api/src/dapi.ts:168` and `:712` both carry `@deprecated The UserDataStorage should not be used. Use {@link UserSettingsStorage} instead`. The deprecation note is real but the codebase has not migrated. Phase 13 D-02 precedent is userDataStorage.

**How to avoid:** Use `grok.dapi.userDataStorage` to match the v1.3 precedent (Phase 13 D-02). Track migration to `grok.userSettings` as a cross-cutting cleanup todo for v1.5+ per memory `feedback_keep_workaround_capture_future.md` — keep the workaround running, capture cleanup as a future-action todo with the release-version trigger.

**Warning signs:** Mixed usage of `grok.dapi.userDataStorage` and `grok.userSettings` in the same package; inconsistent storage keys.

### Pitfall 2 (inherited Pitfall 5 — Per-protein SPC alarm flood)

**What goes wrong:** Implementer expands SPC scope from 4 run-level metrics to per-protein SPC. At 5,000 proteins × 8 Nelson rules, false-alarm rate goes to ~133 alarms per run. Reviewers dismiss SPC.

**Why it happens:** The package's data model is protein-centric. Easy to confuse "sample-level intensity median" with "per-protein median across samples."

**How to avoid:** LOCK scope to 4 run-level metrics (CONTEXT.md D-domain pin). Default Nelson rules to 1+5 only (combined false-alarm rate ~0.4%). Dashboard sidebar shows a live false-alarm-rate disclosure tooltip as the user toggles rules ON. **Test/UAT lever:** Unit test simulates 100 runs of in-control data and asserts < 5 total alarms under default rules.

**Warning signs:** A "by protein" picker on the SPC dashboard; SPC chart shows hundreds of flagged runs.

### Pitfall 3 (inherited Pitfall 6 — Baseline contamination)

**What goes wrong:** Baseline computed from "all runs to date" includes a real special-cause outlier; SD inflates; control limits widen; future drifts no longer trip.

**How to avoid:** D-04 modal flow makes baseline an EXPLICIT operator-curated subset (checkboxes). "Iterate 3σ outlier removal" toggle is ON by default and capped at 2 iterations. Excluded runs recorded in `excluded_run_ids` + `iteration_trace` in baseline JSON. Rebuild requires explicit user action. **Test/UAT lever:** Unit test: baseline of `[normal, normal, normal, normal, normal, 3-sigma outlier]` with iterate toggle ON returns limits computed from the first 5, plus `excluded_run_ids = [outlier_id]`.

**Warning signs:** Rolling-window baseline option in the dashboard; no "exclude this run" affordance.

### Pitfall 4 (inherited Pitfall 7 — Run-order ambiguity)

**What goes wrong:** SPC charts ordered by `computed_at` or file-import time. Multi-instrument runs interleave. Nelson rule "9 points on one side" fires from interleaving artifact, not biology.

**How to avoid:** D-01 captures `(instrument_id, acquisition_datetime)` at annotation time. Dashboard X-axis = `acquisition_datetime` ascending, filtered by `instrument_id`. Backfilled runs sort correctly by acquisition_datetime. **Test/UAT lever:** Unit test inserts run-X (acquisition_datetime = week 0) then run-Y (acquisition_datetime = week -1, backfill); asserts dashboard chart orders Y before X.

**Warning signs:** SPC chart x-axis labeled "Run 1, 2, 3..." instead of acquisition date.

### Pitfall 5 (inherited Pitfall 8 — Storage cost growth)

**What goes wrong:** runs.csv stores full per-protein quant matrices. After 52 weeks runs.csv hits hundreds of MBs and the dashboard takes >5s to open.

**How to avoid:** runs.csv stores ONLY the per-run rollup: 13 columns × ~104 rows/year. Per-protein drill-down NEVER lives in runs.csv; drill-down opens `source_project_id` via `grok.dapi.projects.find(id).open()`. **Test/UAT lever:** Unit test: after N synthetic SPC computations on N distinct (instrument, datetime) keys, runs.csv has exactly N rows × 13 columns.

**Warning signs:** runs.csv column count > 13; "should we store the quant matrix too" review comment.

## Runtime State Inventory

> Phase 16 is a **greenfield** addition (new files + small extensions to two existing files), NOT a rename/refactor/migration. The Runtime State Inventory categories from the GSD framework do not apply.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — verified: no `proteomics.spc_*` tags exist in repo today (grep confirmed). New AppData directory `System:AppData/Proteomics/spc/` is created on first Compute SPC Status invocation; no migration. | None |
| Live service config | None — no n8n / Datadog / Cloudflare integrations | None |
| OS-registered state | None — no scheduled tasks, no daemons | None |
| Secrets/env vars | None — no API keys, no credentials | None |
| Build artifacts | None — no compiled binaries, no egg-info, no Docker tags | None |

**Nothing found in any category:** This is a greenfield phase. Verified by `grep -rn 'spc_\|spc_status\|spc_metrics\|spc_run' src/` — returns zero hits.

## Environment Availability

> Phase 16 has **zero new external dependencies**. It runs entirely against `datagrok-api ^1.25.0` (already in v1.3 baseline). No Python, no R, no Docker, no databases, no CLIs.

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| `datagrok-api` | All Phase 16 code | ✓ | ^1.25.0 (already in v1.3) | — |
| `dayjs` | acquisition_datetime parse/format | ✓ | webpack external (already in v1.3) | — |
| Datagrok server / shell | All UI surfaces + AppData I/O | ✓ | Running instance required for `grok test` and live development | — |

**No missing dependencies. No fallback needed.**

## Validation Architecture

> Nyquist validation is enabled (config.json `workflow.nyquist_validation` absent → treat as enabled). This section is consumed by VALIDATION.md.

### Test Framework

| Property | Value |
|----------|-------|
| Framework | `@datagrok-libraries/test` (existing — used by every `src/tests/*.ts` file) |
| Config file | None — tests live in `src/tests/<file>.ts` and are imported by `src/package-test.ts` |
| Quick run command | `grok test --host localhost --category "SPC"` |
| Full suite command | `grok test --host localhost` |
| Build precondition | `npm run build` MUST run before tests; per memory `feedback_grok_test_skipbuild_stale`, `--skip-build` reuses stale bundles when new tests are added |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SPC-01 (median_intensity) | `computeSpcMetrics(df, groups, runMeta)` returns the correct median of all log2 intensity values across all sample columns, NaN-skipping nulls | unit | `grok test --host localhost --test "SPC:median_intensity"` | ❌ Wave 0 — new `src/tests/spc.ts` |
| SPC-01 (missing_pct) | Returns the correct percentage of (protein × sample) cells that are null across all intensity columns | unit | `grok test --host localhost --test "SPC:missing_pct"` | ❌ Wave 0 |
| SPC-01 (control_corr) | Returns NaN when fewer than 2 control columns exist; returns mean-of-pairwise-Pearson over Group 1 columns; pairs with < 3 jointly-observed proteins excluded from mean | unit | `grok test --host localhost --test "SPC:control_corr"` | ❌ Wave 0 |
| SPC-01 (protein_count) | Returns rowCount minus all-missing-rows | unit | `grok test --host localhost --test "SPC:protein_count"` | ❌ Wave 0 |
| SPC-02 (chart wiring) | `openSpcDashboard()` against a synthetic 8-row runs.csv produces 4 I-chart + 4 MR-chart viewers in the table view's dockManager | integration | `grok test --host localhost --test "SPC:dashboard_renders"` | ❌ Wave 0 |
| SPC-02 (formula lines) | Each I-chart has 3 formula lines (UCL/CL/LCL) with the correct y-values from `baseline.json`; MR-chart has 2 formula lines (UCL/CL) | integration | `grok test --host localhost --test "SPC:formula_lines"` | ❌ Wave 0 |
| SPC-03 (default rules) | `evaluateNelsonRules` with rules 1+5 enabled on a known synthetic series tripping rule 5 returns `status: 'flagged'`, `rulesTripped: ['nelson_5@median_intensity']` | unit | `grok test --host localhost --test "SPC:nelson_default"` | ❌ Wave 0 |
| SPC-03 (false-alarm rate) | 100 in-control synthetic runs (drawn from a single Gaussian) with rules 1+5 enabled produce < 5 total alarms (Pitfall 5 lever) | unit | `grok test --host localhost --test "SPC:false_alarm_rate"` | ❌ Wave 0 |
| SPC-04 (status tags) | `setSpcStatus(df, metrics, ruleResult)` writes `proteomics.spc_status`, `proteomics.spc_metrics` (JSON), `proteomics.spc_rules_tripped` (JSON array), AND the belt-and-braces `spc_metrics_meta` column on the df | unit | `grok test --host localhost --test "SPC:status_tags"` | ❌ Wave 0 |
| SPC-04 (classification) | rule-1 trip → `'out_of_control'`; rule-5 only → `'flagged'`; no rules trip → `'pass'` | unit | `grok test --host localhost --test "SPC:classification"` | ❌ Wave 0 |
| SPC-05 (iterative outlier removal) | Baseline of `[normal, normal, normal, normal, normal, 3-sigma_outlier]` with iterate toggle ON returns mean/SD from the first 5 only; `excluded_run_ids = [outlier_id]`; `iteration_trace` records the iteration that dropped the outlier (Pitfall 6 lever) | unit | `grok test --host localhost --test "SPC:baseline_outlier_removal"` | ❌ Wave 0 |
| SPC-05 (iteration cap) | Pathological "every iteration drops one point" series stops after exactly 2 iterations (cap from CONTEXT.md D-04) | unit | `grok test --host localhost --test "SPC:baseline_iteration_cap"` | ❌ Wave 0 |
| SPC-05 (locked baseline persistence) | After `saveBaseline(instrument, baseline)`, `loadBaseline(instrument)` returns the exact same metric mean/SD values and rule toggles | integration | `grok test --host localhost --test "SPC:baseline_roundtrip"` | ❌ Wave 0 |
| SPC-05 (rebuild overwrites) | Calling `saveBaseline` twice for the same instrument leaves exactly one `baseline-<instrument>.json` file with the second invocation's values | integration | `grok test --host localhost --test "SPC:baseline_rebuild_overwrites"` | ❌ Wave 0 |
| SPC-06 (run identity capture) | `showAnnotationDialog` extended with two new optional inputs; `setRunMeta(df, meta)` persists `proteomics.spc_run_meta` as JSON `{instrument_id, acquisition_datetime}`; `getRunMeta(df)` returns parsed object or null | unit | `grok test --host localhost --test "SPC:run_meta_helpers"` | ❌ Wave 0 |
| SPC-06 (Spectronaut seed) | Spectronaut PG-report header containing `R.RunDate` and `R.InstrumentMethod` columns → parsed df carries `proteomics.spc_run_meta_seed` with earliest-`R.RunDate` and first-`R.InstrumentMethod` (Pitfall 7 lever) | unit | `grok test --host localhost --test "SPC:spectronaut_seed"` | ❌ Wave 0 |
| SPC-06 (backfill ordering) | Backfill-inserted run (acquisition_datetime = week -1) appears BEFORE the previously-inserted week-0 run when the SPC dashboard sorts the instrument's slice (Pitfall 7 lever) | integration | `grok test --host localhost --test "SPC:backfill_ordering"` | ❌ Wave 0 |
| SPC-06 (Candidates refuse) | Compute SPC Status on a df with `proteomics.source === 'spectronaut-candidates'` surfaces the hard-refuse warning and does NOT touch runs.csv | integration | `grok test --host localhost --test "SPC:candidates_refuse"` | ❌ Wave 0 |
| SPC-07 (drill-down resolved) | Click handler on flagged I-chart point with non-null `source_project_id` (and project still exists) calls `grok.dapi.projects.find(id).open()` | integration | `grok test --host localhost --test "SPC:drilldown_resolved"` | ❌ Wave 0 |
| SPC-07 (drill-down null/missing) | Click on flagged point with null `source_project_id` surfaces the correct toast (specific copy in UI-SPEC); click on flagged point with non-null id but `projects.find` returns null surfaces a different toast | integration | `grok test --host localhost --test "SPC:drilldown_missing"` | ❌ Wave 0 |
| SPC-08 (Pareto chart) | With a runs.csv containing 10 runs whose `rules_tripped` JSON arrays sum to {`nelson_1@median_intensity`: 4, `nelson_5@control_corr`: 3, ...}, the Pareto bar chart renders bars in descending count order | integration | `grok test --host localhost --test "SPC:pareto_descending"` | ❌ Wave 0 |
| Idempotency (D-02 key uniqueness) | Two consecutive `Compute SPC Status` invocations on the same df (same instrument_id + acquisition_datetime) leave runs.csv with exactly one row, status reflecting the second invocation | integration | `grok test --host localhost --test "SPC:idempotent_upsert"` | ❌ Wave 0 |
| Storage growth bound (Pitfall 8 lever) | After N synthetic SPC computations on N distinct (instrument, datetime) keys, runs.csv has exactly N rows × 13 columns | unit | `grok test --host localhost --test "SPC:storage_bounded"` | ❌ Wave 0 |
| Schema versioning (D-02 sibling file) | First-ever Compute SPC Status creates `runs-meta.json` with `spc.runs.schema_version: '1'` | integration | `grok test --host localhost --test "SPC:schema_version"` | ❌ Wave 0 |
| Belt-and-braces column (Phase 15 carry) | Setting SPC status writes BOTH the `proteomics.spc_metrics` tag AND `~spc_*` per-field columns | unit | `grok test --host localhost --test "SPC:belt_and_braces"` | ❌ Wave 0 |
| Ensure idempotent column re-write | Re-running Compute SPC Status on the same df does NOT produce duplicate `~spc_*` columns (`ensureFreshFloat` pattern) | unit | `grok test --host localhost --test "SPC:column_idempotent"` | ❌ Wave 0 |

### Sampling Rate

- **Per task commit:** `grok test --host localhost --category "SPC"` — runs every SPC test (~25 tests; expected runtime < 30s for unit + integration since math is pure TS).
- **Per wave merge:** `grok test --host localhost` — full Proteomics suite (existing v1.3 tests + SPC + publish-roundtrip from Phase 15) to verify no regression in upstream tests when shared files (`experiment-setup.ts`, `spectronaut-parser.ts`, `package.ts`) are edited.
- **Phase gate:** `grok test --host localhost --verbose` (full suite green) before `/gsd:verify-work` runs the live `/gsd:uat` flow.

### Wave 0 Gaps

- [ ] `src/tests/spc.ts` — covers SPC-01 through SPC-08 + idempotency + belt-and-braces (~25 tests).
- [ ] No new shared `conftest.py`-equivalent — `src/tests/` files are self-contained per v1.3 convention.
- [ ] Framework install: none — `@datagrok-libraries/test` already a dev dependency.
- [ ] Add `import './tests/spc';` to `src/package-test.ts` so the new tests are discovered.

**Synthetic fixtures Wave 0 must build (inside `src/tests/spc.ts`):**
- `makeInControlSeries(n, mean, sd) → number[]` — Gaussian draws for false-alarm-rate test.
- `makeOutOfControlSeries(n, mean, sd, outlierIndices, outlierMagnitude) → number[]` — controlled outlier injection for Nelson rule tests.
- `makeSyntheticRunsDf(n, instrumentId) → DG.DataFrame` — 13-column runs.csv shape for AppData round-trip tests.
- `makeSyntheticSpectronautHeader({withRunDate, withInstrumentMethod, ...}) → string` — header lines for seed-extraction tests.

## Security Domain

> CONTEXT.md does not declare `security_enforcement` explicitly. Following the GSD default (absent = enabled), this section is included.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | Datagrok platform owns authentication; package operates inside an authenticated session |
| V3 Session Management | no | Platform-owned |
| V4 Access Control | partial | `grok.dapi.projects.find(id).open()` respects the calling user's project ACL — if the user cannot read the source project, the drill-down silently falls back to the toast. No new permission checks required in Phase 16 code. |
| V5 Input Validation | yes | `instrument_id` is free-text — sanitize before using as a filename component (`baseline-<instrument_id>.json`). Restrict to `[A-Za-z0-9._-]+` and reject path traversal characters (`/`, `\`, `..`). Pattern: borrow the slugify rule from Phase 15 D-01 (target slugification at `src/publishing/publish-state.ts`). |
| V6 Cryptography | no | No crypto; no secrets in Phase 16 |

### Known Threat Patterns for Datagrok JS-API + AppData CSV

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Path traversal via `instrument_id` (`../../../etc/passwd`-style baseline filename) | Tampering | Slugify `instrument_id` before passing to `baseline-<slug>.json`. Reuse the Phase 15 slug helper from `src/publishing/publish-state.ts`. |
| CSV injection in `run_label` / `instrument_id` (e.g. `=cmd|...`) when a downstream tool opens runs.csv in Excel | Tampering | runs.csv is consumed by Datagrok's CSV parser, not Excel. No external Excel export path in Phase 16. Document the constraint; defer Excel-safety to v1.5 if a download surface ships. |
| Race condition on runs.csv concurrent writes | Repudiation / Tampering | CONTEXT.md D-02 deferred file-locking under single-writer assumption. If/when concurrent corruption surfaces, ship a `runs.csv.lock` sentinel + 100ms retry. |
| AppData visibility — runs.csv readable by any Datagrok user with AppData read | Information Disclosure | All Cytokinetics proteomics team members have legitimate read access to the SPC log (SPC is operator-collaboration data). No PII / PHI / compound identity in runs.csv. Per-project drill-down still respects each project's ACL. |
| Tag forgery in `proteomics.spc_status` from a malicious user mutating a df before publishing | Tampering | Phase 16 sets the tag; downstream consumers (Phase 17 campaign index, future publish-block-on-QC-fail policy) MUST re-verify by re-reading runs.csv rather than trusting the df tag. Phase 16 documents this contract for Phase 17 to honor. |

**Phase 16 ships no new endpoints, no new external API integrations, no new credentials.**

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Hand-rolled SPC libraries on npm (`process-control-charts`, `spc-limits-js`) | In-package ~390 LOC TS (Shewhart I-MR + Nelson rules 1-8) | 2026-06 (Phase 16 STACK.md verdict) | Avoids the abandoned-looking npm landscape; math is public-domain (Nelson 1984, Wheeler 1992); deterministic enough for client-side TS path |
| `grok.dapi.userDataStorage` | `grok.userSettings: UserSettingsStorage` | (deprecation note in datagrok-api ^1.25.0 at `js-api/src/dapi.ts:168`) | Phase 16 follows v1.3 precedent (userDataStorage) for consistency; v1.5+ cleanup to migrate the whole codebase at once per memory `feedback_keep_workaround_capture_future` |
| `ui.input.dateTime` | `ui.input.date(name, {value: dayjs.Dayjs})` | (UI-SPEC drafted with mistaken assumption; verified at js-api research time) | Phase 16 uses date-only precision for `acquisition_datetime`; weekly cadence makes time-of-day irrelevant. Plan must drop the dateTime input reference in UI-SPEC. |

**Deprecated/outdated to avoid:**
- `process-control-charts` (npm) — attribute charts only (p/np/c/u); wrong type for continuous metrics.
- `spc-limits-js` (npm/GitHub) — repo not resolvable.
- R `qcc` package on the critical path — REQUIREMENTS.md "Future Requirements" parks this for v1.5+.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `viewer.meta.formulaLines.addLine({formula: '${col} = constant', ...})` renders a horizontal line at the constant Y value, independent of which `xColumnName` the lineChart uses | Pattern 4 | If formula engine rejects the constant-Y form, fall back to a sibling DataFrame column holding the constant value at every row and use `${ucl_col}` style formula. Verified at QC dashboard precedent for `${M} = 0` on a scatter plot (`qc-dashboard.ts:88-94`) — line chart formula engine should follow the same shape but live render not done in research. |
| A2 | `grok.dapi.files.writeAsText` is atomic at the filesystem level (no partial-write corruption on concurrent invocation) | Pattern 2 / Pitfall on concurrency | If non-atomic, single-writer assumption breaks; runs.csv could be corrupted. HitTriage's identical pattern at `hit-design-app.ts:1074` has shipped to production without lock, suggesting acceptable single-writer semantics. CONTEXT.md D-02 deferred file-locking explicitly. |
| A3 | Spectronaut PG-report headers reliably include `R.RunDate` and `R.InstrumentMethod` on the standard export | D-01 / SPC-06 | If columns are absent, seed remains empty; operator types both manually. The augmentation is purely additive (no required-columns check change). |
| A4 | `R.RunDate` values are parseable as ISO-8601 or via `dayjs(value)` | D-01 | If they are formatted as Excel-style serials or vendor-specific strings, the parser falls back to setting the seed to the raw string and lets the operator correct via the Annotate Experiment dialog. |
| A5 | The `~spc_*` belt-and-braces column prefix `~` hides those columns from default grid views (HitTriage convention) | Pattern 5 | If the platform's default grid renders `~`-prefixed columns, the published reviewer panel may surface clutter. Verify against the existing v1.3 Filters viewer column filter at plan-check time. |
## Open Questions (RESOLVED)

1. **Q1: Does `DG.Viewer.lineChart`'s formula-lines API support per-line `style: 'dashed'` rendering?**
   - What we know: `FormulaLine.style?: string` is a documented field at `js-api/src/helpers.ts:27`. The `qc-dashboard.ts:88-94` precedent only sets solid M=0 reference lines, so dashed-style rendering on lineChart is unverified.
   - What's unclear: whether `style: 'dashed'` is honored on `DG.Viewer.lineChart` (vs. scatter plot where M=0 ships solid).
   - RESOLVED: Deferred to the Wave-0 spike in Plan 16-01. The spike test runs to completion and writes `16-01-SPIKE-RESULT.md` with `OUTCOME: PASS` (dashed style honored) or `OUTCOME: FAIL` (platform error captured verbatim). Plan 16-06's `createSpcChartPanel` branches on this marker: PASS → use `style: 'dashed'` for UCL/LCL; FAIL → fallback path = constant-column overlay (a second yColumn on the same lineChart holding the constant UCL value at every row). The fallback is documented in 16-06-PLAN.md Task 1 `<behavior>` step 2.

2. **Q2: How does `formulaLines.addLine` interact with `DG.Viewer.lineChart` when `xColumnName` is a datetime column?**
   - What we know: scatter plot precedent at `qc-dashboard.ts:88-94` works against a numeric X column ('A').
   - What's unclear: whether `${median_intensity} = 22.31` renders correctly as a horizontal line when the X-axis is `acquisition_datetime`.
   - RESOLVED: Same Wave-0 spike covers both questions. The spike result file records the platform's behavior verbatim. If broken, alternative is the constant-column overlay path documented in Q1's RESOLVED note.

3. **Q3: Can the SPC Dashboard sidebar persist its rule-toggle state across dashboard re-opens for a SECOND instrument?**
   - What we know: D-04 locks "Rule toggle = persisted (no separate save step)" — writes go to `baseline-<instrument_id>.json`. The dashboard switches instruments by switching baseline files.
   - What's unclear: Whether toggling rules for instrument A and then switching to instrument B (different baseline file) leaves A's toggles intact.
   - RESOLVED: One `baseline-<instrument_id>.json` per instrument (CONTEXT.md D-04). Switching instrument in the sidebar reloads the corresponding baseline.json including its `rules_enabled` map, so per-instrument isolation is structural. Verified by unit test `SPC:rule_toggle_per_instrument` in Plan 16-01.

4. **Q4: Should `Compute SPC Status` show a TaskBarProgressIndicator for the small compute?**
   - What we know: CONTEXT.md "Claude's Discretion" says "good practice but optional". The compute is small (~4 metrics × O(proteins × samples)).
   - What's unclear: Whether the compute crosses the 500ms threshold from Pitfall 14 (every >500ms surface gets a progress indicator) on a large 10k-protein × 50-sample dataset.
   - RESOLVED: Add the indicator unconditionally for Compute SPC Status (Plan 16-05 Task 2). The compute is small but consistent UX > conditional logic. Matches v1.3 codebase convention (`differential-expression.ts:295`, `enrichment.ts:565`, every other long-ish op).

## Sources


### Primary (HIGH confidence)

- **js-api codebase (read at HEAD of worktree):**
  - `js-api/src/dapi.ts:170,716-751,1287-1465` — UserDataStorage + FilesDataSource (readAsText, writeAsText, readCsv, exists)
  - `js-api/src/viewer.ts:243,792-814` — `DG.Viewer.lineChart`, ViewerFormulaLinesHelper
  - `js-api/src/helpers.ts:9-32,83-141` — FormulaLine interface + FormulaLinesHelper.addLine
  - `js-api/src/dataframe/formula-helpers.ts` — DataFrameFormulaLinesHelper (per-df formula lines)
  - `js-api/src/user_settings_storage.ts` — modern API (soft-deprecates userDataStorage)
  - `js-api/ui.ts:831,1044-1098` — `ui.input.string` / `int` / `bool` / `date` factories
- **Proteomics package codebase (in-repo precedents):**
  - `src/analysis/experiment-setup.ts:13-26,29-60,63-97` — setGroups/getGroups tag-helper pattern + seedAnnotationDialogInputs + showAnnotationDialog (D-01 extension point)
  - `src/parsers/spectronaut-parser.ts:12-17,389-406` — Spectronaut PG parser column-resolution map (D-01 augmentation site)
  - `src/viewers/qc-computations.ts:13-32` — groupMean + ensureFreshFloat pattern
  - `src/viewers/qc-dashboard.ts:25,88-94,142-149` — multi-viewer dock pattern + formula-line precedent
  - `src/viewers/enrichment-viewers.ts:170-240` — sibling DataFrame TableView + tv.dockManager.dock pattern
  - `src/analysis/enrichment.ts:574` — grok.shell.addTableView for sibling df
  - `src/publishing/assert-published-shape.ts:64` + `src/publishing/post-open-recovery.ts` — `grok.dapi.projects.find(id).open()` precedent for D-04 drill-down
  - `src/publishing/trim-dataframe.ts` — Phase 15 belt-and-braces column encoding pattern
  - `src/utils/column-detection.ts:6-28` — findColumn / findProteomicsColumns helpers
- **HitTriage package codebase (cross-package precedents):**
  - `packages/HitTriage/src/app/hit-design-app.ts:1056-1074` — read-modify-write CSV pattern under AppData; verified single-writer assumption shipped to production
  - `packages/HitTriage/src/app/hit-triage-app.ts:90,375,420` — `_package.files.readCsv` + `_package.files.writeAsText` precedent + AppData directory layout

### Secondary (MEDIUM confidence — research-doc-citation only, not live-tested in this research session)

- `.planning/research/PITFALLS.md` §"Pitfall 5/6/7/8" — Nelson rule false-alarm rates, baseline iterative removal recommendations, run-order ambiguity
- `.planning/research/STACK.md` §"In-package SPC" — ~390 LOC TS estimate; verified no Shewhart/Nelson in `@datagrok-libraries/statistics`
- `.planning/research/ARCHITECTURE.md` §3/§4 — tag contract extensions + directory structure
- `.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md` — belt-and-braces D-05 pattern + userDataStorage cross-session precedent

### Tertiary (LOW confidence — public-domain references already verified upstream, not re-verified in this research)

- Nelson 1984 — Nelson rules 1-8 (cited in CONTEXT.md canonical_refs).
- Wheeler & Chambers 1992 — d2=1.128, D4=3.267 for n=2 MR-chart constants (cited in CONTEXT.md canonical_refs).

## Project Constraints (from CLAUDE.md)

Phase 16 must honor all of these or fail review:

- **2-space indent, single quotes, semicolons, CRLF line endings, TypeScript strict mode** (root CLAUDE.md + package CLAUDE.md).
- **File naming: kebab-case** — `spc-dashboard.ts`, `spc-storage.ts`, NOT `spcDashboard.ts`.
- **Never edit `.g.ts` / `.api.g.ts`** — auto-generated by `grok api`.
- **Three import lines required at the top of every TS source:**
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as ui from 'datagrok-api/ui';
  import * as DG from 'datagrok-api/dg';
  ```
- **No barrel files (no `index.ts` re-exports).** Imports go directly to the source file.
- **No default exports** — all exports are named.
- **Webpack externals not bundled:** `datagrok-api/dg`, `datagrok-api/grok`, `datagrok-api/ui`, `rxjs`, `rxjs/operators`, `cash-dom`, `dayjs`, `wu`, `openchemlib/full.js`.
- **Server API discipline:** never raw `fetch('/api/...')`; never raw `fetch('https://external.com/...')` for external CORS. Phase 16 has zero external HTTP calls so this is moot — but the discipline must propagate to any future SPC-related extension.
- **`@grok.decorators.func`** for menu items (TS decorator pattern) OR the JSDoc-style metadata-comment pattern from `.claude/rules/function-metadata.md`. Existing `src/package.ts` already uses the decorator pattern with a polyfill; Phase 16 must match.
- **`findColumn` / `SEMTYPE.*` discipline** — never `df.col('Gene names')`; use `findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', ...])`. (Per-package CLAUDE.md.)
- **`grok.dapi.*` over raw HTTP** — runs.csv I/O uses `grok.dapi.files.*`, NEVER `fetch`. Per-package + root CLAUDE.md.
- **Tag namespace prefix:** every new tag MUST be under `proteomics.spc_*`. Verified non-colliding via grep against existing tag set.
- **`ensureFreshFloat` idempotency** when adding the `spc_metrics_meta` column (per-package CLAUDE.md).
- **No new `@datagrok/chem` dep** — Phase 16 has no Chem coupling.
- **Test conventions:** `grok test --host localhost`, file under `src/tests/`, imported from `src/package-test.ts`. Per `.claude/rules/testing.md`.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — every API verified against `js-api/src/` at HEAD and against in-repo precedents (qc-dashboard, enrichment-viewers, publishing/).
- Architecture: HIGH — sibling directories vs new directory resolved by CONTEXT.md + ARCHITECTURE.md research; file-layout contract explicit; tag namespace verified non-colliding.
- Pitfalls: HIGH for codebase-grounded (4 inherited from PITFALLS.md + 1 platform-API-deprecation surfaced in this research); MEDIUM for the platform-API speculations (A1/A2 in Assumptions Log).
- Validation Architecture: HIGH — every requirement mapped to a concrete test; commands runnable from a `grok test` shell; fixtures specified.

**Research date:** 2026-06-08
**Valid until:** 2026-07-08 (30 days; stack is stable v1.3-baseline; will revisit if datagrok-api ^1.25 → ^1.26 lands and deprecates `userDataStorage` more aggressively).
