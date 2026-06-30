# Phase 16: Sample-Level SPC Tracking - Context

**Gathered:** 2026-06-08
**Status:** Ready for planning

<domain>
## Phase Boundary

Per-run computation of four standard QC metrics (median log2 intensity, percent missing, control-replicate Pearson correlation, protein count) written onto each analyzed DataFrame as `proteomics.spc_metrics` + `proteomics.spc_status` tags, **plus** a multi-run SPC dashboard that reads an AppData-backed run log to render Shewhart I-MR charts + Nelson-rules-engine output across runs ordered by `(instrument_id, acquisition_datetime)`. Phase 16 stands alone — it ships a complete SPC story (Compute SPC Status + Dashboard + locked baseline + drill-down + Pareto chart) without depending on Phase 17 campaign storage; Phase 17 later decorates the same AppData log rows with `campaign_id` additively (no migration script). Run-level metrics ONLY — per-protein SPC is explicitly out of scope (Pitfall 5 lock).

Scope is fixed by REQUIREMENTS.md SPC-01..SPC-08 (P1 = SPC-01..SPC-07; P2 = SPC-08 Pareto). This discussion clarifies HOW to implement what's locked, not WHETHER to add new capabilities. Phase 17's campaign data model, Phase 18's cross-run comparison viewer, and Phase 15's reviewer-side publishing surfaces are explicit non-goals here.

**Audience pin (carried from Phase 15):** "Reviewer" / dashboard consumer at Cytokinetics demos may include biologists — the dashboard surface (banner copy, baseline-define modal, drill-down toasts, Pareto labels) gets the same jargon audit Phase 15 applied (banned words: `DataFrame`, `tag`, `semType`, `viewer factory`). The PROTEOMICS EXPERT is the primary operator (they own baseline definition, rule toggles, drill-down decisions) but a biologist-in-the-room must be able to read the dashboard cold.

</domain>

<decisions>
## Implementation Decisions

### Run identity capture (SPC-06, Pitfall 7)
- **D-01:** **Inline in `showAnnotationDialog` (Annotate Experiment), optional.** Two new inputs appended to the existing dialog: `instrument_id` (free-text string, default empty) and `acquisition_datetime` (DateTime input, default empty). Parser-side prefill when metadata is available: the Spectronaut PG-report parser extracts `R.RunDate` (first observed) and `R.InstrumentMethod` (first observed) and stores them on the parsed df as `proteomics.spc_run_meta_seed` (single JSON tag, format `{instrument_id, acquisition_datetime}`); the annotation dialog reads this seed and prefills the two inputs. MaxQuant / FragPipe / Generic / Spectronaut Candidates parsers leave the seed unset; the operator types both fields manually. Persisted on the df as `proteomics.spc_run_meta` (JSON `{instrument_id: string, acquisition_datetime: ISO-8601 string}`). Empty values are TOLERATED at annotation time so the v1.0–v1.3 path is not broken for users who skip SPC.
  - `Analyze → Compute SPC Status` REFUSES to run on a df with empty `spc_run_meta` and surfaces a single-line nudge: `"Open Annotate Experiment to set instrument + acquisition datetime."` Both fields are required for the SPC log row's primary key.
  - One df = one run = one `(instrument_id, acquisition_datetime)` pair. The four metrics aggregate across all samples in the df.
  - **Spectronaut Candidates is hard-refused** for Compute SPC Status — no per-sample intensities exist for the four metrics to compute against. Menu handler surfaces `"SPC requires per-sample intensities. Re-import this analysis from the Spectronaut PG report (not the Candidates report) to compute SPC."`

### Multi-run aggregation surface (SPC-02, SPC-05, SPC-07, Pitfall 6, Pitfall 8)
- **D-02:** **AppData append-only SPC run log: `System:AppData/Proteomics/spc/runs.csv`.** Compute SPC Status on a df appends one row to this CSV. Columns: `run_id` (UUID), `instrument_id`, `acquisition_datetime` (ISO-8601), `run_label` (df.name at compute time), `median_intensity`, `missing_pct`, `control_corr`, `protein_count`, `status`, `rules_tripped` (JSON array), `source_project_id` (null when source df is unpublished), `source_df_name` (df.name), `computed_at` (ISO-8601). Dashboard reads runs.csv, filters by instrument, and renders the charts.
  - **Phase 17 additivity contract:** Phase 17 reads the SAME runs.csv shape and adds a `campaign_id` column when saving a run to a campaign (back-filling existing rows with NULL is fine). The Phase 16 dashboard ignores `campaign_id` if present. No migration script — Phase 17's writer is the only thing that ever adds `campaign_id`.
  - **Storage growth (Pitfall 8):** runs.csv stores ONLY the per-run rollup (one row per run × 13 cols). 52-week × 2-instrument campaign → ~104 rows/year. Per-protein drill-down is NOT in this table.
  - **Re-running Compute SPC Status on the same df:** SECOND invocation overwrites the prior row matched by `(instrument_id, acquisition_datetime)` — idempotent on the run identity primary key. No duplicate-row creep.
  - **Schema versioning:** runs.csv first row column header is taken as authoritative; a one-time `spc.runs.schema_version = '1'` row in a sibling `runs-meta.json` lets future schema bumps detect old CSVs and migrate forward.
  - Locked baseline persistence (D-04 below) lives next to runs.csv in the same `System:AppData/Proteomics/spc/` directory.

### Control-replicate Pearson correlation definition (SPC-01)
- **D-03:** **Mean of pairwise Pearson correlation across Group 1 (Control) intensity columns.** With Group 1 columns `[c_1, c_2, ..., c_k]`, compute `mean(Pearson(c_i, c_j))` for all unordered pairs `(i, j)` where `i < j`. Reuses existing `proteomics.groups` annotation — zero new column concept, zero new annotation step, works on every v1.0–v1.3 dataset. NaN-tolerant pairwise: pairs with fewer than 3 jointly-observed proteins fall back to NaN and are excluded from the mean. If fewer than 2 control columns exist, `control_corr = NaN` and the rules engine treats NaN as "metric unavailable" (skip, not flag).
  - Pooled QC sample as a typed first-class concept (option B from discussion) is deferred to v1.5 — the closest enabling step is Phase 16's `Annotate Experiment` UX, but adding an optional pooled-QC input now risks scope creep and conflicts with the "minimal Annotate UX change" intent of D-01.
  - **Documented limitation:** this metric measures within-run reproducibility of biological controls; it does NOT detect inter-run batch drift (Pitfall 12). The dashboard's pre-demo dress-rehearsal copy must clarify this so a biologist reviewer doesn't over-read a high control-corr as "no batch drift".

### Baseline definition + Nelson rule toggling + persistence (SPC-03, SPC-05, SPC-08, Pitfall 6)
- **D-04:** **SPC Dashboard is the single surface; baseline + rule toggles managed in-dashboard sidebar.** First open of `Visualize → SPC Dashboard...` with no baseline locked: dashboard renders an empty-state banner `"No baseline locked for instrument <N>. [Define baseline...]"`. Click → modal listing all runs for that instrument from `runs.csv` with checkboxes (default: first 7 runs checked; user can adjust). Below the checkboxes: a "Iterate 3σ outlier removal" toggle (default ON; capped at 2 iterations per Pitfall 6) and a per-metric Nelson rule toggle grid (rules 1+5 default ON; 2-4 / 6-8 default OFF). On OK: compute mean+SD per metric over checked-and-non-outlier runs, write `System:AppData/Proteomics/spc/baseline-<instrument_id>.json` with schema:
  ```json
  {
    "instrument_id": "QExactive-01",
    "locked_at": "2026-07-10T...",
    "included_run_ids": ["uuid-a", "uuid-b", ...],
    "excluded_run_ids": ["uuid-x"],
    "iteration_trace": [{"iter": 1, "dropped": ["uuid-x"]}],
    "metrics": {
      "median_intensity": {"mean": 22.31, "sd": 0.42},
      "missing_pct":      {"mean": 12.6,  "sd": 1.1},
      "control_corr":     {"mean": 0.94,  "sd": 0.018},
      "protein_count":    {"mean": 4820,  "sd": 145}
    },
    "rules_enabled": {
      "median_intensity": {"nelson_1": true, "nelson_5": true, ...},
      "missing_pct":      {"nelson_1": true, "nelson_5": true, ...},
      "control_corr":     {"nelson_1": true, "nelson_5": true, ...},
      "protein_count":    {"nelson_1": true, "nelson_5": true, ...}
    }
  }
  ```
  - **Per-instrument baselines.** One baseline.json per `instrument_id`. Dashboard's instrument-picker shows one entry per `(instrument_id, has-baseline?)` pair; switching instrument switches baseline.
  - **Rebuild baseline:** sidebar button "Rebuild baseline..." opens the same modal pre-populated with current selection. On OK, OVERWRITES baseline-<instrument_id>.json (no version history kept — keeping baseline rolls is v1.5+).
  - **Rule toggle:** changes to the per-metric rule grid in the sidebar live-update the chart highlights AND rewrite `rules_enabled` in baseline.json. No separate "save settings" step — toggle = persisted.
  - **Pass / flagged / out_of_control classification (SPC-04, Claude's discretion below):** `pass` = no enabled rule trips; `out_of_control` = Nelson rule 1 (single point beyond 3σ) trips; `flagged` = any other enabled rule trips. This collapses rule-1 as the "hard limit" (matches operator mental model that 3σ excursion is structurally different from a 2-of-3-beyond-2σ pattern). Encoded into the rules-engine return value.
  - **Drill-down on flagged point:** click a flagged point on the I-chart → handler looks up the runs.csv row's `source_project_id`; if non-null AND project still exists, opens that project via `grok.dapi.projects.find(id).open()`. Otherwise surfaces a toast: `"Source for '<run_label>' is not currently available. Open the run's analyzed file and re-run Compute SPC Status."` No automatic file reopen.
  - **Pareto chart (SPC-08):** rendered as a second tab/panel inside the dashboard. Reads each row's `rules_tripped` JSON, counts per (rule, metric), bar-charts descending. Phase 16 P2 — ship if budget allows; the rules-tripped tag is set unconditionally so a deferred Pareto can land in v1.4.1 without re-instrumenting metric code.

### Claude's Discretion
- **`proteomics.spc_metrics` tag shape on the analyzed df.** JSON `{median_intensity: number, missing_pct: number, control_corr: number, protein_count: number, sample_count: number, computed_at: ISO-8601}` — matches D-02 column shape minus the run-identity fields (which live in `spc_run_meta`). Belt-and-braces column encoding follows Phase 15 D-05 pattern: write the same six values into a single-row `spc_metrics_meta` column so a published df's metrics survive Project serialization. Phase 15 publishing helper does NOT exclude SPC tags by default — published projects that carry an SPC-annotated run can surface "this run was QC'd `pass` on 2026-07-10" via the reviewer panel.
- **`proteomics.spc_status` tag values.** Lowercase strings `'pass'` / `'flagged'` / `'out_of_control'`. New SEMTYPE addition deferred to Phase 17 (CAMP-07 in REQUIREMENTS.md already plans `Proteomics-SpcStatus`); Phase 16 just sets the tag and the dashboard reads it directly without semType-based detection.
- **`proteomics.spc_rules_tripped` tag shape.** JSON array of canonical rule names: `['nelson_1@median_intensity', 'nelson_5@control_corr', ...]`. Empty array for `pass`. Used by the Pareto chart AND by drill-down to highlight which axis explains the flag.
- **Protein count above threshold (SPC-01 wording).** Defined as `# proteins with non-null intensity in at least one sample column` — i.e. row count minus all-missing-rows. No user-tunable threshold for v1.4; a configurable minimum-quantified-samples knob is v1.5+ scope (deferred).
- **Run identity coverage on Spectronaut PG-report imports.** `R.RunDate` parsed as ISO-8601; `R.InstrumentMethod` value used verbatim as `instrument_id` until/unless the user overrides it inline. Spectronaut row-per-(protein, sample) means many rows carry the same R.RunDate; parser picks the EARLIEST observed.
- **Compute SPC Status menu position.** `Analyze → Compute SPC Status` (immediate, no `...` suffix per CLAUDE.md convention — uses the existing `proteomics.groups` annotation + `spc_run_meta` annotation; no parameterization needed at the menu surface). Surfacing a `TaskBarProgressIndicator` for the (small) compute is good practice but optional.
- **SPC Dashboard menu position.** `Visualize → SPC Dashboard...` (with `...` because the entry opens a dialog/picker for instrument selection on first invocation; subsequent invocations remember the last-picked instrument in `userDataStorage` per Phase 13 D-02 precedent).
- **Dashboard implementation primitive.** The 4-up I-chart panel uses `DG.Viewer.lineChart` against a per-instrument filtered slice of runs.csv (re-read on dashboard open + on rebuild-baseline). Formula lines render UCL = mean+3SD, CL = mean, LCL = mean-3SD per metric. MR-chart is a sibling line chart against `abs(diff(metric))` columns computed at dashboard-open time from the same slice. Pareto is a `DG.Viewer.barChart`.
- **runs.csv read concurrency.** Single-writer assumption (one Datagrok user per machine appending at a time). If concurrent writes become a Cytokinetics reality, Phase 16 ships a small file-lock (`runs.csv.lock` sentinel) and the writer retries after 100ms — but this is deferred unless a research-time read of how Datagrok serializes AppData writes makes it necessary.

### Folded Todos
None. The only related backlog item is Phase 999.3 (publishing) which was already promoted into Phase 15. No standalone SPC todos exist in `.planning/todos/pending/`.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase 16 requirements + scope (LOCKED)
- `.planning/REQUIREMENTS.md` SPC-01..SPC-08 — the locked requirements for this phase. P1 = SPC-01..SPC-07; P2 = SPC-08 Pareto chart. Do NOT discuss or replan requirement WHATs.
- `.planning/ROADMAP.md` Phase 16 entry — phase goal + 5 numbered success criteria.
- `.planning/PROJECT.md` v1.4 Cross-Team Review section — milestone-level positioning + Cytokinetics audience phrasing.

### Research synthesis (load-bearing for plan-time research)
- `.planning/research/SUMMARY.md` §"Phase 16: Sample-Level SPC Tracking" + §"Recommended Stack" SPC bullet — affirms in-package ~390 LOC TS SPC, no `@datagrok-libraries/statistics` primitive, public-domain math (Nelson 1984; Wheeler's d2=1.128, D4=3.267 for n=2 MR). Phase 16 has NO required plan-time research — math is fully scoped here.
- `.planning/research/ARCHITECTURE.md` §3 "Tag Contract Extensions" (3 new tags: `spc_status`, `spc_metrics`, `spc_rules_tripped`) + §4 "Directory Structure" (math in `src/analysis/spc.ts`, viewer in `src/viewers/spc-dashboard.ts` — NOT a new `src/spc/` directory; resolves STACK vs ARCHITECTURE divergence per SUMMARY).
- `.planning/research/PITFALLS.md` §"Pitfall 5" (per-protein alarm flood — D-domain run-level lock) / §"Pitfall 6" (baseline contamination — D-04 iterative outlier removal + explicit lock) / §"Pitfall 7" (run-order ambiguity — D-01 instrument + datetime capture) / §"Pitfall 8" (storage bloat — D-02 rollup-only). All four pitfall mitigations baked into the decisions above; planner verifies each in the plan's task structure.
- `.planning/research/STACK.md` §"In-package SPC" — Shewhart I-MR + Nelson rules 1–8; ~390 LOC TS; no usable npm; R `qcc` fallback PARKED for v1.5+ per REQUIREMENTS.md "Future Requirements".

### Phase 15 carry-forward (do NOT regress)
- `.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md` — audience pin (biologist-in-the-room demo readiness), belt-and-braces tag+column encoding pattern (Phase 16 reuses for `spc_metrics_meta` column), `userDataStorage` cross-session precedent (D-04 dashboard remembers last-picked instrument).
- `.planning/codebase/ARCHITECTURE.md` §"In-place Mutation" + §"Cross-DataFrame Selection" + §"Anti-Patterns" — confirms why SPC math mutates the source df (tags + spc_metrics_meta column written in place); confirms the drill-down click-handler uses single-producer subscription (NOT `CampaignSelectionBus` — that's Phase 18).
- `.planning/codebase/STRUCTURE.md` §"Where to Add New Code" — file-layout contract (math under `src/analysis/`; viewer under `src/viewers/`; tests under `src/tests/`).
- `.planning/codebase/CONVENTIONS.md` — `findColumn` / `SEMTYPE` discipline carries forward (locate Group 1 columns via `proteomics.groups` JSON, not by `df.col('Group 1')`).

### Datagrok platform reference (verify-at-plan-time)
- `grok.dapi.files.readAsText` / `writeAsText` for `System:AppData/Proteomics/spc/runs.csv` append + read. HitTriage precedent at `packages/HitTriage/src/app/hit-triage-app.ts:375` (also called out in Phase 15 canonical_refs).
- `grok.dapi.projects.find(id).open()` for D-04 drill-down click-handler.
- `grok.userSettings` / `grok.dapi.userDataStorage` for D-04 last-picked-instrument cache (Phase 13 D-02 precedent — confirmed pattern in `13-CONTEXT.md`).
- `DG.Viewer.lineChart` + formula-lines API for D-04 I-chart UCL/CL/LCL rendering. Phase 14 D-04 magenta/cyan/gray color lock does NOT apply here (SPC charts are univariate trend, not bivariate volcano) — use platform default colors.
- `DG.Viewer.barChart` for D-04 Pareto.

### Datagrok platform docs (background)
- `packages/Proteomics/CLAUDE.md` — pipeline tag conventions, `findColumn`/`SEMTYPE` requirements, menu-suffix convention (`...` = dialog; bare = immediate), `ensureFreshFloat` idempotency precedent for any computed columns Phase 16 adds to the source df.

### External SPC literature (verified by research)
- Nelson 1984 — Nelson rules 1–8; false-alarm rate per rule (Nelson, L. S. (1984). Technical Aids. Journal of Quality Technology, 16(4), 238–239).
- Wheeler I-MR chart constants: d2=1.128, D4=3.267 for n=2 (Wheeler, D. J. & Chambers, D. S. (1992). Understanding Statistical Process Control, 2nd ed.). Used in the MR-chart's UCL = D4 × mean(MR) computation.
- SUMMARY.md cites these; planner doesn't need to re-verify.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`src/analysis/experiment-setup.ts:63` `showAnnotationDialog`** — D-01 extension point. Two new optional inputs (`instrument_id` string, `acquisition_datetime` DateTime) appended after Group 2 columns. Existing `seedAnnotationDialogInputs` widened to also read `proteomics.spc_run_meta_seed` for parser prefill.
- **`src/analysis/experiment-setup.ts:13-26` `setGroups`/`getGroups` JSON-tag pattern** — exact precedent for new `setRunMeta(df, meta)` / `getRunMeta(df)` helpers under `src/analysis/spc.ts`.
- **`src/utils/column-detection.ts findColumn`** — D-03 uses this to locate Group 1 columns via the parsed `proteomics.groups` JSON. No raw `df.col('Control')` calls.
- **`src/viewers/qc-computations.ts:27 ensureFreshFloat`** — idempotency precedent. Compute SPC Status writes the `spc_metrics_meta` belt-and-braces column via this pattern (remove-if-exists then add) so re-running is clean.
- **`src/viewers/qc-computations.ts` median / NaN-tolerant aggregation patterns** — the four SPC metric computations follow the same NaN-tolerant pattern (treat `col.isNone(i)` as missing; skip cells, not whole columns).
- **`src/viewers/qc-dashboard.ts:1-30` dashboard composition pattern** — Phase 16 SPC dashboard mirrors this shape (multi-viewer TableView with sidebar controls); Phase 7's per-sample QC dashboard is the SPC dashboard's closest analog (4 panels, formula lines, click-through).
- **`grok.shell.addTableView(...)` pattern** in `src/analysis/enrichment.ts:574` and `src/viewers/enrichment-viewers.ts:191` — D-04 dashboard opens the runs.csv-derived df in its own TableView (sibling-DataFrame rule from ARCHITECTURE.md §7); follows the v1.3 enrichment precedent.
- **`grok.dapi.projects.find(id).open()` pattern** in v1.3 codebase (and the Phase 15 publish-roundtrip test) — D-04 drill-down handler reuses.
- **Phase 15 belt-and-braces column encoding** (`src/publishing/trim-dataframe.ts` if it exists post-Phase-15, or the Phase 15 plan-15-02 spec for it) — D-04 / Claude's-discretion `spc_metrics_meta` column reuses this pattern.

### Established Patterns
- **`proteomics.*` tag namespace for workflow state.** Phase 16 adds 3 new tags on the analyzed df: `proteomics.spc_status` (lowercase string), `proteomics.spc_metrics` (JSON), `proteomics.spc_rules_tripped` (JSON array). PLUS `proteomics.spc_run_meta` (JSON, set by D-01 at Annotate Experiment) and `proteomics.spc_run_meta_seed` (JSON, set by Spectronaut PG parser pre-annotation). All five verified non-colliding via grep against existing tag set (`source`, `groups`, `normalized`, `imputed`, `preNormalized`, `de_complete`, `de_method`, `enrichment*`, `location_acc_hash`, plus the Phase 15 `published*` set).
- **`@grok.decorators.func` for menu items** in `src/package.ts`. New menu entries: `Proteomics | Analyze | Compute SPC Status` (immediate) + `Proteomics | Visualize | SPC Dashboard...` (dialog suffix). Both gate via tag preconditions per ARCHITECTURE.md §6 — Compute SPC Status requires non-empty `spc_run_meta`; SPC Dashboard requires `runs.csv` to exist (else surfaces empty-state banner).
- **Three-level fallback pattern.** Not applicable to SPC math — it's pure TS and deterministic; no R cold-start path. The R `qcc` fallback is REQUIREMENTS.md "Future Requirements" deferred.
- **`ensureFreshFloat` idempotency.** Compute SPC Status re-runs on the same df must produce the same `spc_metrics_meta` column (no duplicate columns). Same pattern, applied to all six values' single-row columns.

### Integration Points
- **`src/package.ts`** — register two new menu items + ensure Compute SPC Status's gate-check uses the new `spc_run_meta` tag.
- **`src/analysis/experiment-setup.ts`** — D-01 modifies in place (add two inputs, two getters, prefill from seed).
- **`src/analysis/spc.ts` (NEW)** — owns the math: `computeSpcMetrics(df, groups, runMeta)`, `setSpcStatus(df, metrics, ruleResult)`, `getSpcStatus(df)`, `evaluateNelsonRules(seriesValue, baselineMean, baselineSd, enabledRules) → {status, rulesTripped}`. Per ARCHITECTURE.md §4: math in `analysis/`, NOT a new `src/spc/` directory.
- **`src/analysis/spc-storage.ts` (NEW)** — owns AppData I/O: `appendRun(row)`, `loadRuns(instrumentId?)`, `loadBaseline(instrumentId)`, `saveBaseline(instrumentId, baseline)`. Sibling to `spc.ts` to keep math pure (no AppData calls from inside `computeSpcMetrics`).
- **`src/viewers/spc-dashboard.ts` (NEW)** — owns the viewer factory + sidebar UI: `openSpcDashboard()` (opens TableView), `createSpcChartPanel(df, baseline, metric)` (renders one I-chart + MR-chart pair with formula lines), `showDefineBaselineDialog(instrumentId, runs)` (modal), `createParetoChart(df)` (P2). Sidebar contains instrument picker + baseline status + per-metric rule toggle grid.
- **`src/parsers/spectronaut-parser.ts`** — Spectronaut PG report parser augmentation to set `proteomics.spc_run_meta_seed` from earliest-observed `R.RunDate` + first-observed `R.InstrumentMethod`. The streaming PG variant (Phase 12) must surface these columns BEFORE filtering — verify they're not dropped by the existing column filter.
- **`src/tests/spc.ts` (NEW)** — unit tests per Pitfall 5/6/7/8 levers (Nelson rule false-alarm rate, baseline iterative outlier removal, backfill-acquisition-datetime sorting, runs.csv row count growth bound).
- **`detectors.js`** — no Phase 16 SEMTYPE additions; `Proteomics-SpcStatus` is deferred to Phase 17 per CAMP-07.

</code_context>

<specifics>
## Specific Ideas

- **AppData log shape is the Phase 17 contract.** `System:AppData/Proteomics/spc/runs.csv` is the file Phase 17 reads from and writes to (with the additional `campaign_id` column). The runs.csv column order + the JSON shape of baseline-<instrument_id>.json IS the inter-phase boundary. Planner must record the columns + schema_version in `spc-storage.ts` JSDoc + reference it from this CONTEXT in any Phase 17 plan.
- **Belt-and-braces (Phase 15 D-05 / Pitfall 3) applies here.** The `spc_metrics_meta` single-row column on each df is required because the SPC tags MIGHT travel into a Phase 15 published Project; the column survives serializer-strip better than the tag (Phase 13 e527d07ba1 evidence).
- **Pareto chart (SPC-08, P2) is genuinely optional** — the rules-tripped tag is set unconditionally so the Pareto can land in a v1.4.1 hotfix if Phase 16 budget runs short. Planner: tag-set is non-negotiable; Pareto VIEWER is conditional.
- **Demo audience model carries from Phase 15.** Cytokinetics demo audience may include biologists. Every reviewer-facing surface (the dashboard's empty-state banner, the baseline-define modal's copy, the drill-down toast, the rule-toggle labels) gets a pre-demo dress rehearsal with a biologist-class user per Pitfall 14 recovery strategy.
- **Nelson rule defaults are LOAD-BEARING (SPC-03).** Rules 1 + 5 ON; rules 2-4 / 6-8 OFF by default. Combined false-alarm rate at this default is ~0.4% (well below the 2.65% all-8 catastrophe per Pitfall 5). The dashboard UI shows the running false-alarm-rate estimate in a tooltip when the user toggles additional rules ON — `"Enabling rules 1+5+2 raises the expected false-alarm rate to ~1.1%."` — this is the disclosure mechanism from Pitfall 5's test/UAT lever.

</specifics>

<deferred>
## Deferred Ideas

- **Pooled-QC sample column as a first-class concept** (option B from the control-corr discussion) — D-03 uses Group 1 pairwise instead. REQUIREMENTS.md "Future Requirements" already defers per-protein curated panel SPC; a pooled-QC story is a sibling deferral; revisit if Cytokinetics surfaces inter-run batch drift as a recurring complaint.
- **Optional R `qcc` SPC fallback path** — REQUIREMENTS.md "Future Requirements" deferred to v1.5+. TS-only client path is the v1.4 scope.
- **CUSUM / EWMA SPC charts** — REQUIREMENTS.md "Future Requirements" deferred. Shewhart I-MR + Nelson rules cover the Cytokinetics ask; CUSUM/EWMA are more-sensitive variants with their own tuning knobs.
- **Per-protein curated panel SPC** — REQUIREMENTS.md "Future Requirements" deferred. Pitfall 5 explicitly bans per-protein SPC at this scope.
- **`Proteomics-SpcStatus` SEMTYPE registration in `src/utils/proteomics-types.ts` + `detectors.js`** — Phase 17 CAMP-07 owns this; Phase 16 sets `proteomics.spc_status` as a raw tag without a SEMTYPE.
- **Cross-instrument comparison view ("compare instrument-1 to instrument-2 trend")** — out of v1.4 scope. Phase 16 ships one chart panel per instrument; switching instrument is a sidebar action.
- **Publish-block-on-QC-fail policy (gating Phase 15 publish OK on `proteomics.spc_status`)** — REQUIREMENTS.md "Future Requirements" deferred. Phase 16 sets the tag; Phase 15 publishing path does NOT gate on it; v1.5+ may introduce the gate.
- **Email / Slack / webhook alerting when a run flags** — REQUIREMENTS.md "Out of Scope" deferred. Operator manually checks the dashboard.
- **Baseline version history / undo "Rebuild baseline"** — D-04 rebuild OVERWRITES the prior baseline.json (one file per instrument). Versioned baselines / undo are v1.5+ scope.
- **runs.csv concurrent-writer file lock** — single-writer assumption in v1.4; ship a `runs.csv.lock` sentinel only if research-time read of Datagrok AppData write semantics shows concurrent corruption risk.
- **Per-protein drill-down from a flagged SPC point** — D-04 drill-down opens the source PROJECT (analyzed run); per-protein investigation lives inside that project's standard viewers. No new SPC-specific per-protein viewer.

</deferred>

---

*Phase: 16-sample-level-spc-tracking*
*Context gathered: 2026-06-08*
