# Phase 16: Sample-Level SPC Tracking - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-08
**Phase:** 16-sample-level-spc-tracking
**Areas discussed:** Run identity capture entry point, Multi-run aggregation surface, Control-replicate Pearson correlation definition, Baseline + Nelson-rule operator workflow

---

## Run identity capture entry point

| Option | Description | Selected |
|--------|-------------|----------|
| Inline in Annotate Experiment, optional | Extend existing dialog with two optional inputs prefilled from parser metadata where available (Spectronaut R.RunDate / R.InstrumentMethod). Empty values OK at annotate time; Compute SPC Status refuses missing identity and surfaces re-annotate nudge. Non-SPC users can ignore the extra fields. | ✓ |
| Separate "Annotate Run" dialog, on-demand | New menu item under Annotate (or auto-opens on first Compute SPC Status). Captures fields + optional control-replicate column. Annotate Experiment stays unchanged. Mild SPC-06 divergence ("at annotate-experiment time") but cleaner separation. | |
| Auto-extract at parse time + confirmation on Compute SPC Status | Spectronaut parsers populate seed from R.RunDate/R.InstrumentMethod; others leave blank. Compute SPC Status pops a confirmation/edit dialog (prefilled when possible) before computing. Highest fit for SPC-06 + Pitfall 7 but most steps. | |

**User's choice:** Inline in Annotate Experiment, optional
**Notes:** Becomes D-01. One df = one run = one (instrument_id, acquisition_datetime). Spectronaut Candidates parser is hard-refused for Compute SPC Status because the report has no per-sample intensities to compute metrics against.

---

## Multi-run aggregation surface (Phase 16 without Phase 17)

| Option | Description | Selected |
|--------|-------------|----------|
| AppData SPC log (CSV, append-only) | One row per Compute SPC Status invocation written to System:AppData/Proteomics/spc/runs.csv. Dashboard reads + filters by instrument. Phase 17 decorates rows with campaign_id additively — no migration. ~104 rows/year per 2 instruments. | ✓ |
| userDataStorage SPC log (per-user, JSON) | Same shape via grok.dapi.userDataStorage (Phase 13 D-02 precedent). Per-user cross-session cache. Per-analyst silos. | |
| Shell aggregation only (no persistence) + warning | Dashboard scans grok.shell.tables for dfs with proteomics.spc_metrics. Baseline lives in userDataStorage. Loses data not currently open. | |
| Compose with Phase 15 published projects | Operator publishes each run; SPC dashboard reads published projects under a chosen Space. Couples SPC to publishing UX — contradicts audience model (publish = biologist; SPC = proteomics expert). | |

**User's choice:** AppData SPC log (CSV, append-only)
**Notes:** Becomes D-02. runs.csv shape + baseline-<instrument_id>.json shape are the Phase 17 contract — additivity, no migration script. Re-running Compute SPC Status on the same (instrument_id, acquisition_datetime) pair overwrites the prior row.

---

## Control-replicate Pearson correlation definition

| Option | Description | Selected |
|--------|-------------|----------|
| Pairwise mean across Group 1 columns | Mean of Pearson(c_i, c_j) for all unordered pairs in Group 1 (Control). Reuses existing group model, zero new annotation, works on every v1.0–v1.3 dataset. Measures within-run reproducibility but misses inter-run batch drift. | ✓ |
| Explicit pooled-QC column (new annotation) | Annotate Experiment gains optional third input "Pooled QC sample" (single column). Control-corr = Pearson(QC, mean(Group 1)). Industry-standard LC-MS batch-effect detection (Pitfall 12). Adds an annotation step; many datasets won't have QC sample. | |
| Pairwise mean across all intensity columns | Pearson across every (log2 intensity) column pair, mean. Conflates within-group + between-group correlation. | |
| Hybrid: pooled-QC if annotated, else Group 1 fallback | Annotate gains an optional Pooled QC input. Filled → Pearson(QC, mean(Group 1)); empty → fallback to pairwise Group 1. Metric definition recorded in spc_metrics so dashboard tracks mode per run. | |

**User's choice:** Pairwise mean across Group 1 columns
**Notes:** Becomes D-03. Documented limitation: this metric does NOT detect inter-run batch drift. Dashboard pre-demo dress rehearsal must surface this so biologists don't over-read a high control-corr as "no batch drift". Pooled QC sample deferred to v1.5.

---

## Baseline + Nelson-rule operator workflow

| Option | Description | Selected |
|--------|-------------|----------|
| Embedded in SPC Dashboard + AppData JSON | Dashboard is single surface. Empty-state banner "No baseline locked. [Define baseline…]". Modal lists runs with checkboxes + iterate 3σ toggle (capped 2 iters). Per-metric Nelson rule toggle grid in sidebar (rules 1+5 default ON). Locked baseline persists as System:AppData/Proteomics/spc/baseline-<instrument_id>.json (mean+SD per metric, included run_ids, rule toggles, iteration trace). Rebuild reopens modal pre-populated. | ✓ |
| Separate menu items | Three menu items: Visualize → SPC Dashboard..., Analyze → Define SPC Baseline..., Analyze → Configure Nelson Rules.... Dashboard read-only. More menu surface; cleaner separation. | |
| Auto-prompt on Nth run + per-baseline persisted toggles | After log has ≥7 runs for one instrument, next Compute SPC Status pops "Enough runs for a baseline. Define now?" → modal. Dashboard renders only after baseline exists. Most opinionated; least menu surface. | |
| Inline edit in dashboard sidebar + userDataStorage | Same UX as option A but grok.dapi.userDataStorage instead of AppData JSON. Per-user-per-machine — each analyst maintains own baseline. Less shareable across Cytokinetics team. | |

**User's choice:** Embedded in SPC Dashboard + AppData JSON
**Notes:** Becomes D-04. Per-instrument baselines (one baseline.json per instrument_id). Rebuild OVERWRITES (no version history — v1.5+). Drill-down opens source_project_id via grok.dapi.projects.find(id).open(); falls back to toast when project absent. Pareto chart is P2 — rules_tripped tag is non-negotiable, viewer is conditional on budget.

---

## Claude's Discretion

The following items were left to Claude's discretion in CONTEXT.md `<decisions>` — `Claude's Discretion`:

- `proteomics.spc_metrics` tag shape on the analyzed df (JSON with the four metrics + `sample_count` + `computed_at`) and the belt-and-braces `spc_metrics_meta` column mirror.
- `proteomics.spc_status` tag values: lowercase `'pass'` / `'flagged'` / `'out_of_control'` (no SEMTYPE — `Proteomics-SpcStatus` is Phase 17 CAMP-07).
- `proteomics.spc_rules_tripped` JSON-array canonical rule-name shape (e.g. `'nelson_1@median_intensity'`).
- "Protein count above threshold" defined as `# proteins with non-null intensity in at least one sample column`. No tunable threshold for v1.4.
- Spectronaut PG-report run-identity coverage: earliest-observed `R.RunDate`, first-observed `R.InstrumentMethod`; row-per-(protein, sample) shape collapsed at parse time.
- Menu positions: `Analyze → Compute SPC Status` (immediate); `Visualize → SPC Dashboard...` (dialog suffix; last-picked instrument cached in `userDataStorage`).
- Dashboard implementation primitive: `DG.Viewer.lineChart` for I-chart + MR-chart with formula lines; `DG.Viewer.barChart` for Pareto.
- runs.csv read/write single-writer assumption (lock-file deferred unless research surfaces a Datagrok AppData-concurrency risk).
- Pass / flagged / out_of_control classification rule: Nelson rule 1 (3σ) trip → `out_of_control`; any other enabled rule trip → `flagged`; none → `pass`.

## Deferred Ideas

- Pooled-QC sample as first-class concept (closest enabling step is D-01 Annotate Experiment UX; revisit if Cytokinetics surfaces inter-run batch drift).
- Optional R `qcc` SPC fallback path (REQUIREMENTS.md "Future Requirements", v1.5+).
- CUSUM / EWMA SPC charts.
- Per-protein curated panel SPC (Pitfall 5 explicit ban at this scope).
- `Proteomics-SpcStatus` SEMTYPE registration (Phase 17 CAMP-07 owns this).
- Cross-instrument comparison view ("compare instrument-1 to instrument-2 trend").
- Publish-block-on-QC-fail policy (REQUIREMENTS.md "Future Requirements").
- Email / Slack / webhook alerting on flag.
- Baseline version history / undo "Rebuild baseline".
- runs.csv concurrent-writer file lock.
- Per-protein drill-down from a flagged SPC point (drill-down opens source PROJECT only).
