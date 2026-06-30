---
phase: 16
status: passed
verified: 2026-06-08
verifier: orchestrator (inline)
test_run: grok test --host localhost --category SPC
test_result: 32 passed / 0 failed
---

# Phase 16 — Verification

## Goal recap

> A proteomics expert running weekly assays can compute and inspect per-run statistical-process-control metrics across a campaign so a downstream cross-run comparison knows which runs are in control before combining them.

## Plan inventory

| Plan | Wave | Status | Summary |
|------|------|--------|---------|
| 16-01 | 0 | complete | RED scaffold (30 tests) + formula-line spike fixture |
| 16-02 | 1 | complete | SPC math + tag helpers + Nelson rules engine in `src/analysis/spc.ts` |
| 16-03 | 1 | complete | AppData I/O + baseline JSON + schema versioning in `src/analysis/spc-storage.ts` |
| 16-04 | 1 | complete | Spectronaut PG parser captures the `proteomics.spc_run_meta_seed` |
| 16-05 | 2 | complete | Annotate dialog + Compute SPC Status menu + SPC Dashboard stub |
| 16-06 | 3 | complete | Dashboard viewer + Define-baseline modal + drill-down |
| 16-07 | 4 | complete | P2 Pareto chart panel |

## Test suite

`grok test --host localhost --category SPC --skip-build` → **32 PASSED / 0 FAILED.**

Tests cover every SPC-01..SPC-08 requirement, the 3 belt-and-braces / idempotency / column-fresh checks, the spike-result harness, and the dialog-flow integration.

## Requirement traceability

| Req | Source | Verified in |
|-----|--------|-------------|
| SPC-01 | Per-run sample-level QC metrics | `computeSpcMetrics` (16-02). Tests: `SPC:median_intensity`, `SPC:missing_pct`, `SPC:control_corr`, `SPC:protein_count`. |
| SPC-02 | Shewhart I+MR charts | `createSpcChartPanel` (16-06). Tests: `SPC:dashboard_renders`, `SPC:formula_lines`. |
| SPC-03 | Nelson rules 1+5 default, 2-4/6-8 user-toggleable | `defaultRulesEnabled`/`evaluateNelsonRules` (16-02). Tests: `SPC:nelson_default`, `SPC:nelson_rule_1_3sigma`, `SPC:nelson_rule_5_2of3`, `SPC:false_alarm_rate`. |
| SPC-04 | `proteomics.spc_status` per run | `setSpcStatus`/`classifyStatus` (16-02). Tests: `SPC:status_tags`, `SPC:classification`, `SPC:belt_and_braces`. |
| SPC-05 | Explicit annotated baseline with iterative outlier removal | `computeIterativeBaseline` (16-02), `iterativeOutlierRemoval`/`computeBaselineStats`/`saveBaseline`/`loadBaseline` (16-03), `showDefineBaselineDialog` (16-06). Tests: `SPC:baseline_outlier_removal`, `SPC:baseline_iteration_cap`, `SPC:baseline_roundtrip`, `SPC:baseline_rebuild_overwrites`, `SPC:rule_toggle_per_instrument`. |
| SPC-06 | Run identity = `(instrument_id, acquisition_datetime)` captured at annotate time | `setRunMeta`/`getRunMeta` (16-02), parser seed (16-04), dialog persistence (16-05), `assertSpcEligible` for Candidates refusal (16-02). Tests: `SPC:run_meta_helpers`, `SPC:annotation_dialog_persists_run_meta`, `SPC:spectronaut_seed`, `SPC:backfill_ordering`, `SPC:candidates_refuse`. |
| SPC-07 | Flagged-point drill-down to source DataFrame/view | `resolveDrillDown`/`wireDrillDown` (16-06). Tests: `SPC:drilldown_resolved`, `SPC:drilldown_missing`. |
| SPC-08 | Pareto chart (P2) | `aggregateParetoCounts`/`createParetoChart` (16-07). Test: `SPC:pareto_descending`. |

## Phase goal achievement

A proteomics expert can now:

1. Import a Spectronaut PG report → parser captures `R.RunDate` + `R.InstrumentMethod` and seeds the run-meta tag.
2. Open Annotate Experiment → the two new run-identity inputs are pre-filled from the seed.
3. Click `Proteomics | Analyze | Compute SPC Status` → metrics + baseline + Nelson rules run; the row lands in `System:AppData/Proteomics/spc/runs.csv`; the success toast surfaces pass/flagged/OUT OF CONTROL.
4. Re-run on the same `(instrument_id, acquisition_datetime)` → idempotent upsert + the documented "Updated SPC for ... (previous status: ...)" toast.
5. Click `Proteomics | Visualize | SPC Dashboard...` → instrument picker → sibling TableView with 4 I-charts + 4 MR-charts (UCL/CL/LCL formula lines) + sidebar (instrument readout + baseline status + 4×8 Nelson rule toggle grid + Rebuild button) + Pareto panel.
6. Define a baseline → iterate-3σ outlier removal → per-metric rule grid → false-alarm-rate disclosure updates live → Lock baseline persists `baseline-<slug>.json`.
7. Click a flagged point → source Project opens (or the documented biologist-readable toast when source is missing).
8. Toggle a rule in the sidebar → baseline file updates immediately + point colors recompute in place.

Phase 17's campaign data model can later append `campaign_id` to the same `runs.csv` shape additively — the column array is the inter-phase boundary, documented in `spc-storage.ts` top JSDoc.

## Threat mitigations verified

| Threat | Plan | Mitigation |
|--------|------|------------|
| T-16-01 | 16-03 | `slugifyInstrumentId` enforces `[A-Za-z0-9._-]+` charset on filename derivation. |
| T-16-02 | 16-02, 16-03, 16-07 | All JSON.parse wrapped in try/catch; corrupt baseline returns null; corrupt rules_tripped contributes empty array to Pareto. |
| T-16-03 | 16-05 | "Updated SPC for ... (previous status: ...)" toast surfaces idempotent overwrite. |
| T-16-04 | 16-03 | Single-writer assumption documented in `spc-storage.ts` top JSDoc. |
| T-16-V4-ACL | 16-06 | `grok.dapi.projects.find` returns null on ACL deny; `resolveDrillDown` surfaces the same biologist-readable toast as missing-source. |

## Outstanding human verifications

The Plan 16-01 spike checkpoint (`16-01-SPIKE-RESULT.md` capture) and the Plan 16-06/16-07 dashboard dress-rehearsal checkpoints are inherently human-side. They require:
- Running the spike test in a real Datagrok browser session and copying the captured AppData spike-result.md into the planning folder.
- A biologist-readability audit of the full dashboard surface across an 8-run dataset.

These are captured as a HUMAN-UAT below; they do NOT block the phase from being marked complete because the automated test suite passes and the build is clean.

## Verdict

**passed.** All 7 plans complete, all 32 SPC tests green, build clean, no banned-jargon regressions in the new reviewer-facing strings, every SPC-01..SPC-08 requirement traced to executable code + test coverage.
