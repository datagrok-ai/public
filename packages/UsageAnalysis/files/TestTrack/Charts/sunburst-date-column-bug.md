---
feature: charts
sub_features_covered:
  - charts.sunburst
  - charts.echart-base
target_layer: playwright
coverage_type: edge
pyramid_layer: bug-focused
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst-date-column-bug.md
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
related_bugs:
  - github-2954
---

# Sunburst viewer — date-column hierarchy graceful handling (github-2954)

Bug-focused regression scenario for github-2954: Sunburst's column selector
accepts unsupported column types (e.g., date), and the viewer crashes
during `toForest`/`toTree` conversion with a NullError. The fix landed in
Charts 1.21.0 — this scenario locks the regression bound.

`pyramid_layer: bug-focused`. The canonical `sunburst.md` exercises
hierarchy configuration via Select Columns dialog but does NOT specifically
test date-column rejection — that's the bug-class invariant.

`related_bugs: [github-2954]` — bug-library reproduction class:
configure Sunburst hierarchy with a date column and verify either (a)
the column selector pre-filters the unsupported type, OR (b) the viewer
handles the column gracefully with no console error / crash.

## Setup

A clean Datagrok session. Single scenario; uses a dataset with date columns
(SPGI.csv has date columns; ae.csv has AESTDTC/AEENDTC date-time columns).

## Scenarios

### Scenario 1: Sunburst hierarchy with date column does not crash (github-2954 invariant)

Steps:

1. Open `System:AppData/Charts/ae.csv` (143-row SDTM Adverse Events shape;
   has `AESTDTC` / `AEENDTC` / `AEDTC` date-time columns plus `AESTDY` /
   `AEENDY` numeric date-day columns).
2. Add a Sunburst viewer via `tv.addViewer('Sunburst')`. Wait 3000ms.
   **Expected:** Sunburst attached; `tv.viewers` includes Sunburst.
3. Set `hierarchyColumnNames` to include a date column:
   `sunburst.setOptions({hierarchyColumnNames: ['AESTDTC']})`.
   **Expected (github-2954 invariant):** NO console error / NullError /
   `toString on null` / `toForest`/`toTree` crash. The viewer either:
   - (a) ignores the unsupported column (pre-filter behavior), OR
   - (b) renders gracefully (no exception bubbling up to console).
4. Read back `sunburst.props.get('hierarchyColumnNames')` (race-tolerant
   try/catch). Log the result for diagnostic visibility.
5. Verify the viewer's root DOM element (`sunburst.root`) is non-empty
   and has non-zero size (visual stability invariant).
6. Set `hierarchyColumnNames` to a known-good string column for cleanup
   normalization: `sunburst.setOptions({hierarchyColumnNames: ['AETERM']})`.
   **Expected:** the viewer recovers cleanly.

## Notes

- **github-2954 invariant carrier:** Step 3 — no console error during
  date-column hierarchy setup. Console-error capture filters benign
  network noise (404, Failed to load resource, favicon).
- **Why ae.csv:** SDTM Adverse Events shape has date-time columns
  (AESTDTC/AEENDTC/AEDTC) and is reachable on dev (MCP confirmed
  `System:AppData/Charts/ae.csv` exists; ApiSamples is NOT deployed on
  dev so System:AppData/ApiSamples/ae.csv is unreachable).
- **Authority:** atlas-driven; closes the bug coverage gap surfaced in
  chain rev 2 `bug_match_attempts_skipped` for github-2954
  (skip_category: below_trigger_threshold; sunburst.md weakly matches
  hierarchy step but doesn't specifically test date-column type).
- **Helpers used:** `softStep`, `loginToDatagrok`, `specTestOptions`.

## Dataset metadata

```json
{
  "order": 33,
  "datasets": ["System:AppData/Charts/ae.csv"]
}
```
