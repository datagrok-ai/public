---
feature: charts
target_layer: playwright
coverage_type: edge
priority: p1
realizes_atlas: [sunburst-unsupported-column-defense]
realizes: [charts.sunburst]
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

Bug-focused regression scenario for github-2954: Sunburst's column
selector accepts unsupported column types (e.g. date), and the viewer
used to crash during hierarchy conversion with a NullError. The fix
landed in Charts 1.21.0 — this scenario locks the regression bound:
configure the Sunburst hierarchy with a date column and verify either
the column selector pre-filters the unsupported type, or the viewer
handles it gracefully with no console error / crash. The canonical
`sunburst.md` exercises hierarchy configuration via the Select Columns
dialog but doesn't specifically test date-column rejection — that's
the gap this scenario closes.

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

- The github-2954 invariant is carried by Step 3: no console error
  during date-column hierarchy setup (console-error capture filters
  out benign network noise like 404s and favicon requests).
- Why `ae.csv`: the SDTM Adverse Events dataset has date-time columns
  (AESTDTC/AEENDTC/AEDTC) and is reachable at
  `System:AppData/Charts/ae.csv` (the `ApiSamples` package isn't
  deployed on dev, so `System:AppData/ApiSamples/ae.csv` isn't
  reachable there).

## Dataset metadata

```json
{
  "order": 33,
  "datasets": ["System:AppData/Charts/ae.csv"]
}
```
