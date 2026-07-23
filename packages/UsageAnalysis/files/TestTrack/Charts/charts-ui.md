---
feature: charts
target_layer: manual
coverage_type: manual
pyramid_layer: ui-only
ui_coverage_responsibility:
  - sunburst-multi-selection
  - sunburst-empty-category-click
  - tree-shift-click-multi-select
  - tree-shift-click-extend-selection
ui_coverage_delegated_to: null
produced_from: extracted
extracted_from:
  - public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst.md
  - public/packages/UsageAnalysis/files/TestTrack/Charts/tree.md
date_created: 2026-05-09
authored_by: charts-remediate-2026-05-09 ui-only extraction
related_bugs: []
---

# Charts — UI-only manual scenarios (canvas-gesture)

Manual-only scenarios for ECharts canvas-rendered gestures that have no
automatable equivalent exercising the same UI invariant. The JS-API
fallbacks (`df.selection.set`, `df.filter.set`) drive the underlying
bitsets but do not replicate the actual canvas-click semantics (segment
hit-testing, drill-down state, multi-select modifier-key handling, etc).

These scenarios require a human operator to drive the actual canvas
click events, and are documented here as the authoritative manual test
catalog. The companion playwright specs (`sunburst-spec.ts`,
`tree-spec.ts`) cover other Sunburst/Tree flows but explicitly do not
cover the canvas gestures listed below.

## Setup

A clean Datagrok session with the test user authenticated. Each
scenario block opens its own dataset.

## Scenarios

### Sunburst — multi-selection on segments (canvas, Click / Ctrl+Click / Ctrl+Shift+Click)

Source: extracted from `sunburst.md` `### Multi-selection behaviour`
section (canvas-rendered ECharts segments — no per-segment DOM, hit
tests must be done by mouse position on canvas).

Dataset: `System:DemoFiles/chem/SPGI.csv` (or `System:AppData/Chem/tests/spgi-100.csv`).

Steps:

1. Open `System:DemoFiles/chem/SPGI.csv` and add a Sunburst viewer.
2. Configure Sunburst hierarchy with at least 2 columns (e.g. via
   the Select Columns dialog: choose **Core** and **R101**).
3. **Click** a single Sunburst segment with the left mouse button.
4. **Verify:** exactly that segment is selected and the corresponding
   grid rows are selected (check the grid row selection count
   matches the segment's row count).
5. Hold **Ctrl** and click an additional segment.
6. **Verify:** both segments are selected; the grid row selection
   reflects the union.
7. Hold **Ctrl + Shift** and click one of the selected segments.
8. **Verify:** that segment is deselected; the grid row selection
   reflects the new union.

JS-API fallback note: `df.selection.setAll(false)` + per-row
`df.selection.set(idx, true)` drives the bitset BUT does NOT
exercise the canvas hit-test, segment-to-row mapping, or modifier-
key dispatching. This scenario must be run by hand.

### Sunburst — empty (null) category click (canvas)

Source: extracted from `sunburst.md` `### Select / filter on empty
category (SPGI_v2.csv)` section (clicking the grey null-segment is
canvas hit-testing on a non-string category bucket).

Dataset: `System:DemoFiles/chem/SPGI.csv`.

Steps:

1. Open SPGI_v2.csv and add a Sunburst viewer.
2. Configure a column known to contain nulls (e.g. **Sampling Time**)
   via Select Columns.
3. **Click** the null (grey) segment.
4. **Verify:** the segment behaves like any other category — the
   corresponding rows (those with null in that column) are selected
   or filtered (per the configured `onClick` action — Select or
   Filter).

JS-API fallback note: `df.selection.set(rowsWithNull, true)` drives
the result, but does NOT verify that the canvas grey-segment-click
is correctly hit-tested and routed to the null-category bucket.

### Tree — Shift+Click branches multi-selection (canvas)

Source: extracted from `tree.md` `### Collaborative filtering` step
1 (Tree branches are canvas inside `_echarts_instance_` — no
per-branch DOM exists; `Shift+Click` modifier dispatch is
non-synthesizable on canvas).

Dataset: `System:DemoFiles/demog.csv` (5850 rows). Hierarchy:
`['CONTROL', 'SEX', 'RACE']`.

Steps:

1. Open demog.csv and add a Tree viewer.
2. Set Tree hierarchy to **CONTROL → SEX → RACE**.
3. In the Tree viewer, hold **Shift** and click to multi-select the
   following three branches:
   - `All → false → F → Asian`
   - `All → false → F → Black`
   - `All → false → M → Asian`
4. **Verify:** the dataframe selection now contains 174 rows
   (`df.selection.trueCount === 174`). The grid row highlight reflects
   this subset.
5. Hold **Shift + Click** on the branch `All → true → F → Black`.
6. **Verify:** the selection extends to 176 rows
   (`df.selection.trueCount === 176`).

JS-API fallback note: `df.selection.set(idx, true)` for the matching
row indices drives the bitset, BUT does NOT exercise the canvas
hit-testing, branch-to-row-set mapping, or Shift-modifier
multi-select semantics that the Tree viewer implements via ECharts.

### Tree — Shift+Click extend selection across non-contiguous parent (canvas)

Source: extracted from `tree.md` `### Collaborative filtering` step
4 (canvas Shift+Click extension into a different parent branch
exercises the cross-parent selection-extension logic, not
synthesizable on canvas).

Setup: continue from "Tree Shift+Click branches multi-selection"
(scenario above) — selection is at 174 rows after the three initial
branch picks.

Steps:

1. In the **Filter Panel**, set **CONTROL = true** filter.
2. **Verify:** the visible row count (filter ∧ selection) is **0**.
3. Hold **Shift** and click the Tree branch `All → true → F → Black`.
4. **Verify:** the visible row count (filter ∧ selection) is **2**.
5. In the **Filter Panel**, clear the **CONTROL = true** filter.
6. **Verify:** the row count returns to **176** (the full
   selection).

JS-API fallback note: same as above — `df.selection` and `df.filter`
bitset operations replicate the cardinality contract, but do NOT
exercise the actual canvas + Filter Panel UI coordination chain.

## Notes

- **Why ui-only?** Canvas-rendered ECharts viewers (Sunburst, Tree,
  Timelines) use a single `<canvas>` element inside the ECharts
  instance container. There are no per-segment / per-branch DOM
  elements to query, hover, or click. Any "click" needs canvas
  pixel coordinates, which depend on viewport size, ECharts layout
  algorithm, and theme — making programmatic synthesis brittle and
  not equivalent to the actual UI invariant.
- **JS-API fallbacks remain in companion specs** (`sunburst-spec.ts`,
  `tree-spec.ts`) for the contract level (filter ∧ selection
  cardinalities, `setOptions`/`props.get` round-trips), but the UI-
  gesture invariant itself is documented here for manual
  verification.
- **Canonical manual-only, by design.** `sunburst-multi-selection`,
  `tree-shift-click-multi-select`, and the Select Columns
  per-column toggle are deliberately manual-only tests — not a
  deferral and not a coverage gap. There is no plan to automate
  them; this file is the authoritative catalog for human QA.

## Dataset metadata

```json
{
  "order": 35,
  "datasets": [
    "System:DemoFiles/chem/SPGI.csv",
    "System:DemoFiles/demog.csv"
  ]
}
```
