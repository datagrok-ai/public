---
feature: charts
realizes_atlas: [charts.cp.click-segment-to-select-or-filter]
realizes: [charts.viewer.sunburst, charts.viewer.tree]
priority: p0
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  The Sunburst and Tree viewers render to a single canvas with no
  per-segment / per-branch DOM elements, so these clicks and modifier
  keys must be driven by hand.
related_bugs: []
---

# Charts — UI-only manual scenarios (canvas-gesture)

These scenarios require a human operator to drive the actual canvas
click events — the gestures below cannot be synthesized by automation.

## Setup
Each scenario block opens its own dataset.

## Scenarios

### Sunburst — multi-selection on segments (canvas, Click / Ctrl+Click / Ctrl+Shift+Click)

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

### Sunburst — empty (null) category click (canvas)

1. Open SPGI.csv and add a Sunburst viewer.
2. Configure a column known to contain nulls (e.g. **Sampling Time**)
   via Select Columns.
3. **Click** the null (grey) segment.
4. **Verify:** the segment behaves like any other category — the
   corresponding rows (those with null in that column) are selected
   or filtered (per the configured `onClick` action — Select or
   Filter).

### Tree — Shift+Click branches multi-selection (canvas)

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

### Tree — Shift+Click extend selection across non-contiguous parent (canvas)

Setup: continue from "Tree Shift+Click branches multi-selection"
(scenario above) — selection is at 174 rows after the three initial
branch picks.


1. In the **Filter Panel**, set **CONTROL = true** filter.
2. **Verify:** the visible row count (filter ∧ selection) is **0**.
3. Hold **Shift** and click the Tree branch `All → true → F → Black`.
4. **Verify:** the visible row count (filter ∧ selection) is **2**.
5. In the **Filter Panel**, clear the **CONTROL = true** filter.
6. **Verify:** the row count returns to **176** (the full
   selection).


---
{
  "order": 35,
  "datasets": ["System:DemoFiles/chem/SPGI.csv", "System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
