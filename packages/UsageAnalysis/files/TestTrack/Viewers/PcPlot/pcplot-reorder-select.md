---
feature: pcplot
realizes:
  - pcplot.cp.reorder-and-select
realizes_atlas:
  - pcplot.cp.reorder-and-select
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - pcplot-reorder-select-spec.ts
related_bugs: []
expected_results:
  - anchor: "Scenario 1 Step 5"
    expectation: "df.selection.trueCount rises above zero after shift+drag rectangle
      over the chart"
  - anchor: "Scenario 1 Step 7"
    expectation: "df.selection.trueCount rises again (additive, not replaced) after
      a second shift+drag band — in PC Plot shift+drag is additive by design, no
      Ctrl modifier is required"
  - anchor: "Scenario 1 Step 9"
    expectation: "df.selection.trueCount returns to zero after clicking empty space
      in the chart"
  - anchor: "Scenario 1 Step 11"
    expectation: "df.currentRowIdx moves off -1 after clicking a polyline"
  - anchor: "Scenario 2 Step 5"
    expectation: "viewer columnNames list reflects the new axis order after drag-reorder"
  - anchor: "Scenario 2 Step 7"
    expectation: "df.selection.trueCount rises above zero after shift+drag selection
      on the reordered chart"
  - anchor: "Scenario 2 Step 9"
    expectation: "df.selection.trueCount returns to zero after clicking empty space
      (round-trip)"
---

# PC Plot — Axis Reorder, Polyline Selection, and Current-Row Sync

## Setup

1. Close all open views.
2. Open `System:DemoFiles/demog.csv`.
3. Add a PC Plot viewer (via the toolbox or Add Viewer ribbon button).
4. In the Context Panel, set the viewer's Column Names to AGE, HEIGHT, WEIGHT (three numeric axes).

## Scenarios

### Scenario 1: Polyline selection and current-row sync (baseline axis order)

Steps:
1. Confirm the PC plot displays three vertical axes (AGE, HEIGHT, WEIGHT) with polylines for each row.
2. In the JS API console, read `grok.shell.t.selection.trueCount` — record the baseline value (expected: 0).
3. Read `grok.shell.t.currentRowIdx` — record the baseline value (expected: -1).
4. Shift+drag a rectangle between the AGE and HEIGHT axes to select the polylines crossing that band.
5. Verify: `grok.shell.t.selection.trueCount` rises above zero.
6. Shift+drag a second narrower band between HEIGHT and WEIGHT (adds to the selection, does not replace — in PC Plot shift+drag is additive by design, so no Ctrl modifier is needed; adding Ctrl instead replaces the selection with a smaller set).
7. Verify: `grok.shell.t.selection.trueCount` rises again (value is strictly greater than after Step 4).
8. Click empty space in the chart area (outside any polyline) to clear the selection.
9. Verify: `grok.shell.t.selection.trueCount` returns to zero (round-trip).
10. Click a visible polyline in the chart.
11. Verify: `grok.shell.t.currentRowIdx` moves off -1 (a specific row is now current).

Expected:
- After shift+drag (Step 4): `df.selection.trueCount > 0`.
- After the second shift+drag (Step 6): `df.selection.trueCount` is strictly larger than the value from Step 4 (additive, not replacing; Ctrl is not required).
- After click empty space (Step 8): `df.selection.trueCount == 0` (round-trip verified).
- After click polyline (Step 10): `df.currentRowIdx >= 0` (not -1).

### Scenario 2: Axis reorder persists in columnNames, then selection still works

Steps:
1. Confirm the PC plot has axes in order AGE, HEIGHT, WEIGHT.
2. Drag the WEIGHT column label to the leftmost position (before AGE) to reorder the axes.
3. Drag the HEIGHT column label to the middle position (between WEIGHT and AGE) if needed to reach WEIGHT, HEIGHT, AGE order.
4. Read `viewer.props.columnNames` (or the Context Panel's Column Names field) — verify the new order is reflected (WEIGHT appears before HEIGHT, HEIGHT before AGE, or the chosen drag order).
5. Verify: The new axis order persists in the viewer's columnNames list (product-state count: still 3 axes).
6. Shift+drag a rectangle over the reordered chart (between the first two axes).
7. Verify: `grok.shell.t.selection.trueCount` rises above zero on the reordered chart.
8. Click empty space to clear the selection.
9. Verify: `grok.shell.t.selection.trueCount` returns to zero (round-trip on the reordered chart).

Expected:
- After axis drag reorder (Step 2-3): `viewer.props.columnNames` reflects the new axis order (3 items, order changed).
- After shift+drag on reordered chart (Step 6): `df.selection.trueCount > 0`.
- After click empty space (Step 8): `df.selection.trueCount == 0` (round-trip).
