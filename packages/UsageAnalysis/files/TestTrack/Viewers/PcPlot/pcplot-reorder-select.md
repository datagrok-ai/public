---
feature: pcplot
realizes:
  - viewers.pc-plot
realizes_atlas:
  - pcplot.cp.reorder-and-select
  - pcplot-area-select-cross-viewer
  - pcplot-current-row-sync
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - pcplot-reorder-select-spec.ts
related_bugs: []
expected_results:
  - anchor: "Setup — confirm three axes AGE, HEIGHT, WEIGHT (DOM axis-slider names)"
    expectation: "After setup, the viewer renders one axis slider per configured
      column — the DOM axis-slider names read AGE, HEIGHT, WEIGHT in order"
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

## Purpose

Verifies the PC Plot's mouse interactions on the chart itself: selecting rows
by shift-dragging a rectangle between axes (additive on repeat), clearing the
selection by clicking empty space, making a row current by clicking its
polyline, and reordering the axes by dragging a column label. Selection,
current row, and axis order are read back from the table and viewer state, so
each gesture is verified by its product effect rather than by how the canvas
looks.

## Setup

1. Close all open views.
2. Open `System:DemoFiles/demog.csv`.
3. Add a PC Plot viewer (via the toolbox or Add Viewer ribbon button).
4. In the Context Panel, set the viewer's Column Names to AGE, HEIGHT, WEIGHT (three numeric axes).

## Scenarios

### Scenario 1: Polyline selection and current-row sync (baseline axis order)

Steps:
1. Confirm the PC plot displays three vertical axes (AGE, HEIGHT, WEIGHT) with polylines for each row.
2. Read the selected-row count — record the baseline value (expected: 0).
3. Read the current-row index — record the baseline value (expected: -1, no current row).
4. Shift+drag a rectangle between the AGE and HEIGHT axes to select the polylines crossing that band.
5. Verify the selected-row count rises above zero.
6. Shift+drag a second narrower band between HEIGHT and WEIGHT (adds to the selection, does not replace — in PC Plot shift+drag is additive by design, so no Ctrl modifier is needed; adding Ctrl instead replaces the selection with a smaller set).
7. Verify the selected-row count rises again (value is strictly greater than after Step 4).
8. Click empty space in the chart area (outside any polyline) to clear the selection.
9. Verify the selected-row count returns to zero (round-trip).
10. Click a visible polyline in the chart.
11. Verify the current-row index moves off -1 (a specific row is now current).

Expected:
- After shift+drag (Step 4): the selected-row count is above zero.
- After the second shift+drag (Step 6): the selected-row count is strictly larger than the value from Step 4 (additive, not replacing; Ctrl is not required).
- After click empty space (Step 8): the selected-row count is zero (round-trip verified).
- After click polyline (Step 10): the current-row index is a real row (not -1).

### Scenario 2: Axis reorder persists in columnNames, then selection still works

Steps:
1. Confirm the PC plot has axes in order AGE, HEIGHT, WEIGHT.
2. Drag the WEIGHT column label to the leftmost position (before AGE) to reorder the axes.
3. Drag the HEIGHT column label to the middle position (between WEIGHT and AGE) if needed to reach WEIGHT, HEIGHT, AGE order.
4. Read the viewer's Column Names (Context Panel) — verify the new order is reflected (WEIGHT appears before HEIGHT, HEIGHT before AGE, or the chosen drag order).
5. Verify the new axis order persists in the viewer's Column Names list (still 3 axes).
6. Shift+drag a rectangle over the reordered chart (between the first two axes).
7. Verify the selected-row count rises above zero on the reordered chart.
8. Click empty space to clear the selection.
9. Verify the selected-row count returns to zero (round-trip on the reordered chart).

Expected:
- After axis drag reorder (Step 2-3): the viewer's Column Names list reflects the new axis order (3 items, order changed).
- After shift+drag on reordered chart (Step 6): the selected-row count is above zero.
- After click empty space (Step 8): the selected-row count is zero (round-trip).

## Automation notes

- The selected-row count is read as `grok.shell.t.selection.trueCount`, the
  current-row index as `grok.shell.t.currentRowIdx`, and the axis order as
  `viewer.props.columnNames`.
- The setup softStep confirms the three configured axes by reading the rendered
  DOM axis-slider names (`axis-slider-AGE`, `axis-slider-HEIGHT`,
  `axis-slider-WEIGHT`) inside the viewer root, in order.
