---
feature: boxplot
realizes_atlas:
  - boxplot.cp.pointer-select-highlight
realizes:
  - viewers.box-plot
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: github-2764
    status: fixed
  - id: github-3066
    status: fixed
  - id: github-3065
    status: fixed
  - id: GROK-19950
    status: fixed
  - id: GROK-20502
    status: fixed
expected_results:
  - anchor: Step 3
    expectation: >-
      df.currentRowIdx equals the rowId delivered by the d4-boxplot-point-click
      event for the clicked marker — both signals agree on the same row.
  - anchor: Step 4
    expectation: >-
      df.selection.trueCount is greater than zero AND countSelectionHuePixels
      against selectedRowsColor grows above its no-selection baseline (the
      default coloring auto-follows Category 1, so its categorical palette
      paints a small non-zero baseline inside the selection-hue range; the
      selection band adds clearly above it).
  - anchor: Step 5
    expectation: >-
      After applying a non-default categorical coloring (marker color column
      switched away from the auto-selected Category 1 column),
      df.selection.trueCount is greater than zero after the second shift-drag
      AND countSelectionHuePixels grows above that coloring's own no-selection
      baseline (github-3066 — selection highlight survives non-default
      coloring).
  - anchor: Step 6
    expectation: >-
      After ctrl-clicking a category label whose rows are DISJOINT from the
      current selection, df.selection.trueCount grows by exactly the count of
      rows in that category (additive union, not a replace).
  - anchor: Step 7
    expectation: >-
      After a plain click on a single category label, df.selection.trueCount
      equals the total row count for that category.
  - anchor: Step 8
    expectation: >-
      After clicking empty space, df.selection.trueCount equals zero (selection
      cleared without resetting the view).
  - anchor: Step 9
    expectation: >-
      After a single click on empty space, no d4-boxplot-reset-view event is
      fired (props.viewport is not populated by pointer interactions on the live
      viewer, so the reset-view event is the observable viewport signal — only a
      double-click calls resetView()).
  - anchor: Step 10
    expectation: >-
      After a double-click on empty space, the d4-boxplot-reset-view event fires
      (resetView is called), distinct from the single click on empty space which
      fires none.
  - anchor: Step 12
    expectation: >-
      With the left-most category filtered out, iterating selection.get(i) over
      all filtered-out row indices returns false for every index — zero
      selection leak into hidden rows (github-2764 coordinate-mapping
      regression).
  - anchor: Step 14
    expectation: >-
      With showSelectedRows=false and a non-empty selection,
      countSelectionHuePixels against selectedRowsColor collapses from the
      selected level back to the categorical-coloring baseline — a non-zero
      remainder, because the categorical palette itself paints pixels inside the
      selection-hue range; that remainder plus the still-set marker color column
      is the evidence that regular color coding is preserved (GROK-19950). The
      coloring legend on this viewer is canvas-drawn — there are no legend DOM
      elements to assert.
  - anchor: Step 15
    expectation: >-
      With rowSource=Selected and showSelectedRows=true and a non-empty
      selection, countSelectionHuePixels against selectedRowsColor returns zero
      — the showSelectedRows dependsOn gate suppresses the selection-color
      override for rowSource Selected/FilteredSelected (GROK-19950 follow-up).
  - anchor: Step 17
    expectation: >-
      A .d4-tooltip element is present in the DOM while hovering a marker point
      (tooltip visible).
  - anchor: Step 18
    expectation: >-
      After setting showMouseOverPoint=false, no hover-highlight overlay is
      drawn over the canvas — the canvas diff between a before-hover and
      after-hover settle-gate snapshot is approximately zero (no-error floor).
  - anchor: Step 20
    expectation: >-
      After hovering a bar on a neighboring Bar chart viewer with distinct
      categorical coloring, the box plot canvas pixel counts are unchanged
      compared to the pre-hover reference snapshot — no cross-viewer color
      repaint (github-3065).
  - anchor: Step 22
    expectation: >-
      After selecting one category's rows and deleting them from the grid, the
      categorical-coloring box plot still renders without errors: the data
      canvas keeps painting, no new shell warnings appear, and the marker color
      column remains applied (GROK-20502 — the coloring legend is canvas-drawn
      on this viewer; a vanished legend would surface as a coloring drop or a
      render fault).
realized_as:
  - boxplot-pointer-select-spec.ts
---

# Box Plot — Pointer Selection and Highlight

## Setup

1. Close all open tables and viewers.
2. Open the demog dataset: navigate to Files > App Data and open
   System:DemoFiles/demog.csv.
3. Add a Box Plot viewer to the table view (toolbar or viewer menu "Add viewer > Box Plot").
4. In the Context Panel > Value, set Value to AGE.
5. In the Context Panel > Data, set Category 1 to RACE.
6. In the Context Panel > Markers, set Marker Size to 10 (ensures individual point markers
   are large enough for click targeting).

## Scenarios

### Scenario 1: Marker click sets the current row and fires a point-click event

Steps:
1. With Value=AGE and Category 1=RACE and Marker Size=10, locate a visible individual
   marker point on the box plot canvas.
2. Register a one-shot listener for the point-click event on the box plot
   viewer before clicking.
3. Click the marker point with a single left-click.
4. Verify that the current row index is set to a valid row (not empty) AND that the
   point-click event was fired with a row identifier matching the current row —
   both signals agree on the same row.

Expected:
- The current row index matches the row identifier delivered by the point-click event
  for the clicked marker — both signals agree on the same row.

### Scenario 2: Shift-drag creates a selection with visible hue

Steps:
1. Hold Shift and drag across a band of marker points within one category column on the canvas.
2. Release the drag.
3. Verify that rows are selected (selection count is greater than zero) and that
   the selection highlight color is visible on the canvas above the coloring's own
   baseline (the marker color auto-follows Category 1, so the default coloring is the
   RACE categorical palette, which itself paints a small amount of selection-like hue).
4. In the Context Panel > Color, set Marker Color Column to SEX (a non-default
   categorical coloring — a deliberate switch away from the auto-selected Category 1
   column).
5. Hold Shift and drag again across a band of marker points in the same category area.
6. Verify that rows are selected after the second shift-drag AND that the selection
   highlight color is still visible on the canvas above that coloring's own baseline —
   selection highlight must survive non-default coloring (regression guard for
   github-3066).
7. In the Context Panel > Color, set Marker Color Column to None (restore default
   coloring).

Expected:
- After the first shift-drag, rows are selected and the selection highlight color is
  visible on the canvas above the coloring's no-selection baseline.
- After applying the non-default categorical coloring and shift-dragging again, rows
  are selected and the selection highlight color is still visible above that coloring's
  own baseline (github-3066 — selection highlight survives non-default coloring).

### Scenario 3: Category label click — plain, ctrl-click, and clear

Steps:
1. Click empty space on the box plot canvas to clear any current selection; verify
   no rows are selected before proceeding.
2. Use shift-drag to build a partial selection inside one category column (fewer rows than
   the full category count). Record the current selection count as the baseline.
3. Identify a second category column whose rows are DISJOINT from the currently selected
   rows (no overlap). Ctrl-click the category label for that second category.
4. Verify that the selection count grew additively by exactly the full row count of
   the ctrl-clicked category — the result equals the baseline count plus that category's row count.
5. Click the category label of a single category (plain left-click, no modifier keys).
6. Verify that the selection count equals the total row count for that category
   (plain category-label click replaces the selection with the full category).
7. Click empty space on the box plot canvas (single click, no modifiers).
8. Verify that the selection count equals zero (selection cleared).
9. Verify that the view was not reset by the single click — no reset-view event fires
   (the viewer does not populate its viewport property from pointer interactions, so
   the reset-view event is the observable signal; only a double-click resets the view).
10. Double-click on empty space on the box plot canvas.
11. Verify that the double-click reset the view — the reset-view event fires, distinct
    from the no-event result of a single click.

Expected:
- After ctrl-clicking a category label whose rows are DISJOINT from the current
  selection, the selection count grows by exactly the count of rows in that category
  (additive union, not a replace).
- After a plain click on a single category label, the selection count equals the
  total row count for that category.
- After clicking empty space, the selection count equals zero (selection cleared
  without resetting the view).
- After a single click on empty space, no reset-view event fires (a single click does
  not reset the view — only a double-click resets it).
- After a double-click on empty space, the reset-view event fires (the view resets),
  distinct from the no-event result of a single click on empty space.

### Scenario 4: No selection leak into filtered-out categories (github-2764)

Steps:
1. In the Filter Panel, open the Category 1 / RACE filter and uncheck the left-most
   category (the category that appears first / leftmost on the box plot X axis) so it is filtered out.
2. Hold Shift and drag across the visible area of the box plot canvas, covering the area
   where the filtered-out category previously appeared.
3. For every row that belongs to the filtered-out category, confirm it is not selected.
4. Verify that the selection only reflects rows in visible (non-filtered-out)
   categories — zero leak into the hidden rows.
5. Reset the RACE filter so all categories are shown again, and close the Filter Panel.

Expected:
- With the left-most category filtered out, none of the filtered-out rows are selected
  after a shift-drag over their former canvas area — zero selection leak into hidden rows
  (github-2764 coordinate-mapping regression).

### Scenario 5: Show Selected Rows and Row Source gate (GROK-19950)

Steps:
1. In the Context Panel > Color, set Marker Color Column to RACE (categorical coloring
   active, legend visible).
2. Click the category label for one of the RACE categories to select its rows.
3. Confirm rows are selected and the selection highlight color is visible on the canvas
   (baseline — selection hue visible).
4. In the Context Panel > Selection, set Show Selected Rows to false.
5. Verify that the selection highlight color drops back to the coloring's baseline level
   on the canvas AND that the RACE categorical coloring itself remains applied (regular
   color coding is preserved, only the selection color override is off; the coloring
   legend on this viewer is drawn on the canvas — there is no legend DOM to inspect).
6. Set Show Selected Rows back to true.
7. In the Context Panel > Data, set Row Source to Selected.
8. Verify that the selection highlight color is no longer visible on the canvas — when
   Row Source is set to Selected, the Show Selected Rows setting suppresses the selection
   color override (GROK-19950 follow-up).
9. In the Context Panel > Data, set Row Source back to All.
10. In the Context Panel > Color, set Marker Color Column to None.

Expected:
- With Show Selected Rows set to false and a non-empty selection, the selection highlight
  color collapses back to the categorical coloring's own baseline level while the
  coloring itself remains applied (GROK-19950 — regular color coding is preserved; the
  coloring legend is canvas-drawn on this viewer).
- With Row Source set to Selected and Show Selected Rows set to true and a non-empty
  selection, the selection highlight color is not visible on the canvas — the
  Show Selected Rows setting suppresses the selection color override for the Selected and
  FilteredSelected row source modes (GROK-19950 follow-up).

### Scenario 6: Hover tooltip and Show Mouse Over Point

Steps:
1. Hover the mouse cursor over an individual marker point on the box plot canvas and
   hold for the tooltip settle time.
2. Verify that a tooltip is visible (a tooltip element appears in the page).
3. In the Context Panel > Misc (or Context Panel > Tooltip), set Show Mouse Over Point to false.
4. Hover the mouse cursor over the same area of the canvas and hold.
5. Capture a canvas snapshot before hover and after the hover settle, and verify the
   canvas is unchanged — no hover-highlight overlay is drawn (overlay suppressed).

Expected:
- A tooltip is visible while hovering a marker point.
- After setting Show Mouse Over Point to false, no hover-highlight overlay is drawn over
  the canvas — the canvas before and after the hover are the same (no-error floor).

### Scenario 7: Cross-viewer highlight must not recolor box plot (github-3065)

Steps:
1. Add a Bar Chart viewer to the same table view alongside the Box Plot.
2. In the Box Plot's Context Panel > Color, set Marker Color Column to RACE (categorical coloring).
3. Capture a reference canvas snapshot of the box plot.
4. Hover a bar inside the Bar Chart viewer (hover over a different category's bar to
   trigger cross-viewer highlighting).
5. After the highlight settle, capture a second canvas snapshot of the box plot.
6. Verify the box plot canvas is unchanged compared to the reference snapshot — the box
   plot canvas color did not repaint due to the cross-viewer hover (github-3065).
7. Remove the Bar Chart viewer (close or remove from layout).
8. In the Box Plot's Context Panel > Color, set Marker Color Column to None.

Expected:
- After hovering a bar on a neighboring Bar chart viewer with distinct categorical
  coloring, the box plot canvas is unchanged compared to the pre-hover reference
  snapshot — no cross-viewer color repaint (github-3065).

### Scenario 8: Categorical coloring survives row deletion (GROK-20502)

Steps:
1. In the Context Panel > Color, set Marker Color Column to RACE (categorical coloring;
   the coloring legend is drawn on the viewer canvas).
2. Click the RACE category label for one specific category to select its rows.
3. Verify the categorical coloring is applied (colored markers and the canvas-drawn
   legend are visible on the viewer).
4. Open the grid viewer and with the selected rows still active, delete the selected
   rows from the grid (Edit > Delete Rows, or the grid context menu delete action).
5. Verify that the box plot still renders its categorical coloring without errors —
   the canvas keeps painting, no warnings appear, and the marker color column remains
   applied; the legend must not disappear after the data mutation (GROK-20502).
6. Restore the deleted rows by pressing Ctrl+Z (undo) so the dataset is intact for any
   downstream teardown.

Expected:
- After selecting one category's rows and deleting them from the grid, the box plot
  still renders its categorical coloring without errors — the canvas keeps painting,
  no new warnings appear, and the marker color column remains applied (GROK-20502 —
  the coloring legend is canvas-drawn on this viewer; a vanished legend would surface
  as a coloring drop or a render fault).

## Automation notes

- target_layer rationale: all selection, hover, and cross-viewer highlight behaviors
  require a live browser with canvas rendering and DOM inspection; Playwright is the
  only layer that can perform shift-drag gestures and read selection state in
  a real Datagrok session.
- countSelectionHuePixels (TestTrack/helpers/viewers.ts) is used POINT-WISE against
  selectedRowsColor — never as a global non-white count. Under any categorical
  coloring the palette itself paints pixels inside the selection-hue range, so the
  honest reading is a DELTA over a measured no-selection baseline (and "suppressed"
  means "returns to that baseline", not "drops to zero"); the absolute == 0 reading
  is valid only where the coloring itself is not rendered (rowSource=Selected gate).
- Setting Category 1 auto-sets the marker color column to the same column, so the spec's
  "default coloring" is the Category 1 categorical palette; a NON-default categorical
  coloring is any deliberate switch of the marker color column away from it.
- The box plot's coloring legend is CANVAS-drawn: no [name="legend"] host and no
  .d4-legend-item elements exist under marker coloring. Legend expectations must use
  canvas/prop channels (render floor + marker color column read), never legend DOM.
- Category-label click band: the clickable category-selection strip sits just below the
  plot area (a narrow band around ~0.88-0.90 of the canvas height); the ink rows
  further below are the statistics table and are click-inert. Locate the band
  empirically — plain-click candidate Y fractions on the first category's center until
  the selection equals that category's row count — and reuse the found Y.
- The ctrl-click category-label precondition requires rows disjoint from the current
  selection: ctrl-click on a fully-selected category DESELECTS it
  (viewer_utils.dart:14-15). Confining the baseline shift-drag to the first category's
  X band guarantees disjointness by construction.
- The double-click reset is implemented at box_plot_core.dart:969-972. props.viewport
  is not populated by pointer interactions, so the observable reset signal is the
  d4-boxplot-reset-view event (count of events fired), not a viewport prop diff.
- Deleting the selected rows: the selected indices are non-contiguous, so use
  df.rows.removeWhereIdx over the captured df.selection.getSelectedIndexes() set;
  rows.removeAt(idx, count) removes a CONTIGUOUS block and deletes the wrong rows.
- Scenario 8 is the final scenario: the spec closes all tables at the end and the
  dataset is re-opened from the CSV on every run, so the automated run needs no undo;
  the Ctrl+Z step applies to manual execution.
- Cross-viewer canvas snapshot diff: use settle-gate (wait for idle/canvas-settle event)
  before capturing both snapshots.
