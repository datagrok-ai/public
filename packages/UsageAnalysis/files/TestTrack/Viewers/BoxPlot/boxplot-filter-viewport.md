---
feature: boxplot
realizes_atlas:
  - boxplot.cp.filter-viewport
realizes:
  - viewers.box-plot
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20171
    status: fixed
realized_as:
  - boxplot-filter-viewport-spec.ts
expected_results:
  - anchor: Scenario 1 Step 4
    expectation: >-
      After narrowing the df.filter on Average Mass, the Box Plot viewport
      narrows (the visible value range on the Y axis contracts to match the
      filtered data) and df.filter.trueCount is less than 100 (the full spgi-100
      row count). The zoomValuesByFilter property reads true.
  - anchor: Scenario 2 Step 6
    expectation: >-
      With zoomValuesByFilter set to false and the same Average Mass filter
      applied, the Box Plot viewport is unchanged from the Step 3 baseline
      capture — the range slider bounds do not move. The filter itself is active
      (df.filter.trueCount < 100), only the viewport response is suppressed.
  - anchor: Scenario 2 Step 8
    expectation: >-
      After switching zoomValuesByFilter back to true, the Box Plot viewport
      narrows to match the current filtered data range. The transition pair is
      confirmed: false suppresses the viewport follow, true restores it.
  - anchor: Scenario 3 Step 5
    expectation: >-
      After applying the viewer-local formula filter set to '${Average Mass} >
      350' on top of the active df.filter, the df.filter.trueCount is UNCHANGED
      from the pre-formula baseline — the formula filter is viewer-local and
      does not affect the DataFrame's global filter. The rendered marker
      population in the viewer canvas drops (canvas floor: settle-gated diff
      shows fewer markers).
  - anchor: Scenario 4 Step 5
    expectation: >-
      With showEmptyCategories off, the box plot drops the empty-valued Series
      category from the X axis: the plot re-lays out around the remaining
      categories (large settle-gated repaint of the data canvas, absent when no
      empty-valued category exists) and the category-label band at the bottom of
      the canvas redraws with the reduced label set.
  - anchor: Scenario 4 Step 7
    expectation: >-
      After re-enabling showEmptyCategories, the empty-valued Series category
      reappears on the X axis: the plot re-lays out to re-add its axis slot
      (settle-gated data-canvas repaint) and the category-label band redraws
      with the restored label set.
  - anchor: Scenario 5 Step 4
    expectation: >-
      After setting Marker Color Column to TPSA with an active filter, the box
      plot repaints, the color-scale footprint appears on the overlay canvas
      (reserved, non-empty scale area — GROK-20171 regression floor), and no new
      console/page error is emitted.
  - anchor: Scenario 5 Step 5
    expectation: >-
      The color scale reflects the filtered (visible) TPSA data range rather
      than the full-dataset range (getVisibleColorRange semantics): with Marker
      Color set to TPSA, widening the filter changes the visible TPSA range and
      the scale drawn on the overlay canvas re-renders its endpoints
      (settle-gated overlay repaint; the overlay is untouched by filter changes
      when no color column is set).
---

# Box Plot — Filter Semantics and Viewport Response

## Setup

1. Close all open tables and viewers.
2. Open the spgi-100 dataset: navigate to Files > App Data and open
   System:AppData/Chem/tests/spgi-100.csv. Confirm the opened table
   contains 100 rows (the in-workspace frame name may differ from the
   file name).
3. Add a Box Plot viewer to the table view via the toolbar (Add viewer
   > Box Plot).
4. In Context Panel > Data, set Value to Average Mass and set Category 1
   to Series.

## Scenarios

### Scenario 1: Zoom Values By Filter on — viewport follows the filter

Steps:
1. Capture the current Box Plot viewport (note the visible Y-axis range
   shown by the range slider handles, before any filter is applied).
2. Verify that all 100 rows are currently passing the filter (no active
   filter is in place — the grid status bar shows 100 of 100 rows).
3. Verify that Zoom Values By Filter in Context Panel > Value
   is enabled (it should be on by default).
4. Open the Filter panel (toolbar Filter button or View > Filter Panel).
   In the Average Mass filter, narrow the range to exclude the lowest and
   highest values (drag the range handles inward so fewer rows pass).
5. Verify that the Box Plot viewport narrowed (the visible Y-axis range
   is smaller than the baseline captured in Step 1) AND that fewer than
   100 rows are now passing (the grid status bar shows a reduced row
   count).

Expected:
- After narrowing the Average Mass filter, the Box Plot viewport
  narrows (visible Y-axis range contracts to the filtered data) and
  the grid status bar confirms fewer than 100 rows pass.

### Scenario 2: Zoom Values By Filter off — viewport stays fixed

Steps:
1. In Context Panel > Value, disable Zoom Values By Filter (set it to false).
2. Reset the Average Mass filter to include all rows (clear or widen it
   fully so the grid status bar returns to 100 of 100 rows).
3. Capture the Box Plot viewport baseline again (note the Y-axis range with
   all rows visible).
4. In the Filter panel, apply the same Average Mass range restriction as in
   Scenario 1 (narrow it so fewer rows pass).
5. Verify that the filter is active — the grid status bar shows fewer than
   100 rows passing.
6. Verify that the Box Plot viewport did NOT change — the range slider
   bounds and the Y-axis range are unchanged from the Step 3 baseline
   (disabling Zoom Values By Filter suppresses the viewport follow).
7. In Context Panel > Value, re-enable Zoom Values By Filter (set it to true).
8. Verify that the Box Plot viewport now narrows to match the current
   filtered data range (switching back to true restores the
   viewport-follow behavior).

Expected:
- With Zoom Values By Filter disabled and an active filter, the Box Plot
  viewport is unchanged from the pre-filter baseline.
- After re-enabling Zoom Values By Filter, the viewport narrows to the
  filtered data bounds.

### Scenario 3: Viewer-local formula filter does not affect the table filter

Steps:
1. Ensure an active Average Mass filter is in place from the previous
   scenario (the grid status bar shows fewer than 100 rows passing).
   Note the current passing row count shown in the status bar.
2. In Context Panel > Value > Filter, type the formula
   '${Average Mass} > 350' and confirm it (the viewer-local formula filter).
3. Observe the Box Plot canvas — note that fewer markers are rendered
   (the visible marker population in the viewer decreases).
4. Check the grid status bar for the current passing row count.
5. Verify that the passing row count shown in the grid status bar is
   UNCHANGED from the value noted in Step 1 — the formula filter is
   viewer-local and must not affect the table's global filter row count.

Expected:
- After applying the viewer-local formula filter '${Average Mass} > 350',
  the passing row count in the grid status bar is unchanged from the
  pre-formula baseline — the formula filter is viewer-local and does not
  modify the global table filter. The rendered marker population in the
  viewer canvas drops.

### Scenario 4: Show Empty Categories — axis structural count

Steps:
1. Clear the Average Mass filter and the formula filter so all rows
   are visible (the grid status bar shows 100 of 100 rows). Reset the
   formula filter in Context Panel > Value > Filter to empty.
2. Inspect the spgi-100 dataset for a Series category value that has no
   Average Mass values (all null). If no such category exists, add a
   calculated column via the grid's Add New Column (plus icon) with a
   formula that produces null for one category of Series, naming the
   column 'AverageMassFixture'; set Value to AverageMassFixture and
   Category 1 to Series. (Record in a note that the fixture column was
   added and must be removed in the finally step.)
3. In Context Panel > Data, verify that Show Empty Categories is enabled
   (on by default). Count the number of category labels visible on the
   Box Plot X axis — note this count.
4. In Context Panel > Data, disable Show Empty Categories (turn it off).
5. Verify that the X-axis category count decreased (the empty-valued
   category disappeared from the axis). The category whose VALUE cells
   are all null is no longer shown on the X axis.
6. In Context Panel > Data, re-enable Show Empty Categories (turn it on).
7. Verify that the X-axis category count returns to the count noted in
   Step 3 (the empty-valued category reappears).
8. If a fixture column was added in Step 2, remove it now: right-click
   the column header in the grid > Remove Column (and in a finally
   context, ensure this cleanup runs even if assertions above failed).

Expected:
- With Show Empty Categories off, the count of X-axis category labels
  decreases by the number of categories with all-null VALUE cells (those
  categories disappear from the axis). With it on, the empty-valued
  category is present on the axis.

### Scenario 5: Bin color recomputes to visible range; legend footprint (GROK-20171)

Steps:
1. Ensure an active Average Mass filter is in place (narrow the Filter
   panel Average Mass range so fewer than 100 rows are passing — confirm
   in the grid status bar).
2. In Context Panel > Color, set Marker Color Column to TPSA (a numerical
   column present in spgi-100).
3. Observe the Box Plot — the color scale is applied to TPSA values.
4. Verify that the color scale legend is present in the viewer and does
   not overlap neighboring chrome elements or produce a console error
   (GROK-20171 regression guard: the color-scale legend keeps a reserved
   footprint; confirm no overlap between the legend area and the adjacent
   axis chrome, and no new console errors since the color step).
5. Verify that the color scale range reflects the filtered (visible) data
   range of TPSA, not the full-dataset range — note the endpoint values
   shown on the color scale, then widen the Average Mass filter back to
   all rows and confirm the scale endpoints change to cover the wider
   visible range (the scale tracks the visible data).

Expected:
- After setting Marker Color to TPSA with an active filter, the color-scale
  footprint is present with no overlap against neighboring chrome
  and no console errors (GROK-20171 regression guard). The color scale
  endpoints reflect the filtered TPSA range, not the full-dataset range,
  and change when the filter range changes.

## Automation notes

- target_layer rationale: the filter-viewport interaction requires driving
  the Filter panel range slider (UI actuation), reading the Box Plot
  viewport bounds (range slider handle positions or bp.viewport prop),
  and checking df.filter.trueCount; the showEmptyCategories structural
  axis count and the formula filter input in the Context Panel also require
  a live browser session. Playwright is the appropriate layer.
- zoomValuesByFilter transition pair (Scenario 2 Steps 7-8): only the
  switch itself produces a signal — the viewport before-vs-after the
  switch is the assertion. Do NOT assert the viewport before the switch
  (that was Scenario 1's concern); assert only the delta on the switch.
- Formula filter channel (Scenario 3): the formula is props.filter. It
  feeds createCombinedRowFilter(look.rowSource, look.filter)
  (box_plot_core.dart:607) and is viewer-internal — it NEVER touches
  df.filter. The assertion is df.filter.trueCount == the pre-formula
  value. The canvas change (fewer markers) is a settle-gated diff floor
  (not a pixel read).
- showEmptyCategories fixture (Scenario 4): spgi-100 may or may not have
  a null-valued category for Average Mass. The fixture addNewCalculated
  path is the deterministic fallback (a calculated copy of Average Mass
  with one category's values nulled). The fixture column must be removed
  in a finally block regardless of pass/fail.
- showEmptyCategories signal (Scenario 4 Steps 5/7): the axis is
  canvas-drawn with no public axis-category-count accessor. The honest
  structural channel is a settle-gated color diff of the data canvas
  (adding/removing an axis slot re-spaces every box — a large re-layout
  repaint; the same toggle with no empty-valued category repaints
  nothing, which makes the diff attributable) plus a second diff scoped
  to the bottom category-label band (the label set itself changes). Guard
  the precondition by computing per-category null counts from the
  dataframe before toggling.
- GROK-20171 floor (Scenario 5 Step 4): the color scale has no DOM — it
  is drawn on the overlay canvas, which stays empty until a color column
  is set. The footprint assertion is the appearance of non-white ink on
  the overlay after setting the color column, plus a zero console/page
  error delta. Do NOT attempt to read the scale's typographic values from
  canvas pixels.
- Visible range color recompute (Scenario 5 Step 5): the Color Min /
  Color Max properties hold only manual bounds and stay empty — the
  computed visible range has no prop or DOM read. The channel is a
  settle-gated color diff of the overlay canvas across a filter change:
  the overlay carries only the color scale (park the pointer so no hover
  ink appears), and a filter change with no color column leaves it
  untouched, so an endpoint re-render is isolated. Guard the trigger by
  computing the visible TPSA min/max from the dataframe under both filter
  states and asserting they differ. This is the getVisibleColorRange
  semantic per column_coloring_mixin.dart.
