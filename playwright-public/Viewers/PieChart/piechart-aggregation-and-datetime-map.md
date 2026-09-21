---
feature: piechart
realizes_atlas:
  - piechart.cp.aggregation-tour-and-datetime-map
realizes:
  - viewers.pie-chart
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-9645
    status: fixed
  - id: GROK-18482
    status: fixed
realized_as:
  - piechart-aggregation-and-datetime-map-spec.ts
expected_results:
  - anchor: "Scenario 1 — angle aggregation tour"
    expectation: "Cycling the angle aggregation through min, max, med, stdev,
      count and avg reads back each value in order, and the whole tour raises
      no page or console error with the Pie chart root element still in the DOM"
  - anchor: "Scenario 1 — angle aggregation tour"
    expectation: "Switching the angle aggregation from avg to sum on AGE
      redistributes the slice angles: the settle-gated per-color canvas
      histogram delta between the two states exceeds 500 pixels"
  - anchor: "Scenario 1 — angle aggregation tour"
    expectation: "Assigning WEIGHT as the Segment Length Column shrinks every
      non-maximal category's outer radius, so the settle-gated canvas ink count
      drops by more than 500 pixels from the full disc, and clearing both the
      angle and length columns restores the full disc (the ink count rises back
      by more than 500 pixels) with both column properties reading back empty"
  - anchor: "Scenario 2 — negative-minimum and all-zero aggregations"
    expectation: "With a scratch column whose min aggregation is negative
      assigned as the angle column, the viewer renders a validation message
      containing 'contains negative values' and the column name; with an
      all-zero scratch column under sum, the message contains 'all values
      are 0'; clearing the angle column removes the message"
  - anchor: "Scenario 3 — Category Map on STARTED"
    expectation: "With Category set to the STARTED date column the Category Map
      defaults to year; switching it to month and then quarter changes the
      legend category set each time (the three label sets are pairwise
      different, month has at most 12 entries, quarter at most 4), and the
      category is reverted to RACE afterwards"
  - anchor: "Scenario 4 — grid color-coding on RACE"
    expectation: "Applying a categorical color to a RACE category through the
      grid column's color coding turns that category's pie legend entry to the
      requested color, and clearing the color coding restores the original
      legend color (round-trip)"
  - anchor: "Scenario 5 — tooltip content follows the configured aggregations"
    expectation: "Hovering a slice with the default count aggregation shows a
      tooltip whose text names a RACE category and includes a row-count digit;
      after setting the angle column to AGE with avg the tooltip text contains
      'avg(AGE)'; after adding WEIGHT as the length column with max it
      additionally contains 'max(WEIGHT)'; after the mouse moves to an empty
      corner and leaves the chart, no new text is rendered into the tooltip
      and it is not visibly shown"
---

# Pie Chart — Aggregation Tour, Validation Messages, and DateTime Category Map

## Purpose

The pie chart always paints a full disc, so most aggregation switches recolor
the canvas without changing how many pixels are inked. What can be checked
honestly is threefold: every aggregation in the tour completes without an
error and reads back its value; the switches that demonstrably change
geometry are confirmed by measuring the canvas — the slices recolor when the
angles redistribute, and the inked area shrinks under length coding; and the
documented validation messages for impossible aggregations appear as readable
text on the viewer. The date Category Map and
the grid color-coding scenarios read the legend, which is real page content:
the category set and the entry colors are asserted directly.

## Setup

1. Close all open views.
2. Open the demog dataset `System:DemoFiles/demog.csv` and wait for the table view.
3. Add a Pie chart viewer from the Toolbox.
4. Set **Category** to RACE and wait for the render.
5. Note the current page-error and console-error counts so each scenario can
   confirm that no new errors appeared.

## Scenarios

### Scenario 1: Angle aggregation tour and length column switching

Steps:
1. Note the current page-error and console-error counts.
2. Set **Segment Angle Column** to AGE.
3. Cycle **Segment Angle Aggr Type** through min, max, med, stdev, count and
   avg, waiting for the render after each switch, and verify each value reads
   back in order.
4. With the aggregation at avg, record the colors drawn on the canvas. Switch
   the aggregation to sum and record them again. Check that the two states
   clearly differ (exact margins in Automation notes) — avg(AGE) is nearly
   uniform across races while sum(AGE) follows the row counts, so the slice
   angles must redistribute even though the disc stays full.
5. Set **Segment Length Column** to WEIGHT with **Segment Length Aggr Type**
   avg, wait for the render, and measure the ink drawn on the canvas. Check it
   is clearly below the full-disc measurement from step 4 — length coding
   shrinks every non-maximal category's outer radius, so ink must decrease.
6. Switch **Segment Length Aggr Type** to min and then max, waiting after each
   — these changes are covered by the scenario's no-new-errors check.
7. Restore both aggregation types to their captured defaults, then clear
   **Segment Angle Column** and **Segment Length Column** (in that order —
   clearing a column while a non-default aggregation is still set blanks the
   chart with a validation message). Check the canvas ink rises back clearly
   above the length-coded measurement (the full disc is restored), both
   column settings read back empty, and no validation message is shown.
8. Check that no new page or console error appeared and that the Pie chart
   viewer element is still present on the page.

Expected:
- Cycling the angle aggregation through min, max, med, stdev, count and avg reads back each value in order, and the whole tour raises no page or console error with the Pie chart root element still in the DOM.
- Switching the angle aggregation from avg to sum on AGE redistributes the slice angles: the settle-gated per-color canvas histogram delta between the two states exceeds 500 pixels.
- Assigning WEIGHT as the Segment Length Column drops the settle-gated canvas ink count by more than 500 pixels from the full disc, and clearing both the angle and length columns restores the full disc — the ink count rises back by more than 500 pixels and both column properties read back empty.

### Scenario 2: Validation messages for negative and all-zero aggregations

The chart refuses to draw when the angle aggregation produces a negative
minimum or an all-zero result, and says so in a message rendered on the
viewer. Demog has no column that triggers either case, so the scenario adds
two scratch columns and removes them afterwards.

Steps:
1. Note the current page-error and console-error counts.
2. Add a scratch column NEG_PROBE holding alternating -5 and 3, and a scratch
   column ZERO_PROBE holding only zeros.
3. Set **Segment Angle Column** to NEG_PROBE with **Segment Angle Aggr Type**
   min. Verify the viewer shows a message containing
   "contains negative values" and the column name NEG_PROBE.
4. Set **Segment Angle Column** to ZERO_PROBE with aggregation sum. Verify the
   viewer shows a message containing "all values are 0" and the column name
   ZERO_PROBE.
5. Restore the default angle aggregation, then clear **Segment Angle Column**.
   Verify the message disappears and the chart renders again.
6. Remove both scratch columns (also on failure), restoring the default
   aggregation and clearing the angle column first so no deliberate error
   state leaks past the scenario.
7. Check that no new page or console error appeared and that the Pie chart
   viewer element is still present on the page.

Expected:
- With a scratch column whose min aggregation is negative assigned as the angle column, the viewer renders a validation message containing "contains negative values" and the column name; with an all-zero scratch column under sum, the message contains "all values are 0"; clearing the angle column removes the message.

### Scenario 3: Category Map on the STARTED date column

Steps:
1. Set **Legend Visibility** to Always so the legend lists the categories.
2. Set **Category** to STARTED (a date column). Verify the **Category Map**
   defaults to year and read the legend category labels.
3. Set **Category Map** to month, wait for the render, and read the legend
   labels again.
4. Set **Category Map** to quarter, wait for the render, and read the legend
   labels a third time.
5. Verify all three label sets are non-empty and pairwise different — the map
   regroups the same dates into different category sets — with at most 12
   month entries and at most 4 quarter entries.
6. Revert: set **Category Map** back to year, **Category** back to RACE and
   **Legend Visibility** back to Auto.

Expected:
- With Category set to the STARTED date column the Category Map defaults to year; switching it to month and then quarter changes the legend category set each time — the three label sets are pairwise different, month has at most 12 entries, quarter at most 4 — and the category is reverted to RACE afterwards.

### Scenario 4: Grid color-coding on RACE reflected by the pie legend

Steps:
1. With **Category** on RACE, set **Legend Visibility** to Always and read the
   current legend color of the Asian entry.
2. Through the RACE grid column's categorical color coding, set the Asian
   category to red, refresh the viewer, and read the Asian legend entry color
   again. Verify it is now red and differs from the original color.
3. Clear the RACE color coding, refresh, and verify the Asian legend entry
   returns to its original color (round-trip).
4. Revert **Legend Visibility** to Auto.
5. Verify that no new page or console error appeared during the scenario and
   that the Pie chart viewer element is still present on the page.

Expected:
- Applying a categorical color to a RACE category through the grid column's color coding turns that category's pie legend entry to the requested color, and clearing the color coding restores the original legend color (round-trip).

### Scenario 5: Tooltip content follows the configured aggregations

Steps:
1. With **Category** on RACE and no angle or length column set, move the mouse
   onto a slice and read the tooltip. Verify its text names one of the RACE
   categories and includes a row-count digit.
2. Set **Segment Angle Column** to AGE with **Segment Angle Aggr Type** avg.
   Hover a slice again — verify the tooltip text now contains "avg(AGE)".
3. Set **Segment Length Column** to WEIGHT with **Segment Length Aggr Type**
   max. Hover a slice again — verify the tooltip text contains both
   "avg(AGE)" and "max(WEIGHT)".
4. Move the mouse to an empty corner of the chart and off the chart. Verify no
   new text is rendered into the tooltip and it is not visibly shown.
5. Restore both aggregation types to their captured defaults, then clear
   **Segment Angle Column** and **Segment Length Column** (standard pie
   restored, no validation message).

Expected:
- Hovering a slice with the default count aggregation shows a tooltip whose text names a RACE category and includes a row-count digit; after setting the angle column to AGE with avg the tooltip text contains "avg(AGE)"; after adding WEIGHT as the length column with max it additionally contains "max(WEIGHT)"; after the mouse moves to an empty corner and leaves the chart, no new text is rendered into the tooltip and it is not visibly shown.

## Automation notes

- Property drives: the category is `pie.props.categoryColumnName`; the angle
  surface is `pie.props.segmentAngleColumnName` /
  `pie.props.segmentAngleAggrType`; the length surface is
  `pie.props.segmentLengthColumnName` / `pie.props.segmentLengthAggrType`
  (clear a column with `''`); the date map is `pie.props.categoryMap`
  (`'year'`, `'month'`, `'quarter'`); the legend is
  `pie.props.legendVisibility`. Each set is followed by a 300-800 ms wait for
  the re-render.
- "Settle-gated canvas ink count" is `v.countCanvasPixels(page, 'Pie chart')
  .total` re-read every 300 ms until two consecutive counts differ by fewer
  than 200 pixels, up to 5 iterations. The avg-to-sum angle switch keeps the
  disc full, so its signal is the per-color histogram delta
  (`v.snapshotCanvasColors` + `v.diffCanvasColors`), which catches recoloring
  that the non-white total misses; the length-coding and clear steps change
  the inked area, so they use the plain ink count. All measured values are
  written to `console.log` on green runs so the fixed 500-pixel thresholds
  can be audited against live numbers; the thresholds are conservative wide
  margins, far below the deltas a real repaint produces. The steps' "clearly
  differ" / "clearly below" / "clearly above" margins are these fixed
  500-pixel thresholds.
- The length-coding drop is guaranteed by construction: the per-category
  length aggregates are min-max scaled, so the smallest category always maps
  to the minimum outer radius regardless of how close the aggregates are.
- Scenario 2 scratch columns are created with
  `df.columns.addNewFloat(name).init(...)` and removed in a `finally` block.
  The validation message renders into the viewer's `.d4-viewer-error`
  element; the assertion is a substring match on its text (message shapes:
  `<aggr>(<col>) contains negative values`, `<aggr>(<col>): all values are
  0`), so wording drift around the core phrase does not break the test.
- Legend labels are read with `v.readLegend(page, 'Pie chart')`
  (`[name="legend"] .d4-legend-item .d4-legend-value` under the viewer root).
  Scenario 4 reads a legend entry's color as the computed `color` style of
  its `.d4-legend-item` element; the grid color coding is applied with
  `df.col('RACE').meta.colors.setCategorical({'Asian': '#ff0000'})` and
  cleared by deleting the `.color-coding-categorical` and
  `.color-coding-type` column tags, with `viewer.invalidate()` after each
  change so the legend repaints.
- Scenario 1's measurements chain within the scenario (the sum-state count is
  the full-disc baseline for the length-coding drop), so the steps must run
  in order within one session.
- State isolation: the default angle and length aggregation types are captured
  from the viewer once, before any scenario changes them, and every clear
  restores them BEFORE clearing the column — clearing an angle/length column
  while a non-default aggregation is still set blanks the chart with a
  validation message. Each scenario opens with a precondition guard asserting
  the viewer's `.d4-viewer-error` text is empty, so a residual error from a
  previous scenario fails loudly at the scenario boundary instead of as a
  misleading downstream signal, and each scenario's teardown restores the
  properties it changed (Scenario 2 also restores the aggregation default and
  clears the angle column in its cleanup block, Scenario 3 reverts the
  Category Map, category and legend visibility, Scenario 4 reverts the color
  coding and legend visibility).
- Scenario 5 hovers by dispatching synthetic `mousemove` events at
  fraction-of-rect points on the pie canvas (the monolith's slice
  coordinates, starting at 0.65/0.4), reading the text of the
  `.d4-tooltip` element after a settle wait. The tooltip element is a
  page-level
  singleton that keeps its last text when hidden, so each hover first moves
  to an empty corner (resetting the hovered-slice state) and blanks the
  tooltip text — any text read afterwards was written by that hover. The
  aggregation substrings are the exact `<aggr>(<col>): <value>` format the
  tooltip builder emits; under the default count aggregation no aggregation
  line is added, so the first check asserts only a category name plus a
  digit. Honest reduction: synthetic mouse events populate the tooltip text
  but do not always flip the element to `display: block`, so the
  hovered-state checks assert on text content only, and the moved-away check
  asserts that nothing is re-rendered after blanking and that the element is
  not `display: block` — per-slice hover identification and the visual
  highlight color stay manual (see the paired manual checklist).

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/demog.csv"]
}
