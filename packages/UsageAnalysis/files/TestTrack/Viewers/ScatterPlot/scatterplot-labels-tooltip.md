---
feature: scatterplot
realizes_atlas:
  - scatterplot.cp.labels-tooltip
realizes:
  - viewers.scatter-plot
priority: p2
target_layer: playwright
coverage_type: regression
realized_as:
  - scatterplot-labels-tooltip-spec.ts
related_bugs:
  - id: GROK-17251
    status: fixed
expected_results:
  - anchor: "Tooltip inherited from the table"
    expectation: "Hovering the aim point — computed from the eleventh row of demog's
      WEIGHT and HEIGHT values, the columns on the X and Y axes — shows the
      tooltip, and its text lists the table's tooltip columns with the values of
      the row the viewer reports as hovered: every value read from the tooltip
      equals the value the table holds for THAT row in the same column, so the
      check depends on neither a fixed column list nor a fixed row"
  - anchor: "Custom tooltip column list"
    expectation: "With Show Tooltip set to the custom mode, Row Tooltip set to AGE
      and SEX, and Data Values set to not add anything, the tooltip text lists
      exactly those two configured values and nothing else"
  - anchor: "Custom tooltip column list"
    expectation: "With Show Tooltip set to the do-not-show mode, hovering the same
      point produces no tooltip text"
  - anchor: "Custom tooltip column list"
    expectation: "Restoring Show Tooltip to inherit from the table and clearing the
      custom row tooltip brings the full inherited tooltip back"
  - anchor: "Labels for the selected rows"
    expectation: "The Shift-drag over the cloud selects rows: the selected-row count
      rises above zero (the product precondition the label check depends on)"
  - anchor: "Labels for the selected rows"
    expectation: "With Label Columns set to AGE and Show Labels For set to Selected,
      selecting a band of points changes the overlay rendering against the
      no-selection baseline (the labels are drawn) and raises no console or page
      error; the label text itself is canvas-only and is not asserted"
  - anchor: "Labels for the selected rows"
    expectation: "Clearing the selection restores the overlay rendering toward the
      baseline, and clearing Label Columns and resetting Show Labels For returns
      the viewer to the setup configuration"
  - anchor: "Labels plus filtering on a datetime axis"
    expectation: "Unchecking every Stereo Category value except the first one in the
      Stereo Category filter card narrows the table: the filtered row count
      drops below the dataset's full row count (the product precondition the
      crash guard depends on)"
  - anchor: "Labels plus filtering on a datetime axis"
    expectation: "On the chemical dataset with X set to the datetime column Whole
      blood assay 2 Date, Y set to Stereo Category, a label column configured
      and Show Labels For set to Selected, applying a categorical filter from
      the filter panel and then selecting rows raises no console or page error —
      in particular no Unsupported operation: Infinity.ceil() error
      (GROK-17251)"
  - anchor: "Labels plus filtering on a datetime axis"
    expectation: "After resetting the filter and the selection, the filtered row
      count returns to the full row count of the dataset"
---

# Scatter Plot — Marker Labels and Tooltip

## Purpose

Verifies the Scatter plot's labelling and tooltip surface: the tooltip
inherited from the table, a custom tooltip restricted to a chosen column list,
the do-not-show mode, marker labels rendered for the selected rows, and the
crash guard that combines labels with filtering on a datetime axis.

The tooltip is real DOM text, so the custom-column check is a direct read of
what the tooltip lists. Marker labels, by contrast, are drawn on the viewer's
overlay canvas and contribute no DOM text at all, so the label steps assert a
rendering change plus an error-free floor and never the label text.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Scatter plot to the current table view via the Toolbox viewer icon.
4. Using the on-viewer column selectors, set X to WEIGHT and Y to HEIGHT.
5. Confirm the starting state: **Show Tooltip** inherits from the table, no
   label columns configured.

## Scenarios

### Scenario 1: Tooltip inherited from the table

The aim point for this file is computed from the eleventh row of demog — its
WEIGHT and HEIGHT values, the columns the setup put on the X and Y axes — so the
same point is hovered again in Scenario 2. demog is dense, so the marker actually
hit there is usually a NEIGHBOURING row, not the eleventh: the viewer's own
reported hovered row is what the tooltip must agree with. Asserting against the
eleventh row literally would fail for the wrong reason.

Steps:
1. Read the eleventh row's values for WEIGHT and HEIGHT and hover the pointer
   over the marker at those two coordinates; wait for the tooltip to appear.
2. Verify the aim landed on a marker: the viewer reports a hovered row, and that
   marker's screen position is within a small fraction of the viewport of the
   point aimed at.
3. Read the tooltip text.
4. Verify it lists the table's tooltip columns with the HOVERED row's values —
   compare each value shown against the value the table holds for that row in the
   same column, so the check depends on neither a fixed column list nor a fixed
   row.

Expected:
- Hovering the aim point shows the tooltip, and every value it lists equals the table's value for the hovered row in the same column

### Scenario 2: Custom tooltip column list

Steps:
1. Open the viewer settings and go to the **Tooltip** section.
2. Set **Show Tooltip** to the custom-tooltip mode.
3. Set **Row Tooltip** to the two columns AGE and SEX.
4. Set **Data Values** so that no data values are added.
5. Hover the reference row's marker as in Scenario 1 and read the tooltip
   text.
6. Verify it lists exactly the two configured values and nothing else.
7. Set **Show Tooltip** to the do-not-show mode, hover the reference row's
   marker again, and verify no tooltip text appears.
8. Revert: set **Show Tooltip** back to inherit from the table, clear
   **Row Tooltip**, and restore **Data Values**.
9. Hover the reference row's marker and verify the full inherited tooltip is
   back.

Expected:
- The custom tooltip lists exactly the configured columns and nothing else
- The do-not-show mode produces no tooltip text
- Restoring the tooltip settings brings the inherited tooltip back

### Scenario 3: Labels for the selected rows

Steps:
1. Note the current browser console-error and page-error counts, and measure
   the plot's overlay rendering with nothing selected — the baseline.
2. Open the viewer settings and, in the **Labels** section, set
   **Label Columns** to AGE.
3. Set **Show Labels For** to **Selected**.
4. Shift-drag a rectangle over a populated part of the cloud to select a band
   of points; verify the selected-row count rises above zero.
5. Verify the overlay rendering changed against the baseline — labels are
   drawn for the selected rows. The label text is painted on the canvas and is
   not asserted.
6. Verify no new console or page error appeared.
7. Revert: click empty background to clear the selection and verify the
   overlay rendering moves back toward the baseline; then clear
   **Label Columns** and set **Show Labels For** back to its default.

Expected:
- The Shift-drag selects rows: the selected-row count rises above zero
- Selecting a band of points with labels configured for Selected changes the overlay rendering and raises no console or page error
- Clearing the selection restores the overlay rendering toward the baseline, and clearing the label settings returns the viewer to the setup configuration

### Scenario 4: Labels plus filtering on a datetime axis

This scenario runs on a chemical dataset that carries a real datetime column,
which is what makes the label placement math degenerate when filtering empties
the drawable range.

Steps:
1. Close all open views and open `System:AppData/Chem/tests/spgi-100.csv`;
   wait for the table view to load. Record the dataset's full row count.
2. Add a Scatter plot via the Toolbox viewer icon.
3. Using the on-viewer selectors, set X to **Whole blood assay 2 Date** (a
   datetime column) and Y to **Stereo Category**.
4. In the viewer settings **Labels** section, set **Label Columns** to **Id**
   and **Show Labels For** to **Selected**.
5. Note the current console-error and page-error counts.
6. Click the **Filters** icon in the Toolbox to open the filter panel and wait
   for the filter cards to build.
7. In the **Stereo Category** filter card — the categorical card for the same
   column that sits on the Y axis — uncheck every category except the first
   one listed, leaving exactly one category selected. Verify the filtered row
   count drops below the dataset's full row count.
8. Shift-drag a rectangle over the plot to select a band of the remaining
   points.
9. Verify no new console or page error appeared across the filtering and
   selection steps — in particular no `Unsupported operation: Infinity.ceil()`
   error (GROK-17251 guard).
10. Revert: reset the filter from the filter panel, clear the selection, and
    verify the filtered row count returns to the dataset's full row count.

Expected:
- Leaving exactly one Stereo Category checked in its filter card drops the filtered row count below the dataset's full row count
- Filtering and then selecting with labels configured on a datetime X axis raises no console or page error, in particular no Infinity.ceil() error
- Resetting the filter and the selection returns the filtered row count to the dataset's full row count

## Automation notes

- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Scatter plot')`; the viewer is
  added via the Toolbox icon `[name="icon-scatter-plot"]`. Resolve the viewer
  root as the `[name="viewer-Scatter-plot"]` element not inside a `.d4-dialog`.
  On-viewer column selectors are lowercase
  (`div-column-combobox-x|y|color|size`) and open on a synthetic `mousedown`,
  never on a synthetic `.click()`; real typing plus `Enter` commits.
- Hovering: aim with
  `sp.worldToScreen(df.col('WEIGHT').get(10), df.col('HEIGHT').get(10))` plus
  the canvas bounding-rect origin (the helper returns canvas-local coordinates).
  Do NOT assume row 10 is the row hovered — demog is dense enough that the hit
  test lands on a neighbour. Read `df.mouseOverRowIdx` and grade the tooltip
  against that row, and separately prove the aim landed on the intended marker by
  checking the hovered marker's screen position against the aim point within a
  small fraction of the viewport. The tooltip-vs-table comparison reads the
  hovered row's value per column from the dataframe, formatted as the tooltip
  formats it, and never hard-codes either the column list or the values. The tooltip is the document-level
  `.d4-tooltip`; `sp.getRowTooltip(rowIdx)` is the JS-side way to read the same
  content without hovering and is a useful cross-check.
- Tooltip rows: `[name="prop-show-tooltip"]`
  (`do not show | inherit from table | show custom tooltip`),
  `[name="prop-row-tooltip"]`, `[name="prop-data-values"]`
  (`Do not add | Data values only | Merge`). Note
  `[name="prop-show-column-names"]` exists in BOTH the Labels and the Tooltip
  categories — scope through the unique `prop-view-*` cells.
- Labels are canvas-only: with label columns configured and rows selected, zero
  new DOM text nodes appear inside the viewer and only the overlay canvas
  rendering moves. Measure `canvas[name="overlay"]` as a settle-gated non-white
  pixel fraction; calibrate the threshold live at spec-authoring time rather
  than carrying a settled constant, and whitelist the browser's
  `willReadFrequently` advisory that the sampling itself triggers. Label rows
  are `[name="prop-label"]` (a multi-column editor whose `...` button opens the
  column picker) and the menu leaves `div-Labels---Label-Columns` /
  `div-Labels---Show-Labels-For---Selected`; the JS fallback is
  `sp.props.labelColumnNames` / `showLabelsFor`.
- Selection is driven with synthetic `MouseEvent`s on `canvas[name="canvas"]`
  (a bubbling `mousedown` + several `mousemove`s + `mouseup` with `shiftKey`),
  which actuate on this viewer; the signal is
  `grok.shell.tv.dataFrame.selection.trueCount`, never a hard-coded count. A
  click on the empty background clears the selection (gated by
  `resetSelectionOnBackgroundClick`); a click on a marker does not.
- Scenario 4: `grok.dapi.files.readCsv` does NOT name the table after the file
  (tables arrive as `Table`, `Table (2)`), so nothing may key on a table name —
  use the view's dataFrame handle and read `df.rowCount` for the full count.
  The Toolbox `[name="div-section--Filters"]` may be present but invisible
  (`offsetParent === null`); guard on that and fall back to
  `tv.getFiltersGroup()`. The filter panel builds slowly on this dataset — poll
  for the `.d4-filter` card count instead of a fixed wait. The target card is
  the categorical filter for `Stereo Category` — locate it by its column name
  among the `.d4-filter` cards rather than by position. Its checkbox list is
  canvas-drawn, so the category labels are NOT DOM text: read the retained
  category from the column itself
  (`df.col('Stereo Category').categories[0]`) and drive the narrowing through
  the filter API or a real coordinate click on the card, never by matching
  label text. The assert is on `df.filter.trueCount` dropping below
  `df.rowCount` — the live reference for this dataset and column is 100 rows
  narrowing to 38, quoted as a reference only, never as the asserted value. The filter is reset with
  `[name="viewer-Filters"] [name="icon-arrow-rotate-left"]` (no confirmation
  dialog). Console guards subtract a baseline and whitelist infra noise
  unrelated to the viewer (the shared dev server's WebSocket reconnect errors,
  the Claude runtime container timeout and the WebGPU `powerPreference`
  advisory); the `Infinity.ceil()` class must always reach the guard.

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
