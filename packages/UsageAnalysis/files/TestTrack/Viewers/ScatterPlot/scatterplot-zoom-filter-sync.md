---
feature: scatterplot
realizes_atlas:
  - scatterplot.cp.zoom-filter-sync
  - scatterplot-axis-type-drives-filter-out-invalid
realizes:
  - viewers.scatter-plot
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - scatterplot-zoom-filter-sync-spec.ts
related_bugs:
  - id: GROK-13141
    status: fixed
  - id: GROK-15496
    status: fixed
  - id: GROK-18379
    status: fixed
  - id: GROK-19153
    status: fixed
expected_results:
  - anchor: "Zoom drives the table filter, reset restores it"
    expectation: "With Zoom and Filter set to 'filter by zoom', zooming into the
      plot drops the table's filtered row count below the full row count, and a
      second zoom step drops it further; resetting the view restores the
      filtered row count to the full row count (round-trip). Every reading is
      polled until the count stops changing, because the filter update trails
      the zoom gesture"
  - anchor: "Pack and Zoom keeps the zoom filtering resettable"
    expectation: "After zoom-filtering the table and then selecting the 'pack and
      zoom by filter' mode, resetting the view returns the filtered row count to
      the full row count — the filtering never becomes stuck (GROK-13141 guard)"
  - anchor: "External filtering drives the viewport in zoom by filter mode"
    expectation: "With Zoom and Filter set to 'zoom by filter', restricting the
      X-axis column through the Filter Panel to a band taken from that column's
      own minimum and maximum drops the filtered row count below the full row
      count and narrows the viewer's viewport well beyond the tolerance that
      counts two rectangles as the same view; removing the filter restores both
      the full row count and the viewport rectangle recorded before the filter.
      The viewport is read as the viewer's own viewport rectangle, never as the
      viewport property, and the baseline rectangle is verified non-empty so the
      narrowing assert cannot pass vacuously"
  - anchor: "Filter Panel reset clears the scatter plot's contribution"
    expectation: "After zoom-filtering the table, the status bar reports filtered
      rows; clicking Reset filter in the Filter Panel returns the filtered row
      count to the full row count AND the status-bar filtered-rows indicator is
      gone, so nothing still reports a scatter plot filter (GROK-15496 guard)"
  - anchor: "Axis type switch on a datetime axis keeps the applied filter"
    expectation: "With STARTED (a datetime column) on the X axis and the Filter
      Panel SEX filter narrowed to the single category F, choosing Logarithmic
      under Properties > X > X Axis Type in the plot's context menu actually
      switches the axis type — the X axis type reads logarithmic after the step
      — AND the filtered row count recorded before the switch is unchanged after
      it, so the switch does not reset the applied filtering (GROK-18379 guard).
      The axis-type change is verified BEFORE the filter count is graded; an
      unchanged filter over an axis that did not switch is not a pass"
  - anchor: "Filter Out Invalid removes the rows a logarithmic axis cannot draw"
    expectation: "With an X column whose values are all strictly positive and a Y
      column carrying a known number of non-positive values, switching the Y
      axis to logarithmic leaves the filtered row count at the full row count
      while Filter Out Invalid is off; turning Filter Out Invalid on drops the
      filtered row count by exactly the number of non-positive Y values, and
      turning it off restores the full row count. With the Y axis back to
      linear, Filter Out Invalid leaves the filtered row count at the full row
      count — the rows are removed only because a logarithmic axis cannot draw
      them. The number of non-positive values is counted from the data, never
      hard-coded"
  - anchor: "Large jitter with a logarithmic axis does not filter rows"
    expectation: "On the spgi-100 dataset with X = CAST Idea ID and Y = Idea ID (the
      live auto-pick; Idea ID is an identifier column whose values are strictly
      positive, so a logarithmic Y axis is meaningful) and with Jitter Size and
      Jitter Size Y raised to sixteen and seventeen, switching Y Axis Type to
      logarithmic leaves the filtered row count at the full row count — jitter
      combined with a logarithmic axis never silently filters data (GROK-19153
      guard)"
---

# Scatter Plot — Zoom and Filter Synchronization

## Purpose

Verifies that the scatter plot's Zoom and Filter coupling moves the table
filter in both directions and always stays reversible: zooming narrows the
filtered rows, resetting the view restores them, the Pack and Zoom mode does
not trap the filtering, and a Filter Panel reset clears the viewer's
contribution together with the indicator that reports it. Two settings that
historically dropped or invented filtering are guarded as well — switching an
axis to a logarithmic type on a datetime axis, and combining large jitter with
a logarithmic axis. The one case where a logarithmic axis is supposed to move
the filter is covered too: the rows it cannot draw are removed only when Filter
Out Invalid is on. Every outcome is read from the table's filtered row count,
which the product computes, and each reading is polled until it settles
because the filter update trails the zoom gesture.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table
   view to load.
3. Add a Scatter plot to the current table view via the Toolbox viewer icon.
4. On the viewer, set the X column selector to WEIGHT and the Y column
   selector to HEIGHT (the auto-pick already matches on demog, but the
   auto-selection order is not contractual).
5. Record the full row count of the table with no filtering applied.

## Scenarios

### Scenario 1: Zoom drives the table filter, reset restores it

Steps:
1. Open the viewer settings (gear icon on the Scatter plot title bar) and, in
   the **Data** section, confirm **Zoom and Filter** reads `filter by zoom`
   (this is the value a freshly added viewer carries); set it explicitly if it
   does not.
2. Read the table's filtered row count — it equals the full row count.
3. Zoom in one step with the mouse wheel over the middle of the plot.
4. Wait for the filtered row count to stop changing, then verify it dropped
   below the full row count.
5. Zoom in a second step with the mouse wheel.
6. Wait for the count to settle again and verify it dropped further.
7. Reset the view with the **Reset View** command in the plot's context menu.
8. Wait for the count to settle and verify it is back at the full row count
   (round-trip).
9. Leave **Zoom and Filter** at `filter by zoom` and the view reset — the
   starting state for the next scenario.

Expected:
- Each zoom step drops the filtered row count below the full row count, the second step drops it further, and resetting the view restores the full row count
- Every count is read only after it stops changing — the filter update trails the zoom

### Scenario 2: Pack and Zoom keeps the zoom filtering resettable

The Pack and Zoom mode used to capture the filtering the scatter plot had
applied, leaving the filtered-out rows unrecoverable. This scenario drives
that sequence and asserts that the filtering still round-trips.

Steps:
1. Zoom in with the mouse wheel until the filtered row count settles below the
   full row count.
2. In the viewer settings **Data** section, set **Zoom and Filter** to
   `pack and zoom by filter`.
3. Reset the view with the **Reset View** command in the plot's context menu.
4. Wait for the filtered row count to settle and verify it is back at the full
   row count — the scatter plot's filtering is not stuck.
5. Set **Zoom and Filter** back to `filter by zoom` and reset the view
   (round-trip revert).

Expected:
- After the Pack and Zoom mode is selected, resetting the view returns the filtered row count to the full row count — the filtering never becomes unrecoverable

### Scenario 3: External filtering drives the viewport in zoom by filter mode

The other modes read the viewport and write the table filter. This one runs the
coupling the other way round: the table filter is the input, and the viewer
answers by fitting its viewport to the rows that survive. The filter is applied
to a column that is on an axis, so the viewport has to move for the scenario to
mean anything.

Steps:
1. In the viewer settings **Data** section, set **Zoom and Filter** to
   `zoom by filter` and confirm the property reads back that value.
2. Record the viewer's viewport rectangle and the filtered row count with
   nothing filtered — the count equals the full row count and the rectangle has
   a non-zero width and height.
3. Open the Filter Panel and restrict WEIGHT — the column on the X axis — to a
   band strictly inside its own minimum and maximum, taken from the column's
   values rather than written down.
4. Wait for the filtered row count to settle and verify it is below the full row
   count and above zero.
5. Verify the viewport narrowed well beyond the tolerance used elsewhere for
   "the same view" — the viewer followed the filter.
6. Remove the filter. Wait for the count to settle and verify it is back at the
   full row count and the viewport is back at the rectangle recorded in step 2.
7. Set **Zoom and Filter** back to `filter by zoom` and close the Filter Panel
   (round-trip revert).

Expected:
- With Zoom and Filter set to `zoom by filter`, restricting the X-axis column through the Filter Panel narrows the viewport well beyond the same-view tolerance while the filtered row count drops below the full row count
- Removing the filter restores both the full row count and the viewport rectangle recorded before the filter was applied

### Scenario 4: Filter Panel reset clears the scatter plot's contribution

Steps:
1. Open the Filter Panel from the Toolbox **Filters** section and wait for its
   filter cards to finish building.
2. Zoom in with the mouse wheel until the filtered row count settles below the
   full row count.
3. Verify the status bar shows the filtered-rows indicator, reporting that the
   table is filtered.
4. Click **Reset filter** in the Filter Panel header (it applies immediately,
   with no confirmation).
5. Wait for the filtered row count to settle and verify it is back at the full
   row count.
6. Verify the status-bar filtered-rows indicator is gone — nothing still
   reports a scatter plot filter on an unfiltered table.
7. Close the Filter Panel and reset the view (round-trip revert).

Expected:
- Zoom-filtering makes the status bar report filtered rows
- Reset filter in the Filter Panel returns the filtered row count to the full row count and removes the filtered-rows indicator from the status bar

### Scenario 5: Axis type switch on a datetime axis keeps the applied filter

On a datetime axis the settings-panel **X Axis Type** row is greyed out — the
axis type is meaningless for a datetime column, and that greying is the fixed
behaviour of a separate guard. The context menu does not mirror that state, so
this scenario actuates the switch from the context menu and confirms the axis
type really changed before grading the filter.

Steps:
1. On the viewer, set the X column selector to STARTED — a datetime column.
2. Open the Filter Panel and, in the SEX filter, restrict the table to the
   single category F.
3. Wait for the filtered row count to settle and record it — it is below the
   full row count.
4. Right-click the plot and choose **Properties... > X > X Axis Type >
   Logarithmic** from the context menu (the settings-panel row for this
   property is greyed out while a datetime column is on the axis).
5. Verify the X axis type now reads logarithmic — the switch actually took
   effect. If it did not, stop: the rest of this scenario cannot be graded.
6. Wait for the render to settle and verify the filtered row count is
   unchanged from the value recorded in step 3 — the axis-type switch does not
   reset the applied filtering.
7. Set **X Axis Type** back to `linear` through the same context-menu path,
   clear the SEX filter, close the Filter Panel, and set the X column selector
   back to WEIGHT (round-trip revert).

Expected:
- Choosing Logarithmic from the context menu actually switches the axis type on a datetime axis — the axis type reads logarithmic after the step
- With the axis type confirmed switched, the filtered row count is unchanged from the value recorded before the switch — the applied filtering survives it

### Scenario 6: Filter Out Invalid removes the rows a logarithmic axis cannot draw

A logarithmic axis cannot draw a row whose value on that axis is zero, negative
or not finite. What happens to such a row is decided by **Filter Out Invalid**:
while it is off the row is simply skipped when drawing, and while it is on the
row is removed from the table filter. This scenario drives that coupling in both
directions on a fixture built for the purpose, because every numeric column of
demog is strictly positive and would never exercise it.

Steps:
1. Add two columns to the demog table for the duration of this scenario: an X
   probe whose values are all strictly positive, and a Y probe whose values are
   strictly positive except for a known share that are negative. Count the
   non-positive Y values from the data.
2. Set the X column selector to the X probe and the Y column selector to the Y
   probe.
3. In the viewer settings **Data** section, confirm **Filter Out Invalid** is
   off, and record the filtered row count — it equals the full row count.
4. In the **Y** section, set **Y Axis Type** to `logarithmic`. Wait for the
   filtered row count to settle and verify it is unchanged — a logarithmic axis
   on its own does not filter anything.
5. Turn **Filter Out Invalid** on. Wait for the count to settle and verify it
   dropped by exactly the number of non-positive Y values counted in step 1.
6. Turn **Filter Out Invalid** off. Wait for the count to settle and verify it is
   back at the full row count.
7. Set **Y Axis Type** back to `linear` and turn **Filter Out Invalid** on again.
   Wait for the count to settle and verify it stays at the full row count — with
   a linear axis the same values are drawable, so nothing is removed.
8. Turn **Filter Out Invalid** off, set the X and Y column selectors back to
   WEIGHT and HEIGHT, and remove both probe columns (round-trip revert).

Expected:
- With Filter Out Invalid off, switching the Y axis to logarithmic leaves the filtered row count at the full row count
- Turning Filter Out Invalid on drops the filtered row count by exactly the number of non-positive Y values, and turning it off restores the full row count
- With the Y axis back to linear, Filter Out Invalid leaves the filtered row count at the full row count — the rows are removed only because a logarithmic axis cannot draw them

### Scenario 7: Large jitter with a logarithmic axis does not filter rows

Jitter offsets the drawn markers; on a logarithmic axis the offset markers
used to be treated as undrawable and the corresponding rows were silently
filtered out. This scenario runs on the compact spgi-100 dataset, where the
original report was raised.

Steps:
1. Open the spgi-100 dataset: `System:AppData/Chem/tests/spgi-100.csv` — wait
   for the table view to load. Do not rely on the table's name; the loader
   does not name a table after its file.
2. Add a Scatter plot to this table view via the Toolbox viewer icon and set
   the X column selector to **CAST Idea ID** and the Y column selector to
   **Idea ID** — the two columns the viewer auto-picks on this dataset. The Y
   column carries strictly positive identifier values, which is what makes a
   logarithmic Y axis meaningful here.
3. Record the full row count and verify the table is unfiltered.
4. Open the viewer settings and, in the **Marker** section, set **Jitter
   Size** to 16 and **Jitter Size Y** to 17.
5. In the **Y** section, set **Y Axis Type** to `logarithmic`.
6. Wait for the render to settle and verify the filtered row count still
   equals the full row count — nothing was filtered out.
7. Set **Y Axis Type** back to `linear`, set both jitter values back to 0, and
   close this table view (round-trip revert).

Expected:
- With large jitter values and a logarithmic Y axis, the filtered row count stays at the full row count — the combination never silently filters rows

## Automation notes

- Narrowing a categorical filter in these steps is driven through the Filter
  Panel's filter-group API rather than by clicking the card's canvas. The guard
  needs a DETERMINISTIC surviving category set, and the card's checkbox list is
  canvas-drawn: a coordinate click can toggle *a* category but cannot choose
  *which* one. Where a guard only needs "exactly one category left, whichever it
  is", the real coordinate click is used instead — see the labels-tooltip
  scenario, which does exactly that. The filter narrowing here is setup for the
  graded signal, never the signal itself.

- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Scatter plot')`; the viewer is
  added via the Toolbox icon `[name="icon-scatter-plot"]` (a synthetic
  `.click()` works). Resolve the viewer root as the
  `[name="viewer-Scatter-plot"]` element that is NOT inside a `.d4-dialog`.
- The single grading signal in this scenario file is
  `grok.shell.tv.dataFrame.filter.trueCount`, compared against
  `dataFrame.rowCount` — never against a hard-coded number, and never a read
  back of a property that was just set.
- LOAD-BEARING TIMING: the filter update is debounced behind the zoom
  gesture — right after a wheel event the count is still the previous value
  and only settles later. Every reading in every scenario must poll until two
  consecutive reads agree (with a generous ceiling); a single fixed wait will
  flake.
- Zoom is a synthetic `WheelEvent` (negative `deltaY`) over
  `canvas[name="canvas"]`; Alt+drag also zooms. Reset goes through the
  `[name="div-Reset-View"]` context-menu leaf, NOT a double-click: a
  double-click resets only when it lands clear of a marker, and over the dense
  centre of the demog cloud it reads as a point double-click and does nothing —
  see the select-and-zoom scenario, which owns that behaviour and grades it.
  `sp.props.viewport` stays `null` — read `sp.viewport` when a viewport
  cross-check is wanted.
- Zoom and Filter is the `[name="prop-zoom-and-filter"]` row (click the view
  cell, set the revealed `select`'s `.value`, dispatch `change`) or the
  `[name="div-Properties...---Data---Zoom-and-Filter---<choice>"]` menu leaf.
  The mode value is spelled `pack and zoom by filter`, and the four choices
  are `no action | filter by zoom | zoom by filter | pack and zoom by filter`.
- Scenario 3 — the filtered column must be one of the axis columns. The viewport
  is fitted to the surviving rows, so a filter on a column that is not drawn
  moves it only as far as that column happens to correlate with the axes, and
  the narrowing assert would rest on an accident of the data. Filtering WEIGHT,
  which is on X, makes the viewport move by construction. The band comes from
  the column's own minimum and maximum through the filter-group range filter
  (`DG.FILTER_TYPE.HISTOGRAM` with `min` / `max`), the same filter-group route
  the categorical filters in this file use. The viewport is `sp.viewport`;
  `sp.props.viewport` is null on this viewer and is not a substitute.
- Scenario 4 — the filter indicator: the `?` icon named in the original report
  does not exist. The readable substitutes are the Filter Panel header's
  `.d4-filter-indicator` text (the count of Filter-Panel filters, which stays
  at its own value while a scatter plot filter is active) and the status bar
  `[name="span-filtered"]`, which carries an "N filtered rows" text and is
  ABSENT from the DOM while the table is unfiltered. Assert presence before
  reading text; absence after the reset is the assertable half of the guard.
  The Filter Panel reset button is
  `[name="viewer-Filters"] [name="icon-arrow-rotate-left"]`; it also resets the
  viewer's viewport. The Toolbox `[name="div-section--Filters"]` entry can be
  present but invisible — guard on `offsetParent` and fall back to
  `tv.getFiltersGroup()`; the panel builds slowly, so poll for its
  `.d4-filter` card count instead of waiting a fixed time.
- The on-viewer column selectors are `[name="div-column-combobox-x"]`,
  `-y`, `-color`, `-size` — lowercase role names, scoped to the viewer root
  (the property panel reuses the same names). They open on a synthetic
  `mousedown` on `.d4-column-selector-column` (never a synthetic `.click()`),
  then real typing plus `Enter` commits; the popup grid is canvas-rendered, so
  column names are not readable as DOM text.
- Property rows: `[name="prop-x-axis-type"]` / `[name="prop-y-axis-type"]`
  (choice), `[name="prop-jitter-size"]` / `[name="prop-jitter-size-y"]`
  (slider textbox — focus, select all, type, `Enter`). `prop-min` / `prop-max`
  are duplicated across X and Y, so scope through the unique
  `[name="prop-view-x-min"]`-style cells when a range is needed.
- Scenario 5 — ACTUATION ROUTE IS LOAD-BEARING. With a datetime column on X,
  the settings-panel row `[name="prop-x-axis-type"]` is greyed out (inline
  `opacity: 0.5` on the row — the only readable disabled signal this panel
  has); driving that row can silently do nothing, and then "the filter did not
  change" passes vacuously over an axis that never switched. The context menu
  does NOT mirror the greying, so the route to use is the menu leaf
  `[name="div-Properties...---X---X-Axis-Type---Logarithmic"]` (leaves are
  clickable without expanding the parent submenu). BEFORE grading the filter
  count, read the axis-type property before and after the step and assert it
  actually flipped to `logarithmic`. If neither the menu route nor any other UI
  route flips it on a datetime axis, then the GROK-18379 guard is not
  verifiable through the UI on this build: record an honest coverage gap and
  ESCALATE it. Do NOT fall back to setting the property from the JS API and
  present that as the UI guard, and do NOT report an unchanged filter over an
  unswitched axis as a passing check. (The property does remain settable from
  JS while the UI control is greyed — which is exactly why the axis-type read
  must gate the assert.)
- Scenario 6 — the fixture carries the whole point. Every numeric column of
  demog is strictly positive, so the coupling cannot be observed on the dataset
  as it ships; both probe columns are added at the start of the scenario and
  removed in a `finally` block, the same way the Legend scenario file builds its
  empty-value probe. The X probe must be free of missing values as well, because
  Filter Out Invalid removes a row that is undrawable on EITHER axis — with a
  stock demog column on X the drop would be larger than the non-positive count
  by the number of that column's missing values. The number of non-positive Y
  values is counted from the data at run time, never written into the spec.
  **Filter Out Invalid** is the checkbox row `[name="prop-filter-out-invalid"]`
  in the **Data** category, whose editor is
  `input[name="prop-view-filter-out-invalid"]` and whose `checked` state mirrors
  the property; the context menu carries the same toggle as
  `[name="div-Properties...---Data---Filter-Out-Invalid"]`. The `Filter` group of
  the context menu does NOT contain it — it holds only Axes Follow Filter, Show
  Filtered Out Points and Show Only Filtered Selected Points.
- Scenario 7 — the column choice is NOT free. The X and Y columns are the live
  auto-pick on spgi-100 (`CAST Idea ID`, `Idea ID`), which is the configuration
  the guard was reconned under, and the guard needs a Y column whose values are
  strictly positive so that a logarithmic Y axis is meaningful rather than
  degenerate. Confirm at spec-authoring time that the chosen Y column's minimum
  is above zero; if it is not, pick another strictly positive numeric column of
  the dataset, name it explicitly in the spec, and say why. A Y column with
  non-positive values would make rows undrawable on a log axis for a legitimate
  reason and would confound the guard.
- Console-error guards used alongside these counts must filter the shared dev
  server's unrelated WebSocket reconnect noise, the WebGPU `powerPreference`
  warning, and the `willReadFrequently` canvas notice that pixel sampling
  itself provokes.

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
