---
feature: pcplot
realizes_atlas:
  - pcplot.cp.transformation-and-filter-integrity
realizes:
  - viewers.pc-plot
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - pcplot-transformation-filter-spec.ts
related_bugs:
  - id: GROK-18091
    status: fixed
  - id: GROK-17306
    status: open
  - id: GROK-18489
    status: fixed
  - id: github-972
    status: fixed
expected_results:
  - anchor: "Scenario 1 (GROK-17306) — Reset filters restores filter, keeps selection"
    expectation: "df.filter.trueCount restores to the full row count (the pre-filter
      value) after Reset filters"
  - anchor: "Reset filters restores filter, keeps selection"
    expectation: "df.selection.trueCount is unchanged — equal to its pre-reset
      value; resetting filters must NOT clear the selection"
  - anchor: "Scenario 2 (GROK-18489) — second filter after DateTime color split works"
    expectation: "df.filter.trueCount changes to a value below the full row count
      (the second filter application takes effect)"
  - anchor: "second filter after DateTime color split works"
    expectation: "df.filter.trueCount changes again when the filter is adjusted a
      second time (the plot does NOT freeze the count)"
  - anchor: "Scenario 3 (github-972) — histogram column change does not reset PC filter"
    expectation: "df.filter.trueCount on the shared DataFrame is unchanged after the
      histogram column change (the PC plot filter is not reset)"
  - anchor: "Scenario 4 (GROK-18091) — aggregation + Filter Panel close, no broken state"
    expectation: "No error is raised and no broken state (empty grid or blank PC
      plot) occurs after closing the Filter Panel; grok.shell.warnings /
      page-error delta is 0"
---

# PC Plot — Transformation and Filter/Selection Integrity

## Purpose

Guard for four filter/selection/transformation integrity bugs in the PC Plot.
Three are fixed (GROK-18091, GROK-18489, github-972); GROK-17306 (reset filters
clears the selection when a transformation is present) is still OPEN, so
Scenario 1 fails loudly until it is fixed. All bugs involve the PC Plot's
interaction with the Filter Panel, row selection, or other viewers sharing the
same DataFrame. Each scenario is judged by the filtered or selected row
count, or by the absence of errors — not by how the canvas looks.

## Setup

1. Close all views.
2. Open the demog dataset (`System:DemoFiles/demog.csv`).
3. Go to its table view.
4. Add a PC Plot viewer.
5. Record the full row count.

## Scenarios

### Scenario 1: Reset filters does not clear row selection (GROK-17306)

Regression guard: with a PC plot transformation present and some rows selected,
resetting the Filter Panel must restore the filtered row count to the full row
count but must NOT change the selected-row count.

Steps:
1. Open the Filter Panel (show the Filters group on the table view).
2. In the Context Panel > Data, set **Transformation** to the aggregation script (a plain
   text field — no dialog needed):
   `[{"#type":"GroupAggregation","aggType":"key","colName":"SEX"},{"#type":"GroupAggregation","aggType":"pivot","colName":"DIS_POP"},{"#type":"GroupAggregation","aggType":"avg","colName":"WEIGHT"}]`.
   This is the GROK-17306 condition. Because the pivot replaces the PC-plot axes, the
   filter is applied on the Filter Panel, not on an in-chart axis slider.
3. Apply a range filter through the **Filter Panel** on the AGE column (narrow it to a
   sub-range, e.g. 30–50). Record the filtered row count.
4. Confirm the filtered row count is below the full count — an active Filter Panel filter
   is in place with the transformation applied, so the reset has something to restore.
5. Select the first ten rows; record the selected-row count.
6. Confirm the selected-row count is above zero — a selection exists. Click "Reset filters"
   in the Filter Panel header, then verify the filtered row count restores to the
   full row count (a round-trip from the filtered value).
7. Verify the selected-row count equals the value recorded in Step 5 — the selection is
   UNCHANGED; the reset did not clobber it. Then clear the Transformation to restore the baseline.

Expected:
- With a transformation present, the filtered row count restores to the full row count after
  Reset filters, having first dropped below it under the range filter (a genuine round-trip).
- The selected-row count is unchanged — equal to its pre-reset value; resetting filters must NOT clear the selection.

### Scenario 2: Second filter after DateTime color split works (GROK-18489)

Regression guard: with a DateTime column set as the Color column, PC Plot
slider-based filtering must keep working after a filter reset — the second
filter application must change the filtered row count.

Steps:
1. Set the PC plot's color column (Context Panel > Color > Color Column) to a DateTime column (STARTED).
2. Drag a per-axis range-slider handle to a value window (AGE axis, for example)
   to apply a first filter — verify the filtered row count drops below the full
   count (the first filter took effect).
3. Reset the filter via the Filter Panel "Reset filters" button; verify the
   filtered row count returns to the full count (filter fully cleared after reset).
4. Drag the same range-slider handle to a different value window to apply a
   second filter on the PC plot.
5. Verify the filtered row count is below the full row count — the second filter
   application takes effect.
6. Verify that repeating the drag/release once more changes the filtered row
   count again (the plot does NOT freeze the count at the previous value).

Expected:
- The filtered row count changes to a value below the full row count (the second filter application takes effect).
- The filtered row count changes again when the filter is adjusted a second time (the plot does NOT freeze the count).

### Scenario 3: Histogram column change does not reset PC plot filter (github-972)

Regression guard: when a PC Plot and a Histogram share the same DataFrame,
changing the histogram's active column must not affect the filtered row count.

Steps:
1. Add a Histogram viewer to the same view.
2. Apply a per-axis range-slider filter on the PC plot (AGE axis) so that the
   filtered row count drops below the full count; record the filtered row count.
3. Confirm the filtered row count is below the full count (active filter in place).
4. Change the Histogram's value column from AGE to HEIGHT.
5. Verify the filtered row count equals the value recorded in Step 2 — the PC plot
   filter is NOT reset by the histogram column change; the shared table filter is unchanged.

Expected:
- The filtered row count on the shared table is unchanged after the histogram column change (the PC plot filter is not reset).

### Scenario 4: Aggregation transformation with Filter Panel open — no broken state (GROK-18091)

Regression guard: with the Filter Panel open, adding an aggregation transformation
on the PC plot then closing the Filter Panel must not clear the grid or blank the
plot. The check is that nothing breaks: no errors, and no blank grid or plot.

Steps:
1. Note the current page-error count.
2. Open the Filter Panel (click the Filters toolbar button or use keyboard shortcut).
3. On the PC plot, set an aggregation transformation (Context Panel > Data >
   Transformation, e.g. an average over AGE).
4. Close the Filter Panel (click its close/X button).
5. Verify no new page errors have been raised since the baseline.
6. Verify the grid is not empty: the table still has rows (the transformation may
   change the visible row set but the table must not become blank).

Expected:
- No error is raised and no broken state (empty grid or blank PC plot) occurs after closing the Filter Panel; the warnings / page-error delta is 0.

## Automation notes

- Views are closed via `grok.shell.closeAll()`; the demog table is opened with
  the shared table-open helper and the viewer is added through the Toolbox
  icon. The full row count is `df.rowCount`.
- "Filtered row count" is `df.filter.trueCount`; "selected-row count" is
  `df.selection.trueCount`. The Scenario 1 selection is made via
  `df.selection.init((i) => i < 10)`.
- The transformation entries in Scenarios 1 and 4 are set via
  `pc.props.transformation` — the same value the Context Panel > Data >
  Transformation text field writes. The Scenario 1 Filter Panel range filter
  is applied via `tv.getFiltersGroup().updateOrAdd({type: 'histogram',
  column: 'AGE', min: 30, max: 50})`, and "Reset filters" is the
  `[name="icon-arrow-rotate-left"]` icon in the Filter Panel header.
- Scenario 2: the color column is set via `pc.props.colorColumnName =
  'STARTED'`; the range filters in Scenarios 2 and 3 are applied by dragging
  the real `axis-slider-AGE` min/max handles with synthetic mouse events.
- Scenario 3: the histogram is added via `tv.addViewer('Histogram')` and its
  column changed via `hist.props.valueColumnName = 'HEIGHT'`.
- Scenario 4: the error baseline counts only page and console errors naming
  the transformation / PC-Plot / aggregation surface (the shared dev server
  emits unrelated console noise); the broken-state check reads `df.rowCount`
  staying above zero and the PC Plot element still present.
