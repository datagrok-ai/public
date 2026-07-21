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
  - anchor: "Scenario 1 Step 6"
    expectation: "df.filter.trueCount restores to the full row count (the pre-filter
      value) after Reset filters"
  - anchor: "Scenario 1 Step 7"
    expectation: "df.selection.trueCount is unchanged — equal to its pre-reset
      value; resetting filters must NOT clear the selection"
  - anchor: "Scenario 2 Step 5"
    expectation: "df.filter.trueCount changes to a value below the full row count
      (the second filter application takes effect)"
  - anchor: "Scenario 2 Step 6"
    expectation: "df.filter.trueCount changes again when the filter is adjusted a
      second time (the plot does NOT freeze the count)"
  - anchor: "Scenario 3 Step 5"
    expectation: "df.filter.trueCount on the shared DataFrame is unchanged after the
      histogram column change (the PC plot filter is not reset)"
  - anchor: "Scenario 4 Step 5"
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
same DataFrame. Each scenario asserts a named df-state signal
(df.filter.trueCount or df.selection.trueCount) or a no-error floor, not a
canvas pixel.

## Setup

1. Close all views: `await grok.shell.closeAll()`.
2. Open the demog dataset: `const df = await grok.data.files.openTable('System:DemoFiles/demog.csv')`.
3. Get the table view: `const view = grok.shell.tableView(df.name)`.
4. Add a PC Plot viewer: `const pc = view.addViewer('PC Plot')`.
5. Record the full row count: `const fullCount = df.rowCount`.

## Scenarios

### Scenario 1: Reset filters does not clear row selection (GROK-17306)

Regression guard: with a PC plot transformation present and some rows selected,
resetting the Filter Panel must restore df.filter.trueCount to the full row count
but must NOT change df.selection.trueCount.

Steps:
1. Open the Filter Panel (show the Filters group on the table view).
2. In the Context Panel > Data, set **Transformation** to the aggregation script (a plain
   text field — no dialog needed):
   `[{"#type":"GroupAggregation","aggType":"key","colName":"SEX"},{"#type":"GroupAggregation","aggType":"pivot","colName":"DIS_POP"},{"#type":"GroupAggregation","aggType":"avg","colName":"WEIGHT"}]`.
   This is the GROK-17306 condition. Because the pivot replaces the PC-plot axes, the
   filter is applied on the Filter Panel, not on an in-chart axis slider.
3. Apply a range filter through the **Filter Panel** on the AGE column (narrow it to a
   sub-range, e.g. 30–50). Capture `const filteredCount = df.filter.trueCount`.
4. Confirm `filteredCount < fullCount` — an active Filter Panel filter is in place with the
   transformation applied, so the reset has something to restore.
5. Select some rows via a JS-API call: `df.selection.init((i) => i < 10)`; capture
   `const selCount = df.selection.trueCount`.
6. Confirm `selCount > 0` — a selection exists. Click "Reset filters" in the Filter Panel
   header, then assert `df.filter.trueCount === fullCount` — the filter restores to the
   full row count (a round-trip from filteredCount).
7. Assert `df.selection.trueCount === selCount` — the selection is UNCHANGED; the reset
   did not clobber it. Then clear the Transformation to restore the baseline.

Expected:
- With a transformation present, df.filter.trueCount restores to the full row count after
  Reset filters, having first dropped below it under the range filter (a genuine round-trip).
- df.selection.trueCount is unchanged — equal to its pre-reset value; resetting filters must NOT clear the selection.

### Scenario 2: Second filter after DateTime color split works (GROK-18489)

Regression guard: with a DateTime column set as the Color column, PC Plot
slider-based filtering must keep working after a filter reset — the second
filter application must change df.filter.trueCount.

Steps:
1. Set the PC plot's color column to a DateTime column (STARTED):
   `pc.setOptions({color: 'STARTED'})`.
2. Drag a per-axis range-slider handle to a value window (AGE axis, for example)
   to apply a first filter — assert `df.filter.trueCount < fullCount` to confirm
   the first filter took effect.
3. Reset the filter via the Filter Panel "Reset filters" button; assert
   `df.filter.trueCount === fullCount` (filter fully cleared after reset).
4. Drag the same range-slider handle to a different value window to apply a
   second filter on the PC plot.
5. Assert `df.filter.trueCount < fullCount` — the second filter application
   takes effect (the count drops below full row count).
6. Assert that repeating the drag/release once more changes df.filter.trueCount
   again (the plot does NOT freeze the count at the previous value).

Expected:
- df.filter.trueCount changes to a value below the full row count (the second filter application takes effect).
- df.filter.trueCount changes again when the filter is adjusted a second time (the plot does NOT freeze the count).

### Scenario 3: Histogram column change does not reset PC plot filter (github-972)

Regression guard: when a PC Plot and a Histogram share the same DataFrame,
changing the histogram's active column must not affect df.filter.trueCount.

Steps:
1. Add a Histogram viewer to the same view: `const hist = view.addViewer('Histogram')`.
2. Apply a per-axis range-slider filter on the PC plot (AGE axis) so that
   `df.filter.trueCount < fullCount`; capture `const filteredCount = df.filter.trueCount`.
3. Confirm filteredCount < fullCount (active filter in place).
4. Change the Histogram's column from AGE to HEIGHT via `hist.setOptions({valueColumnName: 'HEIGHT'})`.
5. Assert `df.filter.trueCount === filteredCount` — the PC plot filter is NOT reset
   by the histogram column change; the shared DataFrame filter is unchanged.

Expected:
- df.filter.trueCount on the shared DataFrame is unchanged after the histogram column change (the PC plot filter is not reset).

### Scenario 4: Aggregation transformation with Filter Panel open — no broken state (GROK-18091)

Regression guard: with the Filter Panel open, adding an aggregation transformation
on the PC plot then closing the Filter Panel must not clear the grid or blank the
plot. Assert a no-error / no-freeze floor.

Steps:
1. Record the page-error count baseline: capture `const errsBefore = (await page.evaluate(() => window.__grokErrors?.length ?? 0))`.
2. Open the Filter Panel (click the Filters toolbar button or use keyboard shortcut).
3. On the PC plot, open the Data context menu and apply an aggregation transformation
   (e.g. right-click > Data > Transformation, enter a GroupAggregation expression).
4. Close the Filter Panel (click its close/X button).
5. Assert page-error delta: `(await page.evaluate(() => window.__grokErrors?.length ?? 0)) === errsBefore`.
6. Assert the grid is not empty: `df.rowCount > 0` (the transformation may change the
   visible row set but the DataFrame must not become blank).

Expected:
- No error is raised and no broken state (empty grid or blank PC plot) occurs after closing the Filter Panel; grok.shell.warnings / page-error delta is 0.
