---
feature: piechart
realizes_atlas:
  - piechart.cp.segment-click-select-filter
  - piechart.segment-click-select-syncs-dataframe
  - piechart.segment-click-filter-multi-category
realizes:
  - viewers.pie-chart
priority: p0
target_layer: playwright
coverage_type: regression
realized_as:
  - piechart-onclick-select-filter-spec.ts
related_bugs:
  - id: github-936
    status: fixed
  - id: GROK-17498
    status: fixed
  - id: GROK-2892
    status: fixed
  - id: GROK-19970
    status: fixed
expected_results:
  - anchor: "slice click selects exactly the category row count"
    expectation: "df.selection.trueCount EQUALS the clicked category's row count
      computed from the dataframe, and every selected row belongs to that one
      category"
  - anchor: "another slice switches the selection"
    expectation: "after clicking a different slice, df.selection.trueCount
      EQUALS the new category's computed row count and the selected category
      differs from the first one"
  - anchor: "clearing returns to zero"
    expectation: "df.selection.trueCount returns to 0 after the selection is
      cleared (round-trip)"
  - anchor: "clicking it selects exactly the rows with missing values"
    expectation: "df.selection.trueCount EQUALS col.stats.missingValueCount of
      the gapped category column and every selected row has a missing value in
      that column"
  - anchor: "Include Nulls off repaints the pie without the slice"
    expectation: "with Include Nulls off, the pie canvas per-color pixel
      histogram changes versus the Include-Nulls-on frame (delta above the
      noise floor) — the missing-values slice left the drawing"
  - anchor: "Include Nulls on restores it"
    expectation: "with Include Nulls back on, clicking the missing-values
      slice again selects exactly the missing-value rows —
      df.selection.trueCount EQUALS the missing-value count (round-trip)"
  - anchor: "slice click filters to the category row count"
    expectation: "in Filter mode, df.filter.trueCount EQUALS the clicked
      category's row count computed from the dataframe, and only that category
      remains among the filtered rows"
  - anchor: "Ctrl+click adds a second category to the sum"
    expectation: "df.filter.trueCount EQUALS the SUM of the two clicked
      categories' computed row counts, and exactly those two categories remain
      among the filtered rows"
  - anchor: "empty-area click restores the full count"
    expectation: "df.filter.trueCount returns to the full row count after a
      click on the empty canvas corner"
  - anchor: "pie click narrows the panel-filtered count further"
    expectation: "with a Filter Panel histogram filter (AGE 30-50) active, a
      slice click in Filter mode drops df.filter.trueCount BELOW the
      panel-only value and EQUALS the count of rows passing the panel filter
      AND belonging to the clicked category — AND-composition, not replacement"
  - anchor: "clearing the pie part returns to the panel-only value"
    expectation: "df.filter.trueCount returns EXACTLY to the panel-only value
      after the pie's click-filter is cleared — the Filter Panel filter
      survives"
  - anchor: "Reset filters restores the full count"
    expectation: "df.filter.trueCount returns to the full row count after the
      Filter Panel Reset filters button clears both filter sources"
  - anchor: "closing the pie chart with an active click-filter restores the full count"
    expectation: "df.filter.trueCount returns to the full row count after the
      pie chart viewer is closed while its click-filter is active — the viewer
      releases its filter on close"
  - anchor: "re-adding the viewer succeeds"
    expectation: "a new Pie chart viewer attaches to the view and the full
      row count is unchanged"
---

# Pie Chart — Segment Click Select and Filter Modes

## Purpose

Verify what a click on a pie slice does to the shared table in both On Click
modes. In Select mode a slice click must select exactly the rows of that
category; in Filter mode it must filter to exactly that category, Ctrl+click
must add a second category, and clicking empty space must restore all rows.
The missing-values slice (Include Nulls) must select exactly the rows with
missing values. A pie click-filter applied on top of a Filter Panel filter
must combine with it (AND), not replace it, and closing the viewer must
release its part of the filter. Every signal is a row count computed from the
data at run time, never a hard-coded number.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Pie chart to the current table view via the Toolbox viewer icon.
4. Set Category to RACE.
5. Set Row Source to All, so every slice stays visible while a click narrows
   the table (with the default Filtered source the pie redraws to only the
   surviving slices and a second category cannot be clicked).

## Scenarios

### Scenario 1: Select mode — slice click selects the category

Steps:
1. Set On Click to Select.
2. Click a slice — the rows of that category become selected; the selected row count equals the number of rows in that category, counted from the data.
3. Click a different slice — the selection switches to the new category; the selected row count equals the new category's row count.
4. Clear the selection — the selected row count returns to zero.

Expected:
- A slice click selects exactly the category row count
- Another slice switches the selection to its own category
- Clearing returns to zero

### Scenario 2: Missing-values slice — Include Nulls

Steps:
1. Add a scratch category column cloned from RACE with every tenth value blanked (no stock demog category column carries missing values) and set Category to it; keep On Click as Select and Include Nulls on.
2. Click the missing-values slice — the selected row count equals the number of rows with a missing value in the column, and every selected row has a missing value.
3. Turn Include Nulls off — the pie repaints without the missing-values slice (the drawn colors change).
4. Turn Include Nulls back on and click the missing-values slice again — it selects exactly the missing-value rows once more.
5. Remove the scratch column and set Category back to RACE.

Expected:
- Clicking the missing-values slice selects exactly the rows with missing values
- Include Nulls off repaints the pie without the slice
- Include Nulls on restores it (the slice is clickable again)

### Scenario 3: Filter mode — click, Ctrl+click, empty-area reset

Steps:
1. Set On Click to Filter.
2. Click a slice — the filtered row count equals that category's row count, and only that category remains visible.
3. Ctrl+click a different slice — the filtered row count equals the sum of the two categories' row counts.
4. Click the empty canvas corner — the filter resets and the filtered row count returns to the full row count.
5. Set On Click back to Select.

Expected:
- A slice click filters to the category row count
- Ctrl+click adds a second category to the sum
- The empty-area click restores the full count

### Scenario 4: Filter mode on top of a Filter Panel filter

Steps:
1. Open the Filter Panel and apply a histogram filter on AGE (window 30-50) — the filtered row count drops below the full count.
2. Set On Click to Filter and click a slice — the filtered row count drops further, to exactly the rows that pass the panel filter AND belong to the clicked category: the two filters combine, the pie click does not replace the panel filter.
3. Click the empty canvas corner — only the pie's part clears; the filtered row count returns exactly to the panel-only value.
4. On the Filter Panel click **Reset filters** — the filtered row count returns to the full count.

Expected:
- The pie click narrows the panel-filtered count further (AND-composition)
- Clearing the pie part returns to the panel-only value — the panel filter survives
- Reset filters restores the full count

### Scenario 5: Closing the viewer releases its filter

Steps:
1. Set On Click to Filter and click a slice — the filtered row count drops below the full count.
2. Close the pie chart viewer while its click-filter is active — the filtered row count returns to the full count: the viewer releases its filter on close.
3. Add the Pie chart viewer again — it attaches and the full row count is unchanged.

Expected:
- Closing the pie chart with an active click-filter restores the full count
- Re-adding the viewer succeeds

## Automation notes

- Every count is computed from the dataframe at run time: category sizes by
  scanning the category column, missing values via
  `col.stats.missingValueCount`, selection via `df.selection.trueCount`,
  filtering via `df.filter.trueCount`. Nothing is hard-coded; the clicked
  categories are the two most frequent ones, derived from the data.
- A slice click is the monolith's event triple — mousedown + mouseup + click
  dispatched on the pie canvas with `clientX`/`clientY` (plus `ctrlKey: true`
  on all three events for the additive Filter-mode click). The core hit-tests
  client coordinates directly, so no mousemove priming is required (verified
  live against dev: the triple alone selects/filters).
- A fixed canvas fraction cannot target a chosen slice (demog RACE is ~90%
  one category), so the click point per target category is computed from the
  data: slice angular spans mirror the pie's layout — per-category counts
  (nulls included only when Include Nulls is on), ascending-count order under
  the default by-value/asc sort, `startAngle` origin, clockwise in screen
  coordinates, radius capped by `maxRadius`. Candidate points inside the
  target span are then CALIBRATED in Select mode: a point is accepted only
  when clicking it selects exactly the target category's rows. Filter-mode
  scenarios click the calibrated points; with Row Source = All the geometry
  ignores the filter state, so points calibrated on the clean state stay
  valid under an active Filter Panel filter.
- Row Source is set to All in Setup: under the default Filtered source the
  pie redraws to only the surviving slices after a click-filter, which makes
  Ctrl+clicking a second category impossible (verified live).
- The missing-values slice (Scenario 2) needs a column with missing values;
  demog has none among its category columns (verified live on dev), so a
  scratch column `RACE_GAPS` is created via the JS API (clone of RACE, every
  10th value blanked) and removed in the step's `finally`. The pie legend
  does NOT list a '(no value)' item (verified live), so slice
  presence/absence is asserted as a per-color canvas pixel-histogram delta
  across the Include Nulls toggle, alongside the exact null-rows click
  selection.
- The empty-area reset click lands at `rect.left + 5, rect.top + 5` (outside
  the pie circle — a no-modifier miss clears the pie's click-filter).
- The Filter Panel histogram filter (Scenario 4) is applied via
  `grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50})`;
  the Reset filters button is
  `.d4-filter-group-header [name="icon-arrow-rotate-left"]`. The AND-composed
  expectation is computed from a snapshot of the panel-only filter mask.
- The viewer close (Scenario 5) uses the viewer's `close()` call — the same
  path the monolith's layout step used; the viewer is re-added through the
  Toolbox icon helper.

---
{
  "order": 6,
  "datasets": ["System:DemoFiles/demog.csv"]
}
