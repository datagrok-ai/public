---
feature: linechart
realizes:
  - viewers.line-chart
realizes_atlas:
  - linechart.cp.multi-axis-and-split
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: github-2904
    status: fixed
  - id: GROK-17835
    status: fixed
  - id: GROK-18484
    status: fixed
  - id: GROK-20033
    status: fixed
realized_as:
  - multi-axis-and-split-spec.ts
expected_results:
  - anchor: "S1: enable Multi Axis with 2 Y columns"
    expectation: >-
      After enabling Multi Axis with 2 Y columns set, no console error or page
      error is raised (warnings delta == 0, page error count unchanged).
  - anchor: "S1: first split column, chart stays non-empty"
    expectation: >-
      After adding a first split column while Multi Axis is active, no console
      error is raised and the chart canvas remains non-empty (github-2904
      partial repro path — chart does not go blank with one split column).
  - anchor: "S1: second split column, chart does not go blank"
    expectation: >-
      After adding a second split column (github-2904 full repro path — two
      split columns combined with multi-axis; the same surface GROK-17835
      regressed on), the chart canvas is still non-empty and no console error
      is raised; a follow-up evaluate resolves, confirming the chart did not
      freeze or go blank. Moving the pointer into the chart canvas (hover over
      the multi-axis + two-splits surface — the literal GROK-17835 repro
      gesture) raises no error either.
  - anchor: "S2: 3 Y columns survive an edit, not reset to 1"
    expectation: >-
      After setting 3 Y columns and triggering a Y-column edit in the property
      panel, the Y-column list length reported by the viewer's props is still 3
      (GROK-18484 regression guard — the edit did not silently reset the
      selection to one column).
  - anchor: "S2: search input is inside the selector bounds"
    expectation: >-
      The multi-axis column-selector search input is present in the DOM and is
      positioned inside the selector bounds (GROK-20033 — the input is not
      clipped or displaced outside the dialog/popup container); or if the search
      input is absent, no console error is raised (no-error floor fallback).
  - anchor: "Hide other charts"
    expectation: >-
      With Multi Axis on and 3 Y columns set, right-clicking the chart area and
      clicking 'Hide other charts' in the per-chart submenu (the group named
      after the chart's Y column) reduces the viewer's yColumnNames to exactly
      one column that was a member of the previous 3-column set (menu-to-prop),
      with no console or page error; the 3-column Y set is then restored so
      Scenario 3 starts from the documented state.
  - anchor: "S3: disable Multi Axis"
    expectation: >-
      After disabling Multi Axis (revert to single-axis), no console error is
      raised (no-error teardown).
  - anchor: "S3: clear all split columns"
    expectation: >-
      After clearing all split columns, no console error is raised (no-error
      teardown confirmed; chart returns to single-series single-axis state).
---

# Line Chart — Multi-Axis and Split

## Setup

1. Close all open views.
2. Open the spgi-100 dataset (use `readDataframe` helper or open via the file browser).
3. Add a Line Chart viewer to the active table view (use `findViewer` helper or
   add via the Toolbox Viewers icon).
4. In the Line Chart property panel, set **X column** to `CAST Idea ID`
   (a numeric column present in spgi-100).
5. Set **Y columns** to `Chemical Space X` and `Chemical Space Y` (two Y columns
   — the minimum required to exercise the Multi Axis path).
6. Record the current console/page error baseline before each action in the
   scenarios below (used as the no-error floor delta reference).

## Scenarios

### Scenario 1: Multi-axis with two split columns (github-2904 repro path)

Steps:
1. In the Line Chart property panel, under **Y axis**, enable **Multi Axis**
   (non-default lead action — toggle from off to on).
2. Confirm the chart renders each Y series on its own axis.
3. Assert: no console error or page error since baseline (Step 3 no-error floor —
   multi-axis enabled with 2 Y columns).
4. Under **Data**, set **Split** to `Stereo Category`
   (a categorical column in spgi-100 — first split column added while multi-axis is active).
5. Verify no console error since the previous checkpoint AND the chart canvas is
   non-empty — the viewer element remains visible and rendered (Step 5 —
   github-2904 partial path, one split column).
6. Under **Data**, append a second split column: add `Series`
   (a second categorical column — github-2904 full repro path: two split columns
   combined with multi-axis).
7. Assert: the chart canvas remains non-empty (a follow-up evaluate resolves
   confirming the viewer element is still visible and rendered), and no console
   error has been raised since the previous checkpoint (Step 7 — github-2904 full
   repro: chart does NOT go blank with two split columns + multi-axis).
8. Move the pointer into the center of the chart canvas (hover over the
   multi-axis + two-splits surface — the literal GROK-17835 repro gesture) and
   assert: no console or page error is raised by the hover handling.

Expected:
- Enabling Multi Axis with 2 Y columns raises no console error (Step 3).
- Adding one split column while multi-axis is active keeps the chart canvas
  non-empty with no console error (Step 5 — github-2904 partial path).
- Adding a second split column still keeps the chart canvas non-empty and no
  console error is raised; the chart does NOT go blank (Step 7 — github-2904
  full repro path).
- Hovering over the chart in that state raises no console or page error
  (Step 8 — GROK-17835 repro gesture).

### Scenario 2: Y-column edit does not reset selection; column-selector search input position

Steps:
1. Starting from the state at the end of Scenario 1 (multi-axis on, two split
   columns), clear the split columns: under **Data**, set **Split** back to none.
   Reset Y columns to `Chemical Space X`, `Chemical Space Y`, and `TPSA`
   (3 Y columns set).
2. Record the current Y-column count from the viewer's properties (expected: 3).
3. In the property panel, open the **Y columns** multi-select input and trigger an
   edit: click into the selector, scroll or type to interact with the search field,
   then close the selector without removing any column (GROK-18484 repro path —
   an incidental edit that should not reset the Y-column selection).
4. Verify the viewer's Y-column count is still 3 — the Y-column list was
   NOT silently reset to 1 by the edit interaction (Step 3 — GROK-18484 regression
   guard; the point is the edit did not shrink the list, not a round-trip of the
   setter value).
5. While the multi-axis column-selector is open (re-open **Y columns** selector),
   locate the search input element inside the selector popup.
   Assert: the search input element is present in the DOM and its bounding rect is
   contained within the selector popup container bounds — `inputRect.top >=
   containerRect.top && inputRect.bottom <= containerRect.bottom` (GROK-20033 —
   the search input is positioned correctly inside the selector, not displaced or
   clipped outside). If the search input is not present, assert no console error
   has been raised (no-error floor fallback for environments where the search input
   is absent).

Expected:
- After an incidental Y-column selector edit, the viewer's Y-column count
  remains 3 — the Y-column list was not reset to 1 (Step 3 — GROK-18484 regression
  guard).
- The multi-axis column-selector search input is present and positioned inside the
  selector popup bounds (Step 5 — GROK-20033); or no console error if absent
  (no-error floor fallback).

### Scenario 2b: Hide other charts via the per-chart context menu

Steps:
1. Starting from the state at the end of Scenario 2 (multi-axis on, 3 Y
   columns, no split columns), confirm Multi Axis is on and record the
   3-column Y set.
2. Right-click the chart area. With Multi Axis on, the context menu contains a
   submenu TITLED WITH A Y COLUMN'S NAME (e.g. **Chemical Space X**). Open that submenu and click **Hide other charts**.
   The submenu only appears while Multi Axis is on with 2+ Y columns.
3. Assert: the viewer's Y column list is reduced to exactly one column, and
   that column is a member of the recorded 3-column set (menu-to-prop signal);
   no console or page error is raised.
4. Restore the recorded 3-column Y set and confirm it reads back, so
   Scenario 3 starts from the documented "multi-axis on, Y columns set" state.

Expected:
- 'Hide other charts' reduces the Y columns to the one chart it was invoked on,
  with no console error; the Y set restore round-trips.

### Scenario 3: Revert multi-axis and split (no-error teardown)

Steps:
1. Starting from the state at the end of Scenario 2 (multi-axis on, Y columns set),
   record the current warning/error baseline.
2. Under **Y axis**, disable **Multi Axis** (toggle off — revert to single-axis).
   Assert: no console error since baseline (Step 2 no-error teardown).
3. Under **Data**, clear all split columns (set **Split** to none if any remain).
4. Assert: no console error since the previous checkpoint (Step 4 — final no-error
   teardown; chart returns to single-series single-axis state with all splits
   cleared).

Expected:
- Disabling Multi Axis raises no console error (Step 2 — no-error teardown).
- Clearing all split columns raises no console error (Step 4 — no-error teardown
  confirmed).

## Automation notes

Setup: the console/page error baseline is the page-error count plus the
`grok.shell.warnings` length; the dataset is opened via the `readDataframe`
helper and the viewer handle obtained via the `findViewer` helper. The Y-column
count is read as `viewer.props.yColumnNames.length`.

Scenario 2: the "Select columns" dialog reached from the Y-columns `...`
editor renders its column LIST on a canvas — there are no DOM rows or checkboxes,
and clicks on the canvas do not toggle columns. A checkbox-toggle + OK mutation
is therefore not reachable in automation; open / search / cancel is the deepest
scriptable editor interaction. The GROK-18484 observable (the edit must NOT
silently reset the selection) is still exercised: the dialog's `3 checked`
identity label confirms the correct row was opened, and the exact `yColumnNames`
set read back after cancel confirms the list was not shrunk to one column. A
checkbox-driven mutation stays a human-side variation of this scenario.

Scenario 2b: right-clicking the chart area appends a per-chart group to the
context menu, named after the Y column the click position resolves to; the
group's 'Hide other charts' item collapses the Y column list to that single
column. The group name is dynamic (a column name), so the menu item is located
by its visible 'Hide other charts' label, not by a name attribute.
