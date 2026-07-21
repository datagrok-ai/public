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
  - id: GROK-18484
    status: fixed
  - id: GROK-20033
    status: fixed
realized_as:
  - multi-axis-and-split-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      After enabling Multi Axis with 2 Y columns set, no console error or page
      error is raised (warnings delta == 0, page error count unchanged).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After adding a first split column while Multi Axis is active, no console
      error is raised and the chart canvas remains non-empty (github-2904
      partial repro path — chart does not go blank with one split column).
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      After adding a second split column (github-2904 full repro path — two
      split columns combined with multi-axis), the chart canvas is still
      non-empty and no console error is raised; a follow-up evaluate resolves,
      confirming the chart did not freeze or go blank.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      After setting 3 Y columns and triggering a Y-column edit in the property
      panel, the Y-column list length reported by the viewer's props is still 3
      (GROK-18484 regression guard — the edit did not silently reset the
      selection to one column).
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      The multi-axis column-selector search input is present in the DOM and is
      positioned inside the selector bounds (GROK-20033 — the input is not
      clipped or displaced outside the dialog/popup container); or if the search
      input is absent, no console error is raised (no-error floor fallback).
  - anchor: "Scenario 3 Step 2"
    expectation: >-
      After disabling Multi Axis (revert to single-axis), no console error is
      raised (no-error teardown).
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      After clearing all split columns, no console error is raised (no-error
      teardown confirmed; chart returns to single-series single-axis state).
---

# Line Chart — Multi-Axis and Split

## Setup

1. Close all open views.
2. Open the SPGI dataset (use `readDataframe` helper or open via the file browser).
3. Add a Line Chart viewer to the active table view (use `findViewer` helper or
   add via the Toolbox Viewers icon).
4. In the Line Chart property panel, set **X column** to `CAST Idea ID`
   (a numeric column present in SPGI).
5. Set **Y columns** to `Chemical Space X` and `Chemical Space Y` (two Y columns
   — the minimum required to exercise the Multi Axis path).
6. Record the current page-error count and `grok.shell.warnings` baseline before
   each action in the scenarios below (used as the no-error floor delta reference).

## Scenarios

### Scenario 1: Multi-axis with two split columns (github-2904 repro path)

Steps:
1. In the Line Chart property panel, under **Y axis**, enable **Multi Axis**
   (non-default lead action — toggle from off to on).
2. Confirm the chart renders each Y series on its own axis.
3. Assert: no console error or page error since baseline (Step 3 no-error floor —
   multi-axis enabled with 2 Y columns).
4. Under **Data**, set **Split** to `Stereo Category`
   (a categorical column in SPGI — first split column added while multi-axis is active).
5. Assert: no console error since the previous checkpoint AND the chart canvas is
   non-empty — a follow-up `page.evaluate` that checks the viewer element is
   visible resolves within timeout (Step 5 — github-2904 partial path, one split
   column).
6. Under **Data**, append a second split column: add `Activity Flag`
   (a second categorical column — github-2904 full repro path: two split columns
   combined with multi-axis).
7. Assert: the chart canvas remains non-empty (a follow-up evaluate resolves
   confirming the viewer element is still visible and rendered), and no console
   error has been raised since the previous checkpoint (Step 7 — github-2904 full
   repro: chart does NOT go blank with two split columns + multi-axis).

Expected:
- Enabling Multi Axis with 2 Y columns raises no console error (Step 3).
- Adding one split column while multi-axis is active keeps the chart canvas
  non-empty with no console error (Step 5 — github-2904 partial path).
- Adding a second split column still keeps the chart canvas non-empty and no
  console error is raised; the chart does NOT go blank (Step 7 — github-2904
  full repro path).

### Scenario 2: Y-column edit does not reset selection; column-selector search input position

Steps:
1. Starting from the state at the end of Scenario 1 (multi-axis on, two split
   columns), clear the split columns: under **Data**, set **Split** back to none.
   Reset Y columns to `Chemical Space X`, `Chemical Space Y`, and `Activity Class`
   (3 Y columns set).
2. Record the current Y-column count from the viewer's property object
   (`viewer.props.yColumnNames.length` — expected: 3).
3. In the property panel, open the **Y columns** multi-select input and trigger an
   edit: click into the selector, scroll or type to interact with the search field,
   then close the selector without removing any column (GROK-18484 repro path —
   an incidental edit that should not reset the Y-column selection).
4. Assert: `viewer.props.yColumnNames.length` is still 3 — the Y-column list was
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
- After an incidental Y-column selector edit, `viewer.props.yColumnNames.length`
  remains 3 — the Y-column list was not reset to 1 (Step 3 — GROK-18484 regression
  guard).
- The multi-axis column-selector search input is present and positioned inside the
  selector popup bounds (Step 5 — GROK-20033); or no console error if absent
  (no-error floor fallback).

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
