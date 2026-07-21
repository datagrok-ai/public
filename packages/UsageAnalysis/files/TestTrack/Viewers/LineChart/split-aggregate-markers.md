---
feature: linechart
realizes:
  - linechart.cp.setup-split-aggregate-markers
realizes_atlas:
  - linechart.cp.setup-split-aggregate-markers
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-17519
    status: fixed
  - id: GROK-19883
    status: fixed
  - id: GROK-20255
    status: fixed
realized_as:
  - split-aggregate-markers-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      After setting X and Y columns, no console error or page error is raised
      (grok.shell.warnings delta == 0, page error count unchanged).
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      After splitting by a categorical column (non-default lead), no console
      error is raised.
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      After setting an aggregation type and whiskers mode, no console error is
      raised.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      After setting marker type and a size-coding column (GROK-19883 repro
      path), no console error is raised.
  - anchor: "Scenario 1 Step 12"
    expectation: >-
      With 2+ split categories, a per-category average line renders (GROK-20255
      repro path) with no console error.
  - anchor: "Scenario 1 Step 14"
    expectation: >-
      After adding several additional split columns (GROK-17519 repro path), the
      page still responds: a follow-up JavaScript evaluate resolves within
      timeout, and no console error is raised.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      After switching X axis to logarithmic, no console error is raised.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      After reverting X axis back to linear, no console error is raised.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      After reverting split to none, no console error is raised (no-error
      teardown confirmed).
---

# Line Chart — Setup, Split, Aggregate, Markers

## Setup

1. Close all open views.
2. Open the SPGI dataset (use `readDataframe` helper or open via the file browser).
3. Add a Line Chart viewer to the active table view (use `findViewer` helper or
   add via the Toolbox Viewers icon).
4. Record the current page-error count and `grok.shell.warnings` baseline before
   each action in the scenarios below (used as the no-error floor delta reference).

## Scenarios

### Scenario 1: Connected-line trend, split, aggregation, whiskers, markers

Steps:
1. In the Line Chart property panel, set **X column** to `CAST Idea ID`
   (a numeric column present in SPGI).
2. Set **Y columns** to `Chemical Space X` (one column to establish a connected-line trend).
3. Verify that the chart builds — Step 4.
4. Assert: no console error or page error since baseline (warnings delta == 0,
   page error count unchanged). The canvas renders a connected trend line.
5. In the property panel under **Data**, set **Split** to `Stereo Category`
   (a categorical column in SPGI — non-default, lead action).
6. Assert: no console error since the previous checkpoint (Step 6 no-error floor).
7. Under **Data**, set **Aggregation** to `avg` (a non-default aggregation for
   duplicate X values).
8. Set **Whiskers** mode to `std err` (an error-bar mode).
   Assert: no console error since the previous checkpoint (Step 8 no-error floor).
9. Under **Markers**, set **Marker type** to `circle`.
10. Under **Markers**, set **Marker size** column to `Chemical Space Y`
    (exercises the size-coding path — GROK-19883 repro path).
    Assert: no console error since the previous checkpoint (Step 10 no-error floor).
11. Add a second Y column: set **Y columns** to `Chemical Space X` and
    `Chemical Space Y` (two Y columns, same split — gives 2+ split categories in
    the rendered chart).
12. Assert: the chart renders a per-category average line for each split category
    (GROK-20255 repro path) — no console error since the previous checkpoint
    (Step 12 no-error floor).
13. Add further split columns: in the property panel under **Data**, append
    split columns `R1`, `R2`, `R3` one at a time (GROK-17519 repro path —
    many splits).
14. After adding each additional split column, assert: a follow-up `page.evaluate`
    resolves (the page still responds) and no console error has been raised since
    the previous checkpoint. Final state: several split columns applied, chart
    stays responsive (Step 14 responsiveness + no-error floor).

Expected:
- The Line Chart builds a connected-line trend with X = `CAST Idea ID`,
  Y = `Chemical Space X`, with no console error (Step 4).
- Split by `Stereo Category` renders per-category series with no console error
  (Step 6).
- Aggregation to `avg` and whiskers mode `std err` are applied with no console
  error (Step 8).
- Marker type and size-coding column set with no console error (Step 10 —
  GROK-19883 path).
- With 2+ split categories, a per-category average line renders with no console
  error (Step 12 — GROK-20255 path).
- With many split columns added, the page still responds and no console error is
  raised (Step 14 — GROK-17519 path).

### Scenario 2: Mid-flow logarithmic X axis toggle and revert teardown

Steps:
1. Starting from the state at the end of Scenario 1 (split applied, markers set),
   record the current warning/error baseline.
2. In the property panel under **X axis**, switch **X axis type** to `Logarithmic`
   (non-default, mid-flow probe).
3. Assert: no console error since baseline (Step 3 no-error floor for log-axis
   toggle).
4. Switch **X axis type** back to `Linear`.
5. Assert: no console error since the previous checkpoint (Step 5 no-error floor
   — round-trip from log back to linear).
6. In the property panel under **Data**, set **Split** back to none (clear all
   split columns).
7. Assert: no console error since the previous checkpoint (Step 7 no-error teardown
   — revert split).

Expected:
- Switching X axis to Logarithmic raises no console error (Step 3).
- Reverting X axis back to Linear raises no console error (Step 5 — round-trip
  confirmed).
- Clearing all split columns (revert) raises no console error (Step 7 — no-error
  teardown).
