---
feature: linechart
realizes:
  - viewers.line-chart
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
  - anchor: "S1: set X and Y columns (connected trend)"
    expectation: >-
      After setting X and Y columns, no console error or page error is raised
      (grok.shell.warnings delta == 0, page error count unchanged).
  - anchor: "split by Stereo Category"
    expectation: >-
      After splitting by a categorical column (non-default lead), no console
      error is raised.
  - anchor: "aggregation avg + whiskers std err"
    expectation: >-
      After setting an aggregation type and whiskers mode, no console error is
      raised.
  - anchor: "marker type + size-coding column"
    expectation: >-
      After setting marker type and a size-coding column (GROK-19883 repro
      path), no console error is raised.
  - anchor: "S1: second Y column with split"
    expectation: >-
      With 2+ split categories, a per-category average line renders (GROK-20255
      repro path) with no console error.
  - anchor: "S1: add R1/R2/R3 splits, page stays responsive"
    expectation: >-
      After adding several additional split columns (GROK-17519 repro path), the
      page still responds: a follow-up JavaScript evaluate resolves within
      timeout, and no console error is raised.
  - anchor: "S2: X axis logarithmic"
    expectation: >-
      After switching X axis to logarithmic, no console error is raised.
  - anchor: "S2: X axis back to linear"
    expectation: >-
      After reverting X axis back to linear, no console error is raised.
  - anchor: "S2: clear all split columns"
    expectation: >-
      After reverting split to none, no console error is raised (no-error
      teardown confirmed).
  - anchor: "S1: set X and Y columns"
    expectation: >-
      After setting X and Y columns the canvas paints a connected trend line:
      the non-white pixel total exceeds the 28000-px render floor (a blank
      frame paints ~0 px; a -1 canvas fault also fails), guarded by
      canvas-dimension checks that protect the floor's calibration.
  - anchor: "S1: split by Stereo Category"
    expectation: >-
      Splitting by Stereo Category repaints the plot from one line to
      per-category lines: the settle-gated per-color canvas diff
      (snapshot, setter, diff) exceeds the repaint floor (deltaPx > 1000).
  - anchor: "S1: aggregation avg + whiskers std err"
    expectation: >-
      Setting aggregation avg with whiskers mode Avg | ±StError adds
      error-bar ink on top of the averaged lines: the settle-gated per-color
      canvas diff exceeds the repaint floor (deltaPx > 300).
  - anchor: "S1: marker type + size-coding column"
    expectation: >-
      Setting circle markers with a size-coding column adds marker ink at the
      data points: the settle-gated per-color canvas diff exceeds the repaint
      floor (deltaPx > 300).
---

# Line Chart — Setup, Split, Aggregate, Markers

## Setup

1. Close all open views.
2. Open the spgi-100 dataset (use `readDataframe` helper or open via the file browser).
3. Add a Line Chart viewer to the active table view (use `findViewer` helper or
   add via the Toolbox Viewers icon).
4. Record the current console/page error baseline before each action in the
   scenarios below (used as the reference for verifying no new errors appear).

## Scenarios

### Scenario 1: Connected-line trend, split, aggregation, whiskers, markers

Steps:
1. In the Line Chart property panel, set **X column** to `CAST Idea ID`
   (a numeric column present in spgi-100).
2. Set **Y columns** to `Chemical Space X` (one column to establish a connected-line trend).
3. Verify that the chart builds — Step 4.
4. Assert: no console error or page error since baseline (no new warnings or
   page errors appear). The canvas renders a connected trend line.
5. In the property panel under **Data**, set **Split** to `Stereo Category`
   (a categorical column in spgi-100 — non-default, lead action).
6. Assert: no console error since the previous checkpoint (Step 6 — verify no
   new errors appear).
7. Under **Data**, set **Aggregation** to `avg` (a non-default aggregation for
   duplicate X values).
8. Set **Whiskers** mode to the standard-error error-bar mode (the property
   enum literal is `Avg | ±StError`).
   Assert: no console error since the previous checkpoint (Step 8 — verify no
   new errors appear).
9. Under **Markers**, set **Marker type** to `circle`.
10. Under **Markers**, set **Marker size** column to `Chemical Space Y`
    (exercises the size-coding path — GROK-19883 repro path).
    Assert: no console error since the previous checkpoint (Step 10 — verify no
    new errors appear).
11. Add a second Y column: set **Y columns** to `Chemical Space X` and
    `Chemical Space Y` (two Y columns, same split — gives 2+ split categories in
    the rendered chart).
12. Assert: the chart renders a per-category average line for each split category
    (GROK-20255 repro path) — no console error since the previous checkpoint
    (Step 12 — verify no new errors appear).
13. Add further split columns: in the property panel under **Data**, append
    split columns `R1`, `R2`, `R3` one at a time (GROK-17519 repro path —
    many splits).
14. After adding each additional split column, verify the page still responds
    (the chart remains responsive) and no console error has been
    raised since the previous checkpoint. Final state: several split columns
    applied, chart stays responsive (Step 14 — chart remains responsive and no
    new errors appear).

Expected:
- The Line Chart builds a connected-line trend with X = `CAST Idea ID`,
  Y = `Chemical Space X`, with no console error (Step 4).
- Split by `Stereo Category` renders per-category series with no console error
  (Step 6).
- Aggregation to `avg` and whiskers mode `Avg | ±StError` (standard error) are
  applied with no console error (Step 8).
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
3. Assert: no console error since baseline (Step 3 — verify no new errors
   appear after the log-axis toggle).
4. Switch **X axis type** back to `Linear`.
5. Assert: no console error since the previous checkpoint (Step 5 — no new
   errors appear; round-trip from log back to linear).
6. In the property panel under **Data**, set **Split** back to none (clear all
   split columns).
7. Assert: no console error since the previous checkpoint (Step 7 — no new
   errors appear during the revert of the split).

Expected:
- Switching X axis to Logarithmic raises no console error (Step 3).
- Reverting X axis back to Linear raises no console error (Step 5 — round-trip
  confirmed).
- Clearing all split columns (revert) raises no console error (Step 7 — no-error
  teardown).

## Automation notes

Scenario 1 (canvas-delta measurement): the split, aggregation/whiskers, and
marker steps each carry a canvas repaint assert in addition to the no-error
floor. A per-color canvas histogram is snapshotted, a settle-precheck diff
proves the previous repaint has drained (near-zero delta), the property is set,
and the measured diff (deltaPx) must exceed a conservative floor — the split
repaints most of the plot ink, whiskers add error-bar ink, and size-coded
markers add ink at every data point. The initial trend render (Step 4) is
asserted as a non-white pixel total above the 28000-px render floor.

Setup / Step 14: the console/page error baseline is the page-error count plus
the `grok.shell.warnings` length; the Step 14 responsiveness check is a
follow-up `page.evaluate` that must resolve within timeout.
