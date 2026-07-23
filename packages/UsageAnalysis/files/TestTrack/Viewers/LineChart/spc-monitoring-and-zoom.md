---
feature: linechart
realizes_atlas:
  - linechart.cp.spc-monitoring-and-zoom
realizes:
  - viewers.line-chart
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20126
    status: fixed
realized_as:
  - spc-monitoring-and-zoom-spec.ts
expected_results:
  - anchor: "S1: enable SPC — page stays responsive, no freeze (GROK-20126)"
    expectation: >-
      After enabling Statistical Process Control, a follow-up page.evaluate
      resolves within timeout (the page has not frozen) and the page/console
      error count is unchanged since baseline.
  - anchor: "S1: disable SPC — no-error teardown"
    expectation: >-
      Reverting the SPC toggle restores the prior state with no new page/console
      error (teardown is asserted, not left unchecked).
  - anchor: "S2: Reset View clears the wheel-zoom but keeps explicit X Min/Max"
    expectation: >-
      Reset View resets the wheel-zoom (the zoomed and reset-view events fire)
      but the explicit X Min/Max stay at the values set in Step 3 — Reset View
      clears the zoom, not the configured axis bounds — and no new page/console
      error is raised.
---

# Line Chart — SPC Monitoring and Zoom

## Setup

1. Close all open views.
2. Open the spgi-100 dataset (use `readDataframe` helper or open via the file browser).
3. Add a Line Chart viewer to the active table view (use `findViewer` helper or
   add via the Toolbox Viewers icon).
4. In the Line Chart property panel, set **X column** to `CAST Idea ID` and
   **Y columns** to `Chemical Space X` (single Y column — required precondition
   for SPC: no split, no multi-axis).
5. Record the current console/page error baseline before Scenario 1.

## Scenarios

### Scenario 1: SPC toggle — no-freeze guard (GROK-20126)

Steps:
1. Verify the current configuration has no split column and multi-axis is off
   (the SPC gating precondition: single Y, no split, no multi-axis).
2. In the property panel, navigate to the **SPC** section.
3. Enable **Statistical Process Control** (toggle SPC on).
4. Immediately after the toggle, run a lightweight scripted check that the page
   still responds.
5. Verify the responsiveness check completes within the configured timeout AND
   no new console or page error has been raised since baseline
   (GROK-20126 repro: the no-freeze guard — Step 5).
6. Disable SPC (revert to off) and verify the no-error teardown: SPC reads back
   off and no new console/page error was raised since the previous checkpoint.

Expected:
- After enabling SPC, the page remains responsive: the responsiveness check
  completes and no new console or page error is raised (Step 5 — GROK-20126
  no-freeze guard).
- Disabling SPC reverts cleanly: SPC reads back off and no new console or page
  error is raised (Step 6 — no-error teardown).

### Scenario 2: Reset View clears the zoom, not the explicit X Min/Max

Reset View resets an interactive zoom; it does NOT clear axis bounds the user
has configured. Those are separate: X Min/X Max are settings, the zoom is view
state.

Steps:
1. Starting from the state after Scenario 1 (SPC off, single Y, no split),
   record the current page/console error baseline.
2. Read the current X-axis full column range: capture the data-derived
   minimum and maximum of `CAST Idea ID` as `fullMin` and `fullMax`.
3. In the property panel under **X axis**, set **X Min** to a value strictly
   greater than `fullMin` (e.g. `fullMin + (fullMax - fullMin) * 0.1`) and
   set **X Max** to a value strictly less than `fullMax`
   (e.g. `fullMax - (fullMax - fullMin) * 0.1`).
4. Zoom the chart in with the mouse wheel over the plot area. Confirm the
   zoom takes effect (the line chart's zoomed event fires).
5. Trigger **Reset View** (context menu > Reset View, or the `pcmdResetView`
   command). Confirm the reset-view event fires.
6. Assert: Reset View reset the zoom (the zoomed and reset-view events both
   fired) but **X Min and X Max are unchanged** — still the values set in
   Step 3, NOT restored to `fullMin`/`fullMax`. No new console/page error was
   raised since the Step 1 baseline.

Expected:
- The wheel zoom takes effect and Reset View resets it (both events fire).
- After Reset View, X Min and X Max still equal the Step 3 values — Reset View
  clears the interactive zoom, not the configured axis bounds.
- No new console or page error is raised during the zoom and reset.

## Automation notes

- The console/page error baseline is the `grok.shell.warnings` count plus the
  page-error count; the Scenario 1 responsiveness check is a
  `page.evaluate(() => true)` round-trip that must resolve within timeout.
- Scenario 2: Reset View is triggered via the context menu or the
  `pcmdResetView` command; the zoom and reset are confirmed through the
  viewer's zoomed / reset-view events.
