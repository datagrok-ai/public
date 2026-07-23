---
feature: linechart
realizes:
  - viewers.line-chart
realizes_atlas:
  - linechart.cp.analytical-overlays
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19732
    status: fixed
  - id: GROK-20218
    status: fixed
  - id: GROK-19943
    status: fixed
  - id: GROK-19949
    status: fixed
realized_as:
  - analytical-overlays-spec.ts
expected_results:
  - anchor: "S1 steps 1-3: enable regression line, no-error floor"
    expectation: >-
      After enabling the regression line (GROK-19732 repro path), no console
      error or page error is raised (grok.shell.warnings delta == 0, page error
      count unchanged).
  - anchor: "S1 steps 4-5: enable rolling average overlay, no split"
    expectation: >-
      After enabling a rolling average overlay without a split applied, no
      console error is raised.
  - anchor: "S1 steps 6-7: enable standard-deviation overlay, no split"
    expectation: >-
      After enabling a standard-deviation overlay without a split applied, no
      console error is raised.
  - anchor: "S1 steps 8-9: apply split, overlays render per category"
    expectation: >-
      After enabling a rolling average overlay with a split applied
      (per-category, GROK-20218 repro path), no console error is raised.
  - anchor: "S1 steps 10-11: revert overlays and split"
    expectation: >-
      After reverting overlays and split to none, no console error is raised
      (no-error teardown confirmed).
  - anchor: "S2 steps 1-5: add a formula line and a formula band"
    expectation: >-
      After adding a formula line and a formula band via the Formula Lines
      dialog, no console error is raised.
  - anchor: "S2 steps 6-8: save layout, clear, reapply — formula lines restored"
    expectation: >-
      After saving a layout with formula lines and reapplying that layout
      (GROK-19943 + GROK-19949 round-trip), the formula-lines configuration is
      still present: the formula lines entry count read back from the viewer
      properties matches the count set before saving, and no console error is
      raised.
---

# Line Chart — Analytical Overlays

## Setup

1. Close all open views.
2. Open the spgi-100 dataset (use `readDataframe` helper or open via the file browser).
3. Add a Line Chart viewer to the active table view (use `findViewer` helper or
   add via the Toolbox Viewers icon).
4. In the Line Chart property panel, set **X column** to `CAST Idea ID` and
   **Y columns** to `Chemical Space X`.
5. Record the current console/page error baseline before each action in the
   scenarios below (used as the no-error floor delta reference).

## Scenarios

### Scenario 1: Regression line, rolling average, and standard-deviation overlays

Steps:
1. In the Line Chart property panel, navigate to the **Lines** section (the
   regression and moving-average toggles live there).
2. Enable the **Show Regression Line** checkbox.
3. Assert: no console error or page error since baseline (warnings delta == 0,
   page error count unchanged). The regression overlay is canvas-rendered — this is
   the no-error floor for the GROK-19732 repro path (Step 3 no-error floor).
4. In the property panel > Lines, enable the rolling-average overlay: check
   **Show Moving Average Line** (no split active at this point — overall dataset).
5. Assert: no console error since the previous checkpoint (Step 5 no-error floor).
6. Enable the standard-deviation overlay: check **Show Moving Average Deviation**
   (property panel > Lines — shades a ±1 standard-deviation band around the
   moving-average line; overall dataset, no split).
7. Assert: no console error since the previous checkpoint (Step 7 no-error floor).
8. In the property panel under **Data**, set **Split** to `Stereo Category`
   (a categorical column in spgi-100). Keep rolling average and std-dev overlays active
   so they render per-category (GROK-20218 repro path).
9. Assert: no console error since the previous checkpoint. The rolling average and
   std-dev overlays render per category with no freeze (Step 9 no-error floor —
   GROK-20218 path).
10. Disable the overlays — uncheck **Show Moving Average Line**, **Show Moving
    Average Deviation**, and **Show Regression Line** (property panel > Lines).
    Set **Split** back to none.
11. Assert: no console error since the previous checkpoint (Step 11 no-error teardown).

Expected:
- Enabling the regression line raises no console error (Step 3 — GROK-19732 path).
- Enabling rolling average overlay (no split) raises no console error (Step 5).
- Enabling standard-deviation overlay (no split) raises no console error (Step 7).
- With a split applied, rolling average and std-dev overlays per category raise no
  console error (Step 9 — GROK-20218 path).
- Reverting all overlays and split raises no console error (Step 11 — no-error
  teardown).

### Scenario 2: Formula lines layout round-trip

Steps:
1. Starting from a clean chart state (no overlays, no split, X = `CAST Idea ID`,
   Y = `Chemical Space X`), record the current warning/error baseline.
2. Open the **Formula Lines** dialog: right-click the chart > **Tools** >
   **Formula Lines...**
3. Add one formula line: set **Expression** to a constant (e.g. `500`) and confirm.
   Add one formula band: set **Expression** to `400..600` or equivalent band syntax
   and confirm.
4. Assert: no console error or page error since baseline (Step 4 no-error floor).
   The formula line and band are canvas-rendered.
5. Note the number of formula-lines entries currently configured (expected: 2 — one
   line, one band).
6. Save the current view layout: top menu **View > Layout > Save to Gallery**.
7. Close the Line Chart viewer (or reset its state), then reapply the saved layout
   (**View > Layout > Open Gallery**, click the saved layout). Assert the
   following after reapply (GROK-19943 + GROK-19949 round-trip, Step 7):
   - The formula-lines count re-read from the viewer equals 2 (the same count set
     before saving — the configuration survived serialize→reapply).
   - No console error or page error is raised.
8. Clean up: remove the formula lines entries and delete the test layout if the
   API supports it.

Expected:
- Adding a formula line and band via the Formula Lines dialog raises no console
  error (Step 4).
- After saving and reapplying the layout, the formula-lines configuration is
  restored: entry count == 2 and no console error is raised (Step 7 — GROK-19943
  + GROK-19949 round-trip).

## Automation notes

Setup: the console/page error baseline is recorded as the page-error count plus
the `grok.shell.warnings` length; every "no console error" assert compares
against this delta.

Scenario 2: the formula line and band are added by assigning the viewer's
`formulaLines` JSON property directly, not through the Formula Lines dialog. The
dialog is a modal whose expression editor and preview are not scriptable (its
interactive contents render on a canvas, matching
the Select-columns dialog pattern documented for the Line chart column editors),
and `formulaLines` is exactly the serialized form the dialog writes and the
layout round-trip persists. The observable — the formula-lines configuration
surviving a layout save → clear → reapply (GROK-19943 + GROK-19949) — is
exercised identically through the prop path. Adding the lines by hand through the
dialog stays a human-side variation of this scenario.
