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
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      After enabling the regression line (GROK-19732 repro path), no console
      error or page error is raised (grok.shell.warnings delta == 0, page error
      count unchanged).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After enabling a rolling average overlay without a split applied, no
      console error is raised.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      After enabling a standard-deviation overlay without a split applied, no
      console error is raised.
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      After enabling a rolling average overlay with a split applied
      (per-category, GROK-20218 repro path), no console error is raised.
  - anchor: "Scenario 1 Step 11"
    expectation: >-
      After reverting overlays and split to none, no console error is raised
      (no-error teardown confirmed).
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      After adding a formula line and a formula band via the Formula Lines
      dialog, no console error is raised.
  - anchor: "Scenario 2 Step 7"
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
5. Record the current page-error count and `grok.shell.warnings` baseline before
   each action in the scenarios below (used as the no-error floor delta reference).

## Scenarios

### Scenario 1: Regression line, rolling average, and standard-deviation overlays

Steps:
1. In the Line Chart property panel, navigate to the **Regression** (or **Overlays**)
   section.
2. Enable the **Regression line** toggle.
3. Assert: no console error or page error since baseline (warnings delta == 0,
   page error count unchanged). The regression overlay is canvas-rendered — this is
   the no-error floor for the GROK-19732 repro path (Step 3 no-error floor).
4. In the property panel, enable the **Rolling average** (moving average) overlay
   (no split active at this point — overall dataset).
5. Assert: no console error since the previous checkpoint (Step 5 no-error floor).
6. Enable the **Standard deviation** overlay (overall dataset, no split).
7. Assert: no console error since the previous checkpoint (Step 7 no-error floor).
8. In the property panel under **Data**, set **Split** to `Stereo Category`
   (a categorical column in spgi-100). Keep rolling average and std-dev overlays active
   so they render per-category (GROK-20218 repro path).
9. Assert: no console error since the previous checkpoint. The rolling average and
   std-dev overlays render per category with no freeze (Step 9 no-error floor —
   GROK-20218 path).
10. Disable rolling average, std-dev, and regression line overlays. Set **Split**
    back to none.
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

Actuation note: the formula line and band are added by assigning the viewer's
`formulaLines` JSON property directly, not through the Formula Lines dialog. The
dialog is a modal whose expression editor and preview are not scriptable (its
interactive contents render on a canvas, matching
the Select-columns dialog pattern documented for the Line chart column editors),
and `formulaLines` is exactly the serialized form the dialog writes and the
layout round-trip persists. The observable — the formula-lines configuration
surviving a layout save → clear → reapply (GROK-19943 + GROK-19949) — is
exercised identically through the prop path. Adding the lines by hand through the
dialog stays a human-side variation of this scenario.

Steps:
1. Starting from a clean chart state (no overlays, no split, X = `CAST Idea ID`,
   Y = `Chemical Space X`), record the current warning/error baseline.
2. Add the formula line and band (in the real UI this is the **Formula Lines**
   dialog; the exact menu path is unverified — the grok-browser refdoc's Tools
   submenu does not list Formula Lines — so the spec drives the serialized
   `formulaLines` property directly, see the actuation note below).
3. Add one formula line: set **Expression** to a constant (e.g. `500`) and confirm.
   Add one formula band: set **Expression** to `400..600` or equivalent band syntax
   and confirm.
4. Assert: no console error or page error since baseline (Step 4 no-error floor).
   The formula line and band are canvas-rendered.
5. Note the number of formula-lines entries currently configured (expected: 2 — one
   line, one band).
6. Save the current viewer layout: use the viewer context menu > **Save Layout**
   (or the platform layout toolbar action), providing a test layout name
   (e.g. `test-overlays-layout`).
7. Close the Line Chart viewer (or reset its state), then reapply the saved layout
   (context menu > **Apply Layout** or select from the layout list). Assert the
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
