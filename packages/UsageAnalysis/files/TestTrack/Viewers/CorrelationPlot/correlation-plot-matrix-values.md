---
feature: correlationplot
realizes_atlas:
  - correlationplot.cp.matrix-values-scope-persist
  - correlationplot.int.numerical-columns-only
realizes:
  - viewers.correlation-plot
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-17480
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 2"
    expectation: >-
      Each off-diagonal cell value equals the runtime-computed Pearson
      correlation between that column pair over the full row set, within numeric
      tolerance (1e-3).
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      A diagonal cell has no numeric correlation value and renders a per-column
      histogram instead (cellType != correlation value cell).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      X Columns and Y Columns list editors show only numerical columns;
      non-numerical columns (e.g. SEX, RACE) are absent from both pickers.
  - anchor: "Scenario 2 Step 2"
    expectation: >-
      After switching Correlation Type to Spearman, each checked cell value
      equals the runtime-computed Spearman coefficient for that pair within
      numeric tolerance (1e-3), and differs from the Pearson value for
      non-monotone-linear pairs.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      After setting Show Pearson R to false, the matrix repaints (settle-gated
      canvas diff exceeds the pre-flip settle noise) while the backing value for
      a spot-checked pair stays readable and unchanged via getCorrelation.
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      After narrowing X Columns to [AGE, HEIGHT, WEIGHT] (in that order) and Y
      Columns to [AGE, HEIGHT], the matrix has exactly 3 X columns and 2 Y
      columns, cell values are correct for the new set and order, and column
      names remain visible (GROK-17480 regression guard).
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      After switching Row Source to Selected, the cell value for a spot-checked
      pair equals the Pearson coefficient recomputed over only the selected row
      indices, within tolerance.
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      After setting filter formula '${AGE} > 40', the spot-checked cell value
      equals the Pearson coefficient over the formula-masked rows, within
      tolerance, AND df.filter.trueCount is unchanged from the value before the
      formula was applied (the formula filter is viewer-local).
  - anchor: "Scenario 4 Step 3"
    expectation: >-
      After a larger Default Cell Font is applied, the matrix visibly repaints
      (settle-gated canvas diff exceeds the pre-change settle noise) with no
      errors raised and the backing value stable.
  - anchor: "Scenario 5 Step 5"
    expectation: >-
      After re-applying the saved layout, the viewer set equals the set recorded
      before editing the layout away, Correlation Type is Spearman, Show Pearson
      R is false, X Columns is [AGE, HEIGHT, WEIGHT], Y Columns is [AGE,
      HEIGHT], and the filter formula is '${AGE} > 40'. A spot-checked cell
      value matches the runtime Spearman reference over the formula mask within
      tolerance.
  - anchor: "Scenario 6 Step 4"
    expectation: >-
      After reopening the saved project, the Correlation plot is present,
      Correlation Type is Spearman, Show Pearson R is false, X Columns is [AGE,
      HEIGHT, WEIGHT], Y Columns is [AGE, HEIGHT], and a spot-checked cell value
      matches the runtime Spearman reference over the full row set within
      tolerance.
realized_as:
  - correlation-plot-matrix-values-spec.ts
---

# Correlation Plot — Matrix Values, Scope, and Persistence

## Setup

1. Open the demog dataset (System:DemoFiles/demog.csv).
2. Add a Correlation Plot viewer to the current Table View via
   Viewers Panel > Correlation Plot (or the Add Viewer toolbar).
3. Note the current number of rows passing the global filter as the baseline for later assertions.
4. Locate the Correlation Plot viewer in the view for subsequent inspection steps.

## Scenarios

### Scenario 1: Default state and numerical-column filter

Steps:
1. Observe the matrix immediately after adding the viewer — no settings have been changed.
2. For three off-diagonal pairs (e.g. AGE vs HEIGHT, AGE vs WEIGHT, HEIGHT vs WEIGHT),
   read each cell's displayed correlation value and verify it matches the expected Pearson
   coefficient for that column pair over the full demog row set, within tolerance 1e-3.
3. Confirm column labels are rendered vertically and the matrix covers the dataset's numerical
   columns (AGE, HEIGHT, WEIGHT, etc.) determined automatically.
4. Inspect a diagonal cell (e.g. AGE vs AGE) — confirm it has no numeric correlation value
   and instead shows a per-column histogram.
5. Right-click the viewer and open Columns > X Columns. Inspect the column list — confirm
   only numerical columns appear; non-numerical columns (SEX, RACE, SUBJ) are absent.
   Close the picker without changes.

Expected:
- Off-diagonal cell values equal the Pearson coefficient within 1e-3 (Step 2).
- Diagonal cell renders a histogram, not a correlation number (Step 4).
- X Columns picker shows only numerical columns; SEX and RACE are absent (Step 5).

### Scenario 2: Correlation type switch, Show Pearson R, and column narrowing

Steps:
1. On the Context Panel, open the Correlation Type dropdown and switch it to Spearman.
   The setting stays at Spearman for all subsequent steps in this scenario.
2. For the same three pairs as Scenario 1, read each cell's displayed value and verify it matches
   the Spearman coefficient for that column pair over the full row set, within tolerance 1e-3.
3. Confirm the values differ from the Pearson values obtained in Scenario 1 Step 2
   (at least one pair should differ, confirming the type switch is reflected in the displayed data,
   not just the label).
4. On the Context Panel, set Show Pearson R to false. The setting stays false.
   Verify the matrix repaints without the in-cell R digits and with narrower value
   columns — capture a settle-gated canvas snapshot before and after the flip and
   assert the two renders differ; verify the backing correlation value for a
   spot-checked pair is still readable via getCorrelation
   (the property affects display, not the underlying data).
5. On the Context Panel, open Columns > X Columns and narrow the selection to
   [AGE, HEIGHT, WEIGHT] in that order. Open Columns > Y Columns and narrow to [AGE, HEIGHT].
   The column sets stay at these values.
6. Confirm the matrix shows exactly 3 X-axis columns and 2 Y-axis columns; verify cell values
   for the now-present pairs match the Spearman coefficient over the full row set within tolerance;
   verify column names are visible in the header row
   (GROK-17480 regression guard — names must not be hidden after reorder or narrow).

Expected:
- Cell values equal the Spearman coefficient within 1e-3 after the switch (Step 2).
- The settle-gated canvas diff shows the in-cell R digits gone and the value columns narrowed
  after disabling Show Pearson R, while getCorrelation still returns the backing value (Step 4).
- Matrix is 3x2, values are correct, and column names are visible after narrowing (Step 6).

### Scenario 3: Row Source and viewer-local filter formula

Steps:
1. Select a band of approximately 20 rows in the grid (e.g. rows 0–19) using the row selector.
2. On the Context Panel > Data, set Row Source to Selected. The setting stays at Selected.
3. For the AGE vs HEIGHT pair, verify the displayed cell value matches the Pearson coefficient
   computed over only the selected rows within tolerance 1e-3.
4. On the Context Panel > Data, set Row Source to Filtered (the default). The displayed value
   for the same pair must recompute back toward the full-set reference.
5. On the Context Panel > Data, expand the Filter section and enter the formula
   `${AGE} > 40` in the filter expression field. Apply the formula.
   Verify the displayed AGE vs HEIGHT value matches the Pearson coefficient for rows where AGE > 40,
   within tolerance 1e-3.
   Confirm the global filter row count is unchanged from the baseline recorded in Setup
   (the formula is viewer-local — it must not alter the global filter).
6. Clear the filter formula field and confirm the displayed value returns toward the full-set
   Pearson reference.

Expected:
- After Row Source = Selected, the cell value equals Pearson over selected rows within tolerance (Step 3).
- After the filter formula, the cell value equals Pearson over formula-masked rows within tolerance,
  AND the global filter row count is unchanged (Step 5).

### Scenario 4: defaultCellFont structural signal

Steps:
1. Capture a settle-gated canvas snapshot of the matrix as the baseline render.
2. On the Context Panel > Style, set Default Cell Font to "bold 16px" (or a size clearly
   larger than the default).
3. Verify the larger cell font visibly repaints the matrix (taller rows / larger digits) —
   capture a settle-gated canvas snapshot before and after the font change and assert the
   two renders differ, with no errors raised.
4. Restore Default Cell Font to the original value and confirm the render returns to the
   baseline snapshot (settle-gated diff is approximately 0) with no errors.

Expected:
- After a larger Default Cell Font, the settle-gated canvas diff confirms the matrix
  repaints (taller rows / larger digits) with no errors raised (Step 3).

### Scenario 5: Layout persistence — save, re-arm, re-apply

Steps:
1. With all settings from Scenarios 2 and 3 active (Spearman, Show Pearson R = false,
   X = [AGE, HEIGHT, WEIGHT], Y = [AGE, HEIGHT], filter formula cleared, Row Source = Filtered),
   save the current view layout via the Save Layout button in the View ribbon
   (record the saved layout name, e.g. "CP-matrix-persist-test").
2. Record the current viewer set (Correlation Plot present, no Scatter Plot).
3. Close the Correlation Plot viewer and add a Scatter Plot viewer to the view
   (this makes the current state differ from the saved layout).
4. Re-apply the saved layout via the Apply Layout option in the View ribbon.
5. Read the current viewer set, Correlation Type, Show Pearson R, X Columns, Y Columns,
   and filter formula from the restored Correlation Plot.
   Spot-check the AGE vs HEIGHT cell value against the Spearman coefficient over the full row set,
   within tolerance.
6. Delete the probe layout via the Layouts panel (cleanup step).

Expected:
- After re-applying the layout, the viewer set equals the saved set (Correlation Plot, no Scatter Plot);
  Correlation Type is Spearman; Show Pearson R is false; X Columns is [AGE, HEIGHT, WEIGHT];
  Y Columns is [AGE, HEIGHT]; filter formula is empty; spot-check value matches Spearman reference (Step 5).

### Scenario 6: Project persistence — save, close all, reopen

Steps:
1. With the settings from Scenario 2 active on the Correlation Plot (Spearman, Show Pearson R = false,
   X = [AGE, HEIGHT, WEIGHT], Y = [AGE, HEIGHT]), save a project named
   "CP-matrix-persist-project" using the ribbon Save button.
2. Use Close All to close all open tables and viewers.
3. Reopen the saved project via File > Recent Projects or the Projects panel.
4. Locate the Correlation Plot in the reopened view.
5. Read Correlation Type, Show Pearson R, X Columns, and Y Columns from the restored viewer.
   Spot-check the AGE vs HEIGHT cell value against the Spearman coefficient over the full row set,
   within tolerance.
6. Delete the probe project via the Projects panel (cleanup step).

Expected:
- After reopening the project, Correlation Plot is present; Correlation Type is Spearman;
  Show Pearson R is false; X Columns is [AGE, HEIGHT, WEIGHT]; Y Columns is [AGE, HEIGHT];
  spot-check value matches Spearman reference within tolerance (Step 5).

## Automation notes

- target_layer rationale: all driving is multi-step UI gestures on the Context Panel, Columns
  pickers, row selector, filter formula field, layout save/apply, and the project round-trip;
  value reads use the exposed getCorrelation(c1, c2) via JS API
  within the Playwright page context. The Playwright layer is the only one that can drive
  both the UI gestures and the in-page JS API reads in the same test run.
- Setup actuation: open demog.csv via readDataframe helper and add to workspace via
  grok.shell.addTableView; obtain the viewer object via findViewer; store a reference to
  correlation(col1, col2) callable for value spot-checks; record df.filter.trueCount as
  the baseline filter count.
- WASM tolerance: Pearson value assertions compare to a runtime-computed Stats.corr / Stats.pearson
  reference, never to a hardcoded number, because the Pearson branch may route through GrokML WASM.
  Tolerance is 1e-3 (conservative); tighten to 1e-6 if live recon shows stable agreement.
- Viewer-local filter invariant (Scenario 3 Step 5): the formula filter is viewer-local
  (props.filter feeds combinedFilter, not df.filter). The test MUST assert df.filter.trueCount
  is unchanged from the baseline; if it changes, the viewer has leaked the formula to the global
  filter — that is the anti-trap.
- Show Pearson R actuation (Scenario 2 Step 4): actuate via cp.props.showPearsonR = false;
  signal is a settle-gated canvas snapshot before/after the flip
  (snapshotCanvasColors/diffCanvasColors from ../../helpers/viewers.ts) asserting the renders
  differ, plus a getCorrelation spot-check confirming the backing value stays readable with the
  display off. Do not read gridCol.showValue or the column width — the inner ColumnGrid is not
  JS-exposed; the showValue/width flips are the Dart-side cause only.
- defaultCellFont actuation (Scenario 4 Step 3): actuate via cp.props.defaultCellFont; signal is
  a settle-gated canvas snapshot before/after the font change
  (snapshotCanvasColors/diffCanvasColors from ../../helpers/viewers.ts) asserting the renders
  differ, plus the no-error floor. Do not read grid rowHeight — the inner ColumnGrid is not
  JS-exposed; rowHeight = parseSize(fontValue) x 1.4 is the documented Dart-side cause only.
- Global filter baseline: record df.filter.trueCount in Setup Step 3; compare in Scenario 3
  Step 5 to confirm viewer-local filter did not alter the global filter.
- Diagonal cell assertion: verified via the backing-data probe (the xCol == yCol skip in
  _refreshValues means no FloatColumn entry exists for that pair); a DOM or canvas check
  is secondary confirmation.
- Row Source default: the actual default is 'Filtered' per look_and_feel.dart:208 and runtime
  confirmation — assertions in Scenario 3 rely on this.
- GROK-17480 guard: the column-name-visible assertion in Scenario 2 Step 6 is the specific
  regression guard for the "column names hidden after reorder" bug (status: fixed).
- Interaction realized: correlationplot.int.numerical-columns-only is asserted in Scenario 1
  Step 5 (X Columns picker shows only numerical columns).
- Finally teardown: probe layout (Scenario 5 Step 6) and probe project (Scenario 6 Step 6)
  must be deleted in finally blocks even on failure to leave the environment clean.
  Scenario 5 teardown: delete via Layouts panel or JS API (grok.dapi.layouts).
  Scenario 6 teardown: delete via Projects panel or JS API (grok.dapi.projects.delete(...)).
- Column rename prohibition: no scenario step renames a column on the viewer — that operation
  is operator-banned per the atlas edge_cases.
