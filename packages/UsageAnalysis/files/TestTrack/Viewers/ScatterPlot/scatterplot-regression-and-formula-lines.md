---
feature: scatterplot
realizes_atlas:
  - scatterplot.cp.regression-and-formula-lines
  - scatterplot-color-drives-per-category-regression
realizes:
  - viewers.scatter-plot
priority: p2
target_layer: playwright
realized_as:
  - scatterplot-regression-and-formula-lines-spec.ts
coverage_type: regression
related_bugs:
  - id: GROK-17282
    status: fixed
  - id: GROK-17116
    status: fixed
  - id: GROK-18381
    status: fixed
  - id: GROK-20458
    status: fixed
  - id: GROK-16214
    status: fixed
  - id: github-2530
    status: fixed
expected_results:
  - anchor: "Regression line on a logarithmic Y axis"
    expectation: "With the Y axis set to logarithmic, turning Show Regression Line
      on changes the overlay rendering against the regression-off baseline
      measured on the same log axis — a line is drawn — and no new console or
      page error appears (GROK-17116 floor); the equation and statistics box is
      painted on the overlay and exposes no readable text, so the coefficients
      are deliberately not asserted"
  - anchor: "Regression line on a logarithmic Y axis"
    expectation: "Turning Show Regression Line off and the Y Axis Type back to
      linear restores the overlay rendering toward the starting baseline"
  - anchor: "Regression Per Category and the axis time unit with the regression line
      on"
    expectation: "With Show Regression Line on and a categorical Color column set,
      turning Regression Per Category on and then off raises no new console or
      page error and the viewer keeps rendering (a no-error floor: the row's
      greying is not assertable in this flow because Show Regression Line is
      already on, see the automation notes)"
  - anchor: "Regression Per Category and the axis time unit with the regression line
      on"
    expectation: "With the regression line on and a datetime column on X, choosing a
      time unit in the on-viewer X selector raises no new console or page error
      and the viewer keeps rendering (GROK-18381)"
  - anchor: "Regression Per Category and the axis time unit with the regression line
      on"
    expectation: "Reverting the time unit, the axis column, Regression Per Category,
      Show Regression Line and Color returns the viewer to the setup
      configuration"
  - anchor: "Formula lines dialog — add a line, edit it across the equals sign, reopen"
    expectation: "A formula whose left side is not a single column is REJECTED by
      the dialog — the editor takes the invalid-input state and OK is disabled —
      while the form stays alive; correcting it back to the supported ${column}
      = expression contract clears the flag and re-enables OK, and after OK plus
      reopen the line is present, loaded into the formula editor, and NOT in the
      invalid-input state (GROK-17282: no dead state)"
  - anchor: "Formula lines dialog — add a line, edit it across the equals sign, reopen"
    expectation: "The viewer's stored formula-line configuration holds exactly one
      line whose formula is the edited one"
  - anchor: "Formula line across an axis-column change, and a band on a logarithmic
      axis"
    expectation: "Changing both axis columns leaves the viewer's stored formula-line
      configuration identical to the value recorded before the change
      (GROK-16214)"
  - anchor: "Formula line across an axis-column change, and a band on a logarithmic
      axis"
    expectation: "With a band configured, switching the Y Axis Type to logarithmic
      and back raises no new console or page error (GROK-20458)"
  - anchor: "Formula line across an axis-column change, and a band on a logarithmic
      axis"
    expectation: "Deleting every formula line and restoring the axis columns returns
      the viewer to the setup configuration with an empty formula-line
      configuration"
  - anchor: "Hover sweep with a formula line present"
    expectation: "Sweeping the pointer across several known points with a formula
      line present shows the tooltip and raises no new console or page error
      (github-2530)"
  - anchor: "Hover sweep with a formula line present"
    expectation: "The probe formula line is deleted afterwards and the formula-line
      configuration is empty again"
  - anchor: "Moving average line, window, per category and deviation"
    expectation: "Each step of the ladder — Show Moving Average Line on, a clearly
      wider Moving Average Window, Moving Average Per Category on with a
      categorical Color column, Show Moving Average Deviation on — changes the
      DATA canvas rendering, with the deviation step the largest change of the
      ladder"
  - anchor: "Moving average line, window, per category and deviation"
    expectation: "The overlay canvas rendering is unchanged across the whole
      moving-average ladder — the moving average draws on the data canvas,
      unlike the regression line, which is what separates the two surfaces"
  - anchor: "Moving average line, window, per category and deviation"
    expectation: "Turning all four moving-average settings back off and restoring
      the window restores the data canvas rendering to the baseline, and
      clearing Color returns the viewer to the setup configuration"
---

# Scatter Plot — Regression Line, Formula Lines, Moving Average

## Purpose

Verifies the Scatter plot's line-drawing surface on the demog dataset: the
regression line on a logarithmic axis, the per-category regression toggle, the
axis time-unit selector while the regression line is on, the Formula Lines
dialog (adding a line, editing it across the equals sign, and its survival
across an axis-column change), a band under a logarithmic axis, hovering with
lines present, and the moving-average feature.

Both the regression line and its equation/statistics box are painted on the
viewer's overlay canvas and expose no readable text or accessor, so the
regression guards are honest floors — a rendering change proving a line is
drawn, plus a clean-console guard — and never a claim about the equation. The
moving average, by contrast, draws on the data canvas, and that separation is
itself asserted.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Scatter plot to the current table view via the Toolbox viewer icon.
4. Using the on-viewer column selectors, set X to WEIGHT and Y to HEIGHT (the
   auto-pick already matches on demog, but the auto-selection order is not
   contractual).
5. Confirm the starting state: **Show Regression Line** off, no Color column,
   no formula lines configured.

## Scenarios

### Scenario 1: Regression line on a logarithmic Y axis

Steps:
1. Open the viewer settings and, in the **Y** section, set **Y Axis Type** to
   logarithmic.
2. Note the current browser console-error and page-error counts, and measure
   the plot's overlay rendering — this is the regression-off baseline on the
   logarithmic axis.
3. In the **Lines** section, turn **Show Regression Line** on.
4. Verify the overlay rendering changed against the baseline — a regression
   line is drawn on the logarithmic axis (GROK-17116 guard).
5. Verify no new console or page error appeared.
6. Turn **Show Regression Line** off and set **Y Axis Type** back to linear.
7. Verify the overlay rendering moved back toward the starting baseline.

Expected:
- With the logarithmic Y axis, turning Show Regression Line on changes the overlay rendering and raises no new console or page error; the equation text has no readable signal and is not asserted
- Turning the regression line off and the axis back to linear restores the overlay rendering toward the baseline

### Scenario 2: Regression Per Category and the axis time unit with the regression line on

Steps:
1. Click the on-viewer **Color** selector and set Color to RACE.
2. In the viewer settings **Lines** section, turn **Show Regression Line** on.
3. Note the current console-error and page-error counts.
4. Turn **Regression Per Category** on; wait for the redraw.
5. Verify no new console or page error appeared and the viewer keeps rendering.
6. Turn **Regression Per Category** off; wait for the redraw and verify the
   same no-error floor holds.
7. Click the on-viewer **X** selector and set X to STARTED (a datetime column)
   — the time-unit selector next to the X axis becomes available.
8. Choose **year** in that time-unit selector.
9. Verify no new console or page error appeared and the viewer keeps rendering
   (GROK-18381 guard).
10. Revert: clear the time unit, set X back to WEIGHT, turn **Show Regression
    Line** off, and clear the **Color** column.

Expected:
- Toggling Regression Per Category on and off with the regression line on raises no new console or page error and the viewer keeps rendering
- Choosing a time unit on a datetime X axis with the regression line on raises no new console or page error and the viewer keeps rendering
- The time unit, axis column, Regression Per Category, Show Regression Line and Color are all returned to the setup configuration

### Scenario 3: Formula lines dialog — add a line, edit it across the equals sign, reopen

The Formula Lines dialog embeds its own preview scatter plot, so every step
below acts on the controls inside the dialog; readings of the main plot happen
only after the dialog is closed.

Steps:
1. Right-click the plot area and choose **Tools > Formula Lines...** — the
   Formula Lines dialog opens on its **Viewer** tab.
2. Click **ADD NEW** and choose **Line** — the formula editor is pre-filled
   with an equation over the two axis columns.
3. Edit the formula so the second column moves across the equals sign — turn
   `Y = X + 1` into `Y - X = 1`. This form is not supported: the left side must be
   a single column. Verify the editor takes the invalid-input state, that **OK** is
   disabled, and that the dialog stays open and editable rather than dying.
4. Correct the formula to the supported form `Y = X + 1`, verify the invalid-input
   state clears and **OK** becomes available again, then click **OK** to commit and
   close the dialog.
5. Reopen **Tools > Formula Lines...**.
6. Verify the edited line is present and loaded into the formula editor, and
   that the editor is NOT in the invalid-input state (GROK-17282 guard).
7. Verify the viewer's stored formula-line configuration holds exactly one
   line, carrying the edited formula.
8. Close the dialog with **OK**, leaving the line in place for the next
   scenario.

Expected:
- After OK and reopen, the formula line is present, loaded into the editor, and not flagged as invalid
- The viewer's formula-line configuration holds exactly one line whose formula is the edited one

### Scenario 4: Formula line across an axis-column change, and a band on a logarithmic axis

Steps:
1. Record the viewer's stored formula-line configuration (the line left by
   Scenario 3).
2. Using the on-viewer selectors, change X to AGE and Y to WEIGHT.
3. Verify the stored formula-line configuration is identical to the recorded
   value (GROK-16214 guard).
4. Open **Tools > Formula Lines...**, click **ADD NEW** and choose
   **Band > Horizontal**, then click **OK**.
5. Note the current console-error and page-error counts.
6. In the viewer settings **Y** section, set **Y Axis Type** to logarithmic;
   wait for the redraw, then set it back to linear.
7. Verify no new console or page error appeared across the switch and back
   (GROK-20458 guard).
8. Revert: open **Tools > Formula Lines...**, delete every line with the
   **Delete** button, click **OK**, and restore X to WEIGHT and Y to HEIGHT.

Expected:
- Changing both axis columns leaves the stored formula-line configuration identical
- Switching the Y axis to logarithmic with a band configured, and back, raises no new console or page error
- All formula lines are deleted and the axis columns are restored

### Scenario 5: Hover sweep with a formula line present

Steps:
1. Open **Tools > Formula Lines...**, click **ADD NEW**, choose **Line**, and
   click **OK** — one plain formula line is configured.
2. Note the current console-error and page-error counts.
3. Sweep the pointer across the plot over several known data points, pausing
   on each long enough for the tooltip to appear.
4. Verify the tooltip appeared during the sweep and no new console or page
   error was raised (github-2530 guard).
5. Revert: open **Tools > Formula Lines...**, delete the line, and click
   **OK** — the formula-line configuration is empty again.

Expected:
- The hover sweep with a formula line present shows the tooltip and raises no new console or page error
- The probe formula line is deleted and the formula-line configuration is empty

### Scenario 6: Moving average line, window, per category and deviation

The moving average renders on the data canvas, whereas the regression line
renders on the overlay. This scenario asserts both the ladder of rendering
changes and that separation.

Steps:
1. Click the on-viewer **Color** selector and set Color to RACE. Confirm
   **Show Regression Line** is off.
2. Measure both the data canvas and the overlay rendering — the baseline.
3. In the viewer settings **Lines** section, turn **Show Moving Average Line**
   on; wait for the redraw and verify the data canvas rendering changed.
4. Set **Moving Average Window** to a clearly larger value than its default
   and verify the data canvas rendering changed again.
5. Turn **Moving Average Per Category** on — one line per color category —
   and verify the data canvas rendering changed again.
6. Turn **Show Moving Average Deviation** on and verify the data canvas
   rendering changed again, by the largest step of the ladder.
7. Verify the overlay rendering stayed unchanged across steps 3-6.
8. Revert: turn **Show Moving Average Deviation**, **Moving Average Per
   Category** and **Show Moving Average Line** off and restore **Moving
   Average Window** to its default.
9. Verify the data canvas rendering returned to the baseline, then clear the
   **Color** column.

Expected:
- Each moving-average step changes the data canvas rendering, the deviation step being the largest change of the ladder
- The overlay canvas rendering is unchanged across the whole ladder — the moving average draws on the data canvas, unlike the regression line
- Turning every moving-average setting back off restores the data canvas rendering to the baseline, and clearing Color returns the viewer to the setup configuration

## Automation notes

- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Scatter plot')`; the viewer is
  added via the Toolbox icon `[name="icon-scatter-plot"]` (synthetic `.click()`
  works). Resolve the viewer root as the `[name="viewer-Scatter-plot"]` element
  that is NOT inside a `.d4-dialog` — the Formula Lines dialog embeds its own
  preview scatter plot with its own `canvas`, `overlay`, `div-column-combobox-*`
  and sliders, so a bare `[name="viewer-Scatter-plot"]` query is ambiguous while
  that dialog is open. On-viewer column selectors are lowercase
  (`div-column-combobox-x|y|color|size`); they open on synthetic `mousedown`
  (never on synthetic `.click()`), then real typing plus `Enter` commits.
- Two canvases matter and must not be confused: `canvas[name="canvas"]` (data —
  markers, axes, grid lines, moving-average line and deviation band) and
  `canvas[name="overlay"]` (regression line, equation/statistics box, labels,
  drop lines). Both readings are settle-gated non-white pixel fractions; all
  thresholds are calibrated live at spec-authoring time and never carried over
  as settled constants, since they depend on the viewer size and the device
  pixel ratio. Sampling `getImageData` triggers the browser's
  `willReadFrequently` console advisory — whitelist it in the console guards.
- Scenario 1 (GROK-17116): the regression equation and statistics box has no
  DOM text and the viewer exposes no equation/slope/fit accessor, so the
  original "equation is not NaN + NaN*X" assertion is impossible. The honest
  floor is the overlay pixel delta plus a clean console; coefficient
  correctness stays an uncovered gap. Regression settings are the
  `[name="prop-show-regression-line"]` checkbox row and the mirrored menu leaf
  `[name="div-Tools---Show-Regression-Line"]`; the axis type is
  `[name="prop-y-axis-type"]` (choice editor) or
  `[name="div-Properties...---Y---Y-Axis-Type---Logarithmic"]`.
- Scenario 2 carries only a no-error floor over the **Regression Per Category**
  toggle. The row's disabled state IS readable — inline `opacity: 0.5` on
  `[name="prop-regression-per-category"]`, the same regime as every other row —
  but it is gated solely by Show Regression Line, which this flow switches on at
  step 2, before the toggle. The row therefore sits at opacity 1 for the whole of
  steps 4-6 and no greying transition exists here to assert. The drawn
  per-category split itself is an overlay-canvas outcome with no countable signal
  of its own, so the toggle is exercised and guarded against errors rather than
  graded on a rendering claim.
- Scenario 2, time unit: the on-viewer selector is
  `[name="input-aggr-selector-x-map"]`, `display:none` for numeric and
  categorical columns and `inline-block` for a datetime column. It requires a
  REAL option selection — setting `.value` plus a synthetic `change` does not
  commit; the JS fallback `sp.props.xMap = 'year'` commits and re-renders
  cleanly with the regression line on.
- Scenarios 3-5, Formula Lines: open via a synthetic `contextmenu` on
  `canvas[name="canvas"]` then `[name="div-Tools---Formula-Lines..."]` (the leaf
  is clickable synthetically even while its submenu is unexpanded). The dialog
  is `[name="dialog-Formula-Lines"]` with `[name="button-Add-new"]`,
  `[name="button-Delete"]`, `[name="button-OK"]`, `[name="button-CANCEL"]` and
  the ADD NEW leaves `div-Line`, `div-Line---Vertical`, `div-Line---Horizontal`,
  `div-Band---Vertical`, `div-Band---Horizontal`, `div-Region---*`. The formula
  editor is a `textarea.ui-input-editor` with NO `name=` — locate it as the
  textarea whose value matches `/\$\{/`. The invalid-input state is the class
  `d4-forced-invalid` on that textarea. The list of existing lines is a
  canvas-rendered grid inside the dialog, so line titles are not DOM text; the
  product-state signal is `sp.props.formulaLines`, a JSON string that is parsed
  for the round-trip and the identity assert.
- Scenario 5 aims the hover with `sp.worldToScreen(xValue, yValue)` plus the
  canvas bounding-rect origin (the helper returns canvas-local coordinates);
  the tooltip is the document-level `.d4-tooltip` and is populated by a
  synthetic `mousemove`.
- Scenario 6 drives `[name="prop-show-moving-average-line"]`,
  `[name="prop-moving-average-window"]`,
  `[name="prop-moving-average-per-category"]`,
  `[name="prop-show-moving-average-deviation"]` (or the mirrored
  `div-Properties...---Lines---*` leaves). The deviation band is the largest
  data-canvas delta of the ladder and the round-trip back to the baseline is
  exact, so the restore may be asserted as equality with a small tolerance; the
  overlay reading must be captured with the SAME sampling stride throughout so
  the "overlay unchanged" claim is meaningful.
- All console guards subtract a baseline and whitelist infra noise unrelated to
  the viewer (the shared dev server's WebSocket reconnect errors, the Claude
  runtime container timeout, the WebGPU `powerPreference` advisory and the
  `willReadFrequently` advisory caused by the pixel sampling itself).

---
{
  "order": 6,
  "datasets": ["System:DemoFiles/demog.csv"]
}
