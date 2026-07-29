---
feature: scatterplot
realizes_atlas:
  - scatterplot.cp.axes-and-encode
realizes:
  - viewers.scatter-plot
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-18945
    status: fixed
  - id: GROK-19334
    status: fixed
  - id: GROK-13110
    status: fixed
  - id: GROK-20395
    status: fixed
  - id: GROK-18411
    status: fixed
realized_as:
  - scatterplot-axes-and-encode-spec.ts
expected_results:
  - anchor: "Set the axes and the encodings through the on-viewer selectors"
    expectation: "Each pick made through an on-viewer selector is reflected by the
      viewer's own column settings (xColumnName, yColumnName, colorColumnName,
      sizeColumnName) read back AFTER a later unrelated action, and by the
      matching selector's on-viewer label text; no new browser console or page
      errors appear across the whole sequence"
  - anchor: "Set the axes and the encodings through the on-viewer selectors"
    expectation: "Clicking the Color selector LABEL TEXT (not the triangle) opens
      the column popup — the column-selector popup is present in the DOM
      (GROK-18411 guard)"
  - anchor: "One column on both axes, then renamed"
    expectation: "After renaming the column that serves as both X and Y through the
      grid column header, the viewer's xColumnName and yColumnName both read the
      NEW column name (GROK-19334 guard) and no new console or page error
      appears, including with a formula line configured on the viewer"
  - anchor: "Logarithmic and inverted axis with a reversed range window"
    expectation: "Setting X Axis Type to logarithmic, turning Invert X Axis on, and
      then entering X Min 60 with X Max 20 (minimum above maximum) produces NO
      'Wrong range' console error and no page error (GROK-13110 guard); the
      viewer stays attached and its viewport settles to a finite rectangle"
  - anchor: "Axis type control disabled for a datetime axis"
    expectation: "With a datetime column on X the property panel's X Axis Type row
      is shown disabled (dimmed) while the X time-unit row is active; with a
      numeric column on X the pair is reversed — the X Axis Type row is active
      and the time-unit row is dimmed (GROK-20395 guard, asserted on the
      property-panel rows only)"
  - anchor: "Layout and project persistence at peak configuration"
    expectation: "After saving the view layout, perturbing the view (adding another
      viewer and clearing Color) and re-applying the saved layout, the view's
      viewer set equals the SAVED set (a Scatter plot is present AND the
      later-added viewer is absent) AND the restored Scatter plot carries the
      recorded colorColumnName, xColumnName, yColumnName, sizeColumnName and
      markersColumnName (GROK-18945 guard); the probe layout is deleted
      afterwards even on failure"
  - anchor: "Layout and project persistence at peak configuration"
    expectation: "After saving the view as a project through the real ribbon Save
      button, Close All, and reopening the saved project, a Scatter plot viewer
      is present in the reopened view AND its axis and encoding settings equal
      the recorded values — a cross-session round-trip; the probe project is
      deleted afterwards even on failure"
---

# Scatter Plot — Axes, Encodings, Persistence

## Purpose

Verifies the Scatter plot's daily-run configuration surface on the demog
dataset: choosing the X and Y columns and the Color, Size and Marker encodings
through the real on-viewer selectors, the axis-type / inversion / range-window
controls including the reversed range window and the datetime type guard, the
rename of a column that serves both axes, and the survival of the fully
configured viewer across a saved layout and a saved project. Markers and axes
are painted on canvas, so every check reads a product-state signal instead:
the viewer's own column settings read back after an independent action, the
selector label text, the property-panel row state, and an error-free console.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table
   view to load.
3. Add a Scatter plot to the table view via the Toolbox viewer icon and wait
   for it to attach and render.
4. Record the current browser console-error and page-error counts as the
   baseline for the clean-console checks below.

## Scenarios

### Scenario 1: Set the axes and the encodings through the on-viewer selectors

The viewer auto-picks two axis columns when it attaches; the auto-pick order is
not contractual, so every column used later is set explicitly here.

Steps:
1. Click the on-viewer **X column selector** (bottom center of the plot) — the
   column popup opens; type AGE and press Enter to commit.
2. Click the on-viewer **Y column selector** (vertical, left edge) and pick
   HEIGHT the same way.
3. Click the on-viewer **Color selector** (vertical, right edge) directly on
   its **label text**, not on the dropdown triangle — verify the column popup
   opens (GROK-18411 guard: the label text itself must be a hit target). Pick
   RACE.
4. Click the on-viewer **Size selector** (top center) and pick WEIGHT.
5. Open the viewer settings (gear icon on the Scatter plot title bar) and, in
   the **Marker** section, set **Marker Column** to SEX.
6. Switch the X column to WEIGHT through the on-viewer X selector, then switch
   it back to AGE — an independent action between the picks and the read-back.
7. Read the viewer's column settings back: X is AGE, Y is HEIGHT, Color is
   RACE, Size is WEIGHT, Marker Column is SEX; the X, Y, Color and Size
   selector labels on the viewer show the same names.
8. Verify no new console or page error appeared across the sequence.

This configuration is the peak state the persistence scenario later saves, so
it is deliberately left in place; the scenarios in between revert only their
own changes.

Expected:
- Each pick through an on-viewer selector is reflected by the viewer's own column setting read back after a later unrelated action and by the selector's on-viewer label
- Clicking the Color selector label text (not the triangle) opens the column popup (GROK-18411)
- No new console or page errors across the whole sequence

### Scenario 2: One column on both axes, then renamed

A column may serve as both X and Y. Renaming it must propagate to both axes,
and it must do so even when the viewer carries a formula line.

Steps:
1. Note the current console-error and page-error counts.
2. Set both the X and the Y column to AGE through the on-viewer selectors.
3. Right-click the plot area, open **Tools > Formula Lines...**, click
   **ADD NEW** and choose **Line**, then click **OK** to close the dialog with
   one line configured.
4. In the grid, open the AGE column header context menu, choose
   **Column Properties...**, type a distinct probe name into the dialog's
   **New name** field and confirm with **OK**.
5. Verify the viewer's X column and Y column settings BOTH read the new name.
6. Verify no new console or page error appeared across the rename (the rename
   used to throw with a formula line present).
7. Revert: rename the column back to AGE, delete the formula line through the
   **Formula Lines...** dialog, and set the Y column back to HEIGHT.

Expected:
- Renaming a column used on both axes leaves both the X and the Y column settings tracking the new name (GROK-19334)
- The rename raises no new console or page error with a formula line configured

### Scenario 3: Logarithmic and inverted axis with a reversed range window

Steps:
1. Note the current console-error and page-error counts.
2. Confirm X is AGE, then in the property panel's **X** section set
   **X Axis Type** to logarithmic.
3. Turn **Invert X Axis** on.
4. Enter **X Min** 60 and **X Max** 20 — a window whose minimum is above its
   maximum on an inverted logarithmic axis.
5. Verify no 'Wrong range' console error and no page error appeared, the
   viewer is still attached, and its viewport settles to a finite rectangle.
6. Revert: clear X Min and X Max, turn **Invert X Axis** off and set
   **X Axis Type** back to linear.

Expected:
- A reversed range window on an inverted logarithmic axis produces no 'Wrong range' console error and no page error (GROK-13110), and the viewer stays attached with a finite viewport

### Scenario 4: Axis type control disabled for a datetime axis

A logarithmic scale is meaningless for a datetime axis, so the axis-type
control is disabled rather than erroring.

Steps:
1. Set the X column to STARTED (a datetime column) through the on-viewer X
   selector.
2. In the property panel's **X** section, verify the **X Axis Type** row is
   shown disabled — the row is dimmed — while the **X** time-unit row is
   active.
3. Set the X column back to AGE (numeric).
4. Verify the pair reverses: the **X Axis Type** row is active and the
   time-unit row is dimmed.
5. Leave X on AGE (the peak configuration of Scenario 1 is restored).

Expected:
- With a datetime column on X the X Axis Type row is disabled (dimmed) and the time-unit row is active; with a numeric column the pair is reversed (GROK-20395)

### Scenario 5: Layout and project persistence at peak configuration

The persistence tail runs at the peak configuration reached above — X=AGE,
Y=HEIGHT, Color=RACE, Size=WEIGHT, Marker Column=SEX — with nothing reverted
before saving.

Steps:
1. Record the current axis and encoding settings of the viewer.
2. Save the view layout.
3. Perturb the view: add another viewer (a Histogram) to the same view and
   clear the **Color** column on the Scatter plot.
4. Re-apply the saved layout and wait for the view to settle.
5. Verify the view's viewer set equals the SAVED set: a Scatter plot is
   present AND the later-added Histogram is absent.
6. Verify the restored Scatter plot kept its configuration — Color is RACE
   again (GROK-18945 guard: a configured Color must survive a layout
   re-apply), together with the recorded X, Y, Size and Marker columns.
7. Delete the probe layout.
8. Save the view as a project through the ribbon **Save** button, dismiss the
   Share dialog that follows, then Close All and reopen the saved project.
9. Verify a Scatter plot viewer is present in the reopened view and its axis
   and encoding settings equal the recorded values — a cross-session
   round-trip.
10. Delete the probe project.

Expected:
- Re-applying the saved layout restores the SAVED viewer set (Scatter plot present, the later-added Histogram absent) with the recorded colorColumnName, xColumnName, yColumnName, sizeColumnName and markersColumnName (GROK-18945)
- Reopening the saved project restores a Scatter plot with the same axis and encoding settings
- The probe layout and project are deleted even when a verification fails — they never leak

## Automation notes

Setup: the viewer handle is
`const sp = grok.shell.tv.viewers.find(v => v.type === 'Scatter plot');` views
are closed via `grok.shell.closeAll()`. The viewer is added via the Toolbox
icon `[name="icon-scatter-plot"]` (a synthetic `.click()` works; the viewer
attaches in about 1.5 s). Resolve the viewer root as
`[...document.querySelectorAll('[name="viewer-Scatter-plot"]')].find(e => !e.closest('.d4-dialog'))`
— the name is NOT unique while the Formula Lines dialog is open, because that
dialog embeds its own preview scatter plot with its own canvases and column
selectors. Do not key any step on the table name: `grok.dapi.files.readCsv`
does not name the table after the file, so the table arrives as `Table` /
`Table (2)`; use the view's dataFrame handle instead.

Scenario 1 — on-viewer selectors: the popup opens on a real (trusted) click or
on a synthetic `mousedown` on `.d4-column-selector-column` inside
`[name="div-column-combobox-x"]`, `-y`, `-color`, `-size` (all LOWERCASE role
suffixes; a synthetic `.click()` does NOT open it). Scope the query to the
viewer root — the property-panel Color editor reuses the same
`div-column-combobox-color` name. Then type the column name with REAL keyboard
events (the popup grid is canvas-rendered, so the column names are not DOM
text) and press Enter, which commits the match on its own. The GROK-18411
guard is expressed by dispatching the `mousedown` on the
`.d4-column-selector-column` label element specifically and waiting for
`.d4-column-selector-backdrop`. The props-echo mirrors are
`sp.props.xColumnName` / `yColumnName` / `colorColumnName` / `sizeColumnName` /
`markersColumnName` plus each selector's `.d4-column-selector-column` text.
A selector with no column assigned is `visibility: hidden` while still present
with a non-zero rect — assert `getComputedStyle(el).visibility`, never
presence. Marker Column is set from the property-panel Marker row.

Scenario 2 — rename: the grid column header context menu has NO `Rename...`
leaf. The working UI route is `[name="div-Column-Properties..."]`, which opens
`[name="dialog-<COLUMN>"]`; type into `input[name="input-New-name--"]` and
confirm with `[name="button-OK"]`. Drive it with a real right-click at the
header's computed coordinates and real typing. The axis propagation is read
from `sp.props.xColumnName` / `yColumnName`, and the configured formula line is
rewritten to the new column name as well. The JS-API equivalent
(`df.col(name).name = newName`) is NOT a substitute here: the guard is about the
UI path, so a failure of the dialog route is a finding to report, not something
to route around.
Formula lines: right-click the data canvas (a synthetic `contextmenu` event
opens the menu) and click `[name="div-Tools---Formula-Lines..."]`; named menu
leaves can be clicked synthetically even while their submenu is unexpanded.
Inside `[name="dialog-Formula-Lines"]` use `[name="button-Add-new"]` then
`div-Line`, and read the round-trip state from `sp.props.formulaLines` (a JSON
string).

Scenario 3: `[name="prop-x-axis-type"]` is a choice row (click the
`[name="prop-view-x-axis-type"]` cell, set the revealed select's `.value`,
dispatch `change`); `[name="prop-invert-x-axis"]` is a checkbox row. The X and
Y min/max rows collide on `name="prop-min"` / `prop-max`, so scope them via
the unique view cells `[name="prop-view-x-min"]` / `prop-view-x-max` and
`.closest('.property-grid-item')`. The viewport is read as `sp.viewport`
(`{x, y, width, height}`) — `sp.props.viewport` is always `null` and must never
be used.

Scenario 4 — disabled state: the ONLY readable disabled signal in the property
grid is an inline `opacity: 0.5` on the row `<tr>` (there is no `disabled`
attribute, no `aria-disabled` and no class), so assert
`document.querySelector('[name="prop-x-axis-type"]').style.opacity`. The value
updates reactively when the governing column changes, with no panel rebuild.
The context menu does NOT mirror this state — `div-Properties...---X---X-Axis-Type`
and its children look identical for numeric and datetime columns — so the
assert must stay on the property-panel row. Note also that the property itself
remains settable from JS while the UI control is greyed; the guard is about
the control, not the property.

Scenario 5: the layout is saved and re-applied via the JS API
(`tv.saveLayout()` / `grok.dapi.layouts.save` / `tv.loadLayout`) — the layout
round-trip needs about 1.5 s to settle after the save and about 4 s after the
re-apply. Clearing a column sets the prop to the EMPTY STRING, not `null` —
compare accordingly when checking the perturbation took effect. The project is
saved through the real ribbon Save button (`[name="button-Save"]` →
`[name="dialog-Save-project"]` → `input#name` → OK) via the
`saveProjectViaUI` helper from `helpers/projects.ts`, because only the UI Save
captures the view layout; the follow-up Share dialog is dismissed via its
CANCEL button. Probe layout and project names carry a `Date.now()` suffix so
concurrent runs never collide, and both are deleted in `finally` teardowns
(`grok.dapi.layouts.delete` / `grok.dapi.projects.delete`) so they are removed
even when an assertion fails.

Console guards: filter the known infra noise of the shared dev server
(WebSocket reconnect errors, the Claude-runtime Docker timeout, the WebGPU
`powerPreference` notice) and the `Canvas2D ... willReadFrequently` notice
that pixel sampling itself provokes; the project-save window additionally
emits the benign cloned-iframe pair whitelisted per canon.

---
{
  "order": 2,
  "datasets": ["System:DemoFiles/demog.csv"]
}
