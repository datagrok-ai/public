---
feature: matrixplot
realizes_atlas:
  - matrixplot.cp.configure-axes-inner-type
  - matrixplot-axes-drive-inner-grid
realizes:
  - viewers.matrix-plot
realized_as:
  - matrixplot-configure-axes-inner-type-spec.ts
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-16473
    status: fixed
  - id: GROK-20438
    status: fixed
  - id: GROK-10925
    status: fixed
expected_results:
  - anchor: "Default state — the demog auto-pick"
    expectation: "After adding the viewer through the Toolbox icon, 16 inner cell
      canvases render (a 4x4 matrix), the X and Y column sets both read AGE,
      HEIGHT, WEIGHT, STARTED (the demog auto-pick), and Cell Plot Type reads
      its default, Density plot"
  - anchor: "Change the X column set through the Select columns dialog"
    expectation: "After unchecking WEIGHT and STARTED in the X Select columns dialog
      and clicking OK, the matrix re-tiles to 2x4 (8 inner cell canvases) and
      the column label texts above the matrix list exactly AGE and HEIGHT
      (GROK-20438 guard) — the labels are read from the DOM, not from the
      property"
  - anchor: "Cycle the column sets — console-error guard"
    expectation: "Cycling the X and Y column sets through several configurations and
      back to the full 4x4 raises no new browser console or page error (error
      delta == 0, GROK-16473 guard) and the matrix ends re-tiled to 16 cells"
  - anchor: "Switch the Cell Plot Type"
    expectation: "Switching Cell Plot Type from Density plot to Scatter plot changes
      an off-diagonal cell's settle-gated pixel measurement (render delta above
      the noise floor), and switching back to Density plot changes it again"
  - anchor: "Layout and project persistence at peak configuration — layout round-trip"
    expectation: "After saving the view layout, adding a Scatter plot viewer, and
      re-applying the saved layout, the view's viewer set equals the SAVED set
      (a Matrix plot is present AND the later-added Scatter plot is absent) AND
      the restored Matrix plot carries the configured xColumnNames,
      yColumnNames, and cellPlotType; the probe layout is deleted afterwards
      even on failure"
  - anchor: "Layout and project persistence at peak configuration — project save /
      Close All / reopen"
    expectation: "After saving the view as a project through the ribbon Save button,
      Close All, and reopening the saved project, a Matrix plot viewer is
      present in the reopened view AND its xColumnNames, yColumnNames, and
      cellPlotType equal the persisted values (GROK-10925 guard, a cross-session
      round-trip); the probe project is deleted afterwards even on failure"
---

# Matrix Plot — Column Sets, Cell Plot Type, Persistence

## Purpose

Verifies the Matrix plot's main configuration surface on the demog dataset:
the default auto-picked column sets, changing the X column set through the
real Select columns dialog with the rendered column labels as the readable
outcome, error-free re-tiling while the column sets are cycled, switching the
inner cell plot type between Density plot and Scatter plot, and the survival
of the configured viewer across a saved layout and a saved project. The inner
cells draw to canvas, so their repaint is measured as a settle-gated pixel
change; the column labels are plain DOM text and are asserted directly.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Matrix plot to the current table view via the Toolbox viewer icon.

## Scenarios

### Scenario 1: Default state — the demog auto-pick

Steps:
1. Wait for the Matrix plot to attach and render.
2. Count the inner cell canvases — 16 cells render (a 4x4 matrix).
3. Read the X and Y column sets: both are AGE, HEIGHT, WEIGHT, STARTED (the
   numerical and datetime demog columns the viewer picks by itself).
4. Read **Cell Plot Type** — it is at its default, Density plot.

Expected:
- 16 inner cell canvases render, the X and Y column sets both read AGE, HEIGHT, WEIGHT, STARTED, and Cell Plot Type reads Density plot

### Scenario 2: Change the X column set through the Select columns dialog

Steps:
1. Open the viewer settings (gear icon on the Matrix plot title bar) and find
   the **X** row in the Data section.
2. Click the **...** button on the X row — the **Select columns** dialog
   opens, listing the four eligible demog columns with a checked checkbox per
   column and a "4 checked" counter.
3. Uncheck **WEIGHT** and **STARTED**, leaving AGE and HEIGHT checked.
4. Click **OK** to commit.
5. Verify the matrix re-tiled to 2x4 — 8 inner cell canvases.
6. Read the column label texts rendered above the matrix and verify they list
   exactly AGE and HEIGHT — the labels must match the chosen set (GROK-20438
   guard).

Expected:
- The matrix re-tiles to 8 inner cell canvases and the DOM column labels above the matrix list exactly AGE and HEIGHT

### Scenario 3: Cycle the column sets — console-error guard

Steps:
1. Note the current browser console-error and page-error counts.
2. Change the X column set to AGE, HEIGHT, WEIGHT; wait for the re-tile.
3. Change the Y column set to AGE, HEIGHT; wait for the re-tile.
4. Change both sets back to the full AGE, HEIGHT, WEIGHT, STARTED; wait for
   the re-tile — the matrix is 4x4 (16 cells) again.
5. Verify no new console or page error appeared across steps 2–4
   (GROK-16473 guard — the historical NoSuchMethodError class).

Expected:
- Cycling the X and Y column sets raises no new console or page error and the matrix ends re-tiled to 16 cells

### Scenario 4: Switch the Cell Plot Type

Steps:
1. Measure a settle-gated pixel reading of an off-diagonal cell (for example
   the X=HEIGHT, Y=AGE cell) in the Density plot state.
2. Set **Cell Plot Type** to Scatter plot; wait for the redraw.
3. Measure the same cell again and verify the reading changed — the
   off-diagonal cell repainted as a scatter plot. Diagonal cells keep
   rendering histograms throughout; this is incidental and is not asserted
   per cell.
4. Set **Cell Plot Type** back to Density plot; wait for the redraw and
   verify the cell's reading changed again.

Expected:
- The off-diagonal cell's settle-gated pixel measurement changes on the switch to Scatter plot and changes again on the switch back to Density plot

### Scenario 5: Layout and project persistence at peak configuration

This scenario persists the configuration exactly as the previous scenarios
left it (nothing is reverted first), so the round-trips restore concrete
non-default values.

Steps:
1. Configure the known peak state: X columns AGE, HEIGHT; Y columns AGE,
   HEIGHT, WEIGHT; **Cell Plot Type** Scatter plot. Record these values.
2. Layout round-trip: save the view layout, add a Scatter plot viewer, then
   re-apply the saved layout.
3. Verify the view's viewer set equals the SAVED set: a Matrix plot is
   present AND the later-added Scatter plot is absent.
4. Verify the restored Matrix plot kept its configuration: the recorded X
   columns, Y columns, and Cell Plot Type.
5. Delete the probe layout.
6. Project round-trip: save the view as a project through the ribbon **Save**
   button, dismiss the Share dialog, Close All, then reopen the saved
   project.
7. Verify a Matrix plot viewer is present in the reopened view AND its X
   columns, Y columns, and Cell Plot Type equal the recorded values — a
   cross-session round-trip (GROK-10925 guard).
8. Delete the probe project.

Expected:
- Re-applying the saved layout restores the SAVED viewer set (Matrix plot present, Scatter plot absent) with the configured xColumnNames/yColumnNames/cellPlotType
- Reopening the saved project restores a Matrix plot with the same xColumnNames/yColumnNames/cellPlotType
- The probe layout and project are deleted even when a verification fails — they never leak

## Automation notes

- The viewer is added via the Toolbox icon `[name="icon-matrix-plot"]`
  (synthetic `.click()` works); the handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Matrix plot')`. The re-tile
  signal is the count of
  `[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer` elements.
- Scenario 2 drives the real dialog: the `...` button on the
  `[name="prop-x"]` row opens `[name="dialog-Select-columns..."]` only on a
  REAL click (synthetic clicks are swallowed or hit the shell view-selector).
  The per-column checkboxes are canvas-drawn — toggling needs a REAL
  coordinate click on the checkbox cell near the grid's right edge
  (`page.mouse.click`, row height ~28 CSS px under a ~24 px header;
  re-measure the geometry at spec time). The label-All / label-None buttons
  and OK / CANCEL are DOM and synthetic-clickable; All / None ignore the
  active search filter; toggles live-apply before OK.
- WARNING: NEVER leave the selection empty. Clicking None with an empty
  selection throws a RangeError and bricks the viewer (props read back but
  the matrix never re-tiles again; only close + re-add recovers). This is a
  known open ticket and is deliberately NOT exercised — every set change in
  this scenario keeps at least one column checked.
- The GROK-20438 label read: the X labels are unnamed leaf `<div>`s in the
  strip above the cells — collect leaf divs inside
  `[name="viewer-Matrix-plot"]` whose trimmed text equals a column name and
  compare the text set to the chosen columns.
- Scenario 3 cycles the sets via the JS API
  (`mp.props.xColumnNames = [...]` / `yColumnNames`) — the same
  onLookChanged re-tile path the dialog commits through; the dialog itself
  is exercised once with real clicks in Scenario 2. The error baseline
  subscribes `page.on('pageerror')` and counts console errors, filtering the
  dev server's unrelated WebSocket reconnect noise.
- Cell Plot Type is set via the `[name="prop-view-cell-plot-type"]` row: a
  click reveals a `select` editor; set `.value` and dispatch `change`
  (synthetic works). The pixel measurement is a settle-gated non-white pixel
  count over the target cell canvas, re-read until two consecutive readings
  agree; cell canvas backing-store pixels differ from CSS px (device pixel
  ratio) — scale coordinates by `rect.width / canvas.width`.
- Scenario 5: the layout is saved and re-applied via the JS API
  (`tv.saveLayout()` / `grok.dapi.layouts.save` / `tv.loadLayout`); the
  project is saved through the real ribbon Save button
  (`[name="button-Save"]`) because only the UI Save captures the view
  layout; the Share dialog is dismissed via its CANCEL button. Probe names
  carry a `Date.now()` suffix; the probe layout and project are deleted in
  `finally` teardowns (`grok.dapi.layouts.delete` /
  `grok.dapi.projects.delete`).
- The project-publish preview clones the live view into an offscreen iframe;
  that clone emits "Unable to find element in cloned iframe" plus the Dart
  "NullError: method not found ... on null" + "Stack trace" pair against the
  detached ProjectMeta view. The console guard whitelists this class only
  inside the ribbon project-save window (a flag around `saveProjectViaUI`);
  everywhere else — notably the GROK-16473 set-cycling scenario — the same
  Dart error class is the regression signal and must reach the guard.

---
{
  "order": 12,
  "datasets": ["System:DemoFiles/demog.csv"]
}
