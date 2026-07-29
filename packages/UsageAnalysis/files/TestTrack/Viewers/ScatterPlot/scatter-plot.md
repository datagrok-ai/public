---
feature: scatterplot
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas:
  - scatterplot-lines-by-overrides-color-split
realizes: []
realized_as:
  - scatter-plot-spec.ts
related_bugs:
  - id: GROK-13533
    status: fixed
expected_results:
  - anchor: "Axis histograms"
    expectation: "Turning Show X Histogram on, then Show Y Histogram on, then
      changing Histogram Bins each changes the data canvas rendering, and
      turning both histograms back off restores the rendering toward the
      baseline"
  - anchor: "Grid lines, axes and selector visibility"
    expectation: "Turning Show Vertical Grid Lines and Show Horizontal Grid Lines
      off changes the data canvas rendering, and turning Show X Axis and Show Y
      Axis off changes it again — the axes and grid lines are drawn on the
      canvas"
  - anchor: "Grid lines, axes and selector visibility"
    expectation: "Turning Show X Selector and Show Y Selector off makes the
      on-viewer X and Y column selectors invisible while the elements stay in
      place — the check reads the computed visibility, not element presence,
      because a presence check would pass regardless (GROK-13533)"
  - anchor: "Grid lines, axes and selector visibility"
    expectation: "Turning all four visibility settings and both grid-line settings
      back on restores the selectors to visible and the data canvas rendering
      toward the baseline"
  - anchor: "Whiskers"
    expectation: "Setting X Whisker Min and Max, then Y Whisker Min and Max, changes
      the data canvas rendering and raises no console or page error; clearing
      all four whisker columns restores the rendering toward the baseline"
  - anchor: "Context menu"
    expectation: "A right-click on the plot area opens the context menu and it
      carries the Reset View, Lasso Tool, Tools and Properties entries;
      dismissing the menu closes it and leaves the viewer alive with no console
      or page error"
  - anchor: "Title and description"
    expectation: "Setting the Title renders that text on the viewer and setting the
      Description renders that text inside the viewer element; clearing both
      removes the rendered texts"
  - anchor: "Lines By overrides the color column when splitting connecting lines"
    expectation: "The Lines By row is greyed out while no Lines Order column is set
      and becomes enabled once one is"
  - anchor: "Lines By overrides the color column when splitting connecting lines"
    expectation: "With Lines Order set and Color set to RACE, pointing Lines By at
      RACE leaves the data canvas rendering unchanged while pointing it at SEX,
      which has fewer categories, changes the rendering"
  - anchor: "Lines By overrides the color column when splitting connecting lines"
    expectation: "Clearing Lines By returns the rendering to the color-split
      rendering, and clearing Lines Order and Color returns it toward the
      baseline"
---
# Scatter plot tests

## Purpose

Verifies the Scatter plot's secondary settings surface on the demog dataset:
the axis histograms, grid lines and the visibility toggles of the axes and the
on-viewer column selectors, whiskers, the context menu, the title and
description, and the way connecting lines are split into series. The focused
scenarios in this folder own the primary workflows
(axes and encoding, selection and zoom, zoom-driven filtering, the legend,
regression and formula lines, labels and tooltip); nothing here duplicates
them.

Most of what this file touches is painted on the plot canvas, so those
settings are checked as a settle-gated rendering change plus an error-free
floor rather than judged by appearance. The one exception is the pair of
selector-visibility toggles, which have a readable computed-visibility signal.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Scatter plot to the current table view via the Toolbox viewer icon.
4. Using the on-viewer column selectors, set X to WEIGHT and Y to HEIGHT.

## Scenarios

### Scenario 1: Axis histograms

Steps:
1. Measure the plot's data canvas rendering — the baseline.
2. Open the viewer settings and, in the **Axes** section, turn
   **Show X Histogram** on; verify the rendering changed.
3. Turn **Show Y Histogram** on; verify the rendering changed again.
4. Set **Histogram Bins** to 20; verify the rendering changed again.
5. Revert: turn both histograms off and restore **Histogram Bins**; verify the
   rendering moved back toward the baseline.

Expected:
- Each histogram step changes the data canvas rendering
- Turning both histograms back off restores the rendering toward the baseline

### Scenario 2: Grid lines, axes and selector visibility

Steps:
1. Measure the plot's data canvas rendering — the baseline.
2. In the viewer settings **X** section turn **Show Vertical Grid Lines** off,
   and in the **Y** section turn **Show Horizontal Grid Lines** off; verify the
   rendering changed.
3. Turn **Show X Axis** and **Show Y Axis** off; verify the rendering changed
   again.
4. Turn **Show X Selector** and **Show Y Selector** off.
5. Verify the on-viewer X and Y column selectors became invisible while the
   elements themselves stay in place — read the computed visibility of the
   selectors, not their presence (GROK-13533 guard; a presence check would
   pass either way and would be false green).
6. Revert: turn **Show X Selector**, **Show Y Selector**, **Show X Axis**,
   **Show Y Axis** and both grid-line settings back on.
7. Verify the selectors are visible again and the rendering moved back toward
   the baseline.

Expected:
- Turning the grid lines off changes the data canvas rendering, and turning the axes off changes it again
- Turning the selector visibility settings off makes the on-viewer X and Y selectors invisible while their elements remain in place
- Turning every setting back on restores the selectors to visible and the rendering toward the baseline

### Scenario 3: Whiskers

Steps:
1. Note the current browser console-error and page-error counts, and measure
   the plot's data canvas rendering — the baseline.
2. In the viewer settings **X** section set **X Whisker Min** to AGE and
   **X Whisker Max** to WEIGHT.
3. In the **Y** section set **Y Whisker Min** to HEIGHT and **Y Whisker Max**
   to WEIGHT.
4. Verify the rendering changed against the baseline and no new console or
   page error appeared.
5. Revert: clear all four whisker columns and verify the rendering moved back
   toward the baseline.

Expected:
- Configuring the X and Y whisker columns changes the data canvas rendering and raises no console or page error
- Clearing all four whisker columns restores the rendering toward the baseline

### Scenario 4: Context menu

Steps:
1. Right-click the plot area — the context menu opens.
2. Verify the menu carries the **Reset View**, **Lasso Tool**, **Tools** and
   **Properties...** entries.
3. Dismiss the menu.
4. Verify the menu closed, the viewer is still alive, and no console or page
   error appeared.

Expected:
- The context menu opens with the Reset View, Lasso Tool, Tools and Properties entries
- Dismissing the menu closes it and leaves the viewer alive with no console or page error

### Scenario 5: Title and description

Steps:
1. Open the viewer settings and set **Title** to "Test Plot"; verify that text
   renders on the viewer.
2. Set **Description** to "Test description"; verify that text renders inside
   the viewer element.
3. Revert: clear both the **Title** and the **Description**; verify both
   rendered texts disappear.

Expected:
- Setting the Title renders that text on the viewer and setting the Description renders that text inside the viewer element
- Clearing both removes the rendered texts

### Scenario 6: Lines By overrides the color column when splitting connecting lines

Connecting lines are drawn in the order given by the **Lines Order** column and
split into series by the **Lines By** column when one is set, and by the color
column when it is not — so **Lines By** takes over from the color column. Both
renderings are painted on the canvas, so the two split rules are told apart by
comparing them against each other rather than by inspecting a line.

Steps:
1. In the viewer settings **Data** section, verify the **Lines By** row is greyed
   out while no **Lines Order** column is set.
2. Measure the plot's data canvas rendering — the baseline.
3. In the **Color** section set the color column to RACE, and read the number of
   categories of RACE and of SEX from the data — SEX must have fewer.
4. Set **Lines Order** to AGE; verify the **Lines By** row is no longer greyed
   out and the rendering changed — the connecting lines are now drawn, split by
   the color column.
5. Measure the rendering — the color-split reading.
6. Set **Lines By** to RACE, the same column the color uses; verify the
   rendering is unchanged against the color-split reading.
7. Set **Lines By** to SEX; verify the rendering changed against the color-split
   reading.
8. Revert: clear **Lines By** and verify the rendering returned to the
   color-split reading; then clear **Lines Order** and the color column and
   verify the rendering moved back toward the baseline.

Expected:
- The Lines By row is greyed out until a Lines Order column is set
- With Lines Order set, pointing Lines By at the color column leaves the rendering unchanged, and pointing it at a column with fewer categories changes it
- Clearing Lines By restores the color-split rendering, and clearing Lines Order and the color column restores the rendering toward the baseline

## Automation notes

- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Scatter plot')`; the viewer is
  added via the Toolbox icon `[name="icon-scatter-plot"]` (synthetic `.click()`
  works). Resolve the viewer root as the `[name="viewer-Scatter-plot"]` element
  not inside a `.d4-dialog`. The settings panel opens from the gear icon found
  through `root.closest('.panel-base')`, never a fixed parent-hop count.
- Rendering readings are settle-gated non-white pixel fractions over
  `canvas[name="canvas"]` (the data canvas — markers, axes, grid lines,
  histograms, whiskers). Thresholds are calibrated live at spec-authoring time,
  never carried over as settled constants: they depend on the viewer size and
  the device pixel ratio. The pixel sampling itself triggers the browser's
  `willReadFrequently` console advisory — whitelist it in the console guards.
- Scenario 2, GROK-13533: `showXSelector` / `showYSelector` set to false only
  flip the selector element's computed `visibility` to `hidden` — the element
  stays in the DOM with a non-zero rect and a live `offsetParent`, so presence
  and `offsetParent` are both false-green signals. Assert
  `getComputedStyle(el).visibility` on
  `[name="div-column-combobox-x"]` / `-y` inside the viewer root. Note the same
  selector names are reused by the property-panel column editors, so scope to
  the viewer root. The axis and grid-line rows are
  `[name="prop-show-x-axis"]`, `[name="prop-show-y-axis"]`,
  `[name="prop-show-vertical-grid-lines"]`,
  `[name="prop-show-horizontal-grid-lines"]`, with the mirrored menu leaves
  `div-Controls---Show-X-Axis` / `Show-Y-Axis`.
- The property panel duplicates `prop-min` / `prop-max` across X and Y — scope
  every axis-scoped row through the unique `prop-view-*` cells and
  `.closest('.property-grid-item')`. Column-valued rows (whiskers, label
  columns) embed a `div-column-combobox-<role>` that opens on a synthetic
  `mousedown` and commits on real typing plus `Enter`.
- The context menu opens from a synthetic `contextmenu` MouseEvent on
  `canvas[name="canvas"]`; scope the assertions to the `.d4-menu-popup`
  container, because a bare `.d4-menu-item` query also matches the application
  menubar. Every leaf carries a `name="div-<Path---With---Separators>"`
  attribute, so the entries are checked by name
  (`div-Reset-View`, `div-Lasso-Tool`, `div-Tools`, `div-Properties...`).
- Scenario 6: the connecting-line rows live in the **Data** category, not in
  Lines, and their column comboboxes carry a double-dash role token —
  `[name="prop-lines-order"]` embeds `div-column-combobox-lines--order` and
  `[name="prop-lines-by"]` embeds `div-column-combobox-lines--by`. The greyed
  state of the Lines By row is the same inline `opacity: 0.5` regime the rest of
  the property grid uses. Pointing Lines By at the color column is what makes
  the check discriminating: it isolates the change of the split key from the
  mere act of setting a property, because the same split key must produce the
  same picture.
- Title and description are DOM text reads scoped to the viewer element and its
  title bar. Note the scatter-plot title bar carries no close icon — close the
  viewer through the hamburger menu or the JS API if a scenario ever needs it.
- Console guards subtract a baseline and whitelist infra noise unrelated to the
  viewer (the shared dev server's WebSocket reconnect errors, the Claude
  runtime container timeout, the WebGPU `powerPreference` advisory and the
  `willReadFrequently` advisory caused by the pixel sampling).

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/demog.csv"]
}
