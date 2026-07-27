---
feature: densityplot
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - density-plot-spec.ts
related_bugs: []
expected_results:
  - anchor: "Show/Hide Color Scale"
    expectation: "Toggling Show Color Scale off and back on round-trips the
      showColorScale prop (true -> false -> true) with no page or console
      errors; the color scale itself is canvas-drawn, so the check is the prop
      round-trip plus an error-free repaint, not a DOM assert"
  - anchor: "Axis Visibility"
    expectation: "Hiding and re-showing the X and Y axes round-trips the showXAxis /
      showYAxis props with no page or console errors; the axis lines and tick
      labels are canvas-drawn, so there is no DOM signal beyond the prop
      round-trip and the error-free floor"
  - anchor: "Show/Hide Selectors and Bin Slider"
    expectation: "Hiding the X and Y selectors switches their DOM elements' computed
      visibility to hidden (the elements stay in DOM with a nonzero rect);
      re-showing restores computed visibility to visible; the bin slider input
      stays present in DOM through its own hide/show toggle with the
      showBinSelector prop round-tripping"
  - anchor: "Title and Description"
    expectation: "The title renders in the panel titlebar and the description
      renders inside the viewer element at the configured position; clearing
      both removes them"
  - anchor: "Row Source Filtering"
    expectation: "Setting Filter to ${AGE} > 30 repaints the plot (fewer rows are
      binned) or holds an error-free floor, and clearing it restores; switching
      Table to spgi-100 re-binds the viewer to the new table's numeric columns
      with no page or console errors, and switching back to demog restores"
---

# Density plot tests

## Purpose

Verifies the Density Plot's secondary display surface: showing and hiding the
color scale, the axes, the on-viewer column selectors and the bin slider,
title and description, and row-source filtering with a table switch. The core
workflows — axis column picking, binning, color mapping, persistence,
selection, zoom, bin-to-range and axis scale — are owned by the focused
scenarios in this folder and are not repeated here. Most of these settings
only change how the canvas is painted, so they are checked as prop
round-trips with an error-free floor; the selector and bin-slider toggles are
the ones with a readable DOM outcome.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Density plot to the current table view via the Toolbox viewer icon.

## Show/Hide Color Scale

1. Go to the Context Panel > Style and set **Show Color Scale** to false — the
   color scale disappears from the plot
2. Set **Show Color Scale** to true — the color scale is back
3. Verify the prop round-trips and no page or console errors appeared

## Axis Visibility

1. Go to the Context Panel > X and set **Show X Axis** to false — the X axis
   disappears
2. Go to the Context Panel > Y and set **Show Y Axis** to false — the Y axis
   disappears
3. Set **Show X Axis** to true and **Show Y Axis** to true — both axes are back
4. Verify both props round-trip and no page or console errors appeared

## Show/Hide Selectors and Bin Slider

1. Go to the Context Panel > X and set **Show X Selector** to false — the
   on-viewer X column selector is hidden
2. Go to the Context Panel > Y and set **Show Y Selector** to false — the
   on-viewer Y column selector is hidden
3. Go to the Context Panel > Misc and set **Show Bin Selector** to false — the
   on-viewer bin slider is hidden
4. Set **Show X Selector**, **Show Y Selector**, and **Show Bin Selector**
   back to true — all three controls are visible again
5. Verify the selectors' visibility actually changed in the DOM in both
   directions

## Title and Description

1. Go to the Context Panel > Description, enable **Show Title** and set
   **Title** to "Density Distribution" — the title appears in the panel
   titlebar
2. Set **Description** to "AGE vs HEIGHT density" — the description text
   appears inside the viewer
3. Set **Description Visibility Mode** to Always and **Description Position**
   to Bottom — the description stays visible at the bottom
4. Clear both fields — the title and the description disappear

## Row Source Filtering

1. Go to the Context Panel > Data and set **Filter** to `${AGE} > 30` — the
   plot re-bins over the matching rows only
2. Clear **Filter** — the full row set is binned again
3. Open spgi-100
4. Go back to the demog table and click the density plot
5. On the Context Panel > Data set **Table** to spgi-100 — the viewer
   re-binds to the new table's numeric columns; verify no errors
6. Set **Table** back to demog — the original binding is restored
7. Close All

## Automation notes

Setup: the viewer handle is
`const dp = grok.shell.tv.viewers.find(v => v.type === 'Density plot');` views
are closed via `grok.shell.closeAll()`; the viewer is added via the Toolbox
icon `[name="icon-density-plot"]`.

Show/Hide Color Scale: the scale is drawn on the single
`canvas[name="canvas"]` — the viewer's DOM element count does not change when
it is toggled and the scale's min/max numbers are not DOM text, so the
assertable signal is the `showColorScale` prop round-trip plus the error-free
floor (a settle-gated pixel delta is optional and must be live-calibrated).
The checkbox row is `[name="prop-show-color-scale"]`.

Axis Visibility: the axis lines and tick labels are canvas-drawn with no DOM
counterpart, so this section is a prop round-trip
(`showXAxis` / `showYAxis`, rows `[name="prop-show-x-axis"]` /
`[name="prop-show-y-axis"]`) over the no-error floor.

Show/Hide Selectors and Bin Slider: hiding a selector keeps its element in
DOM with a nonzero rect — only `visibility: hidden` changes — so the DOM
assert reads `getComputedStyle(el).visibility` on
`[name="div-column-combobox-x"]` / `[name="div-column-combobox-y"]` in both
directions; presence or offsetParent checks are false-green. The bin slider
is the `<input type="range">` inside the viewer root; it likewise stays in
DOM through the toggle, so its signal is the `showBinSelector` prop
round-trip. The rows are `[name="prop-show-x-selector"]`,
`[name="prop-show-y-selector"]`, `[name="prop-show-bin-selector"]`.

Title and Description: these are base-viewer properties shared by all viewers
and sit outside the density-plot atlas denominator (hence `realizes_atlas` is
empty on this file). The title is read from the panel titlebar
(`.panel-titlebar-text`) and the description from the element inside the
viewer root; both are DOM-readable.

Row Source Filtering: the filter formula is driven via `dp.props.filter` and
the table switch via `dp.props.table` (row `[name="prop-table"]`); spgi-100
is opened from `System:AppData/Chem/tests/spgi-100.csv`. The filter's render
effect is a settle-gated pixel delta candidate but may honestly stay a
no-error floor; the table switch is asserted as no new page or console errors
with the viewer root still attached. Console guards filter infra noise (the
shared dev server logs unrelated `WebSocket ... 503` reconnect errors).

---
{
  "order": 11,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
