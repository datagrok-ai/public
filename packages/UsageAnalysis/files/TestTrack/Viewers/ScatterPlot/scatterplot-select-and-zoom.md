---
feature: scatterplot
realizes_atlas:
  - scatterplot.cp.select-and-zoom
realizes:
  - viewers.scatter-plot
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - scatterplot-select-and-zoom-spec.ts
related_bugs:
  - id: GROK-19147
    status: fixed
expected_results:
  - anchor: "Select points by dragging, add a band, deselect, clear"
    expectation: "df.selection.trueCount rises above zero after the first Shift-drag
      rectangle over the cloud, rises again after a second Shift-drag over a
      different band (Shift ADDS to the selection), drops after a
      Ctrl+Shift-drag over part of the selected region (deselect), and returns
      to zero after a click on empty plot background (round-trip)"
  - anchor: "Selection survives a jitter change"
    expectation: "With Jitter Size 20 and Jitter Size Y 15 a Shift-drag selects rows
      (df.selection.trueCount above zero); changing Jitter Size to 30 leaves
      df.selection.trueCount UNCHANGED (GROK-19147 guard — the jitter re-render
      must not detach the selection from its rows), and a following
      Ctrl+Shift-drag over the selected region still lowers the count"
  - anchor: "Zoom, pan and reset the viewport"
    expectation: "With Zoom and Filter set to 'no action', an Alt-drag over a
      sub-region shrinks sp.viewport (both width and height fall below the
      baseline rectangle), a plain drag translates it (x moves while width stays
      constant), a mouse-wheel step in shrinks it further, and the Reset View
      context-menu command restores sp.viewport to the baseline full-data
      rectangle (round-trip)"
  - anchor: "Zoom, pan and reset the viewport"
    expectation: "Dragging the X range slider's minimum handle moves sp.viewport
      (its x rises and its width falls), and the Reset View command restores the
      baseline rectangle"
  - anchor: "Zoom, pan and reset the viewport"
    expectation: "After a wheel zoom, a double-click on a spot inside the plot field
      that carries no marker — confirmed by df.mouseOverRowIdx reading -1 with
      the pointer there — restores sp.viewport to the baseline full-data
      rectangle"
  - anchor: "Keyboard selection and view shortcuts"
    expectation: "After a real click on the plot canvas gives it focus: a real
      Ctrl+A sets df.selection.trueCount to the FILTERED row count, not the
      table's total row count, while a categorical filter is applied; a real
      Ctrl+Shift+A returns it to zero; a real Escape after a fresh Shift-drag
      also returns it to zero; a real H after a wheel zoom restores sp.viewport
      to the baseline rectangle"
  - anchor: "Keyboard selection and view shortcuts"
    expectation: "A real L turns the Lasso Tool on (the viewer's lassoTool setting
      reads on), a Shift-drag closed polygon then raises df.selection.trueCount
      above zero, and a second real L turns the Lasso Tool back off"
---

# Scatter Plot — Point Selection and Viewport Navigation

## Purpose

Verifies the Scatter plot's interactive surface on the demog dataset:
selecting points by dragging a rectangle over the cloud, adding to and
subtracting from that selection, clearing it, keeping it intact while the
markers are re-jittered, and navigating the viewport by zooming, panning,
dragging a range-slider handle and resetting the view. Markers are painted on
canvas, so the gestures are graded by their product effect instead — the
table's selected-row count and the viewer's own viewport rectangle — never by
the pixels they move.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table
   view to load.
3. Add a Scatter plot to the table view via the Toolbox viewer icon and wait
   for it to attach and render.
4. Set the X column to WEIGHT and the Y column to HEIGHT through the on-viewer
   column selectors (the auto-pick already matches on demog, but the
   auto-selection order is not contractual).
5. Record the current browser console-error and page-error counts as the
   baseline for the clean-console checks below.

## Scenarios

### Scenario 1: Select points by dragging, add a band, deselect, clear

Steps:
1. Read the selected-row count — record the baseline value (expected: 0).
2. Shift-drag a rectangle across the middle of the plot area, from about
   three tenths to about seven tenths of the plot in both directions.
3. Verify the selected-row count rose above zero.
4. Shift-drag a second, narrower band over a different part of the cloud
   (lower left).
5. Verify the selected-row count rose again — a Shift-drag ADDS to the
   existing selection rather than replacing it.
6. Ctrl+Shift-drag a rectangle over part of the already selected region.
7. Verify the selected-row count dropped — the covered rows were deselected.
8. Click an empty area of the plot background (a corner with no markers).
9. Verify the selected-row count returned to zero (round-trip). A click ON a
   marker does not clear the selection; only a background click does.

Expected:
- The first Shift-drag raises the selected-row count above zero
- The second Shift-drag raises it again (Shift adds to the selection)
- The Ctrl+Shift-drag lowers it (deselect)
- The background click returns it to zero

### Scenario 2: Selection survives a jitter change

Jitter offsets the drawn markers away from their exact coordinates. A change
of the jitter amount re-renders every marker, and that re-render must not
detach the selection from the underlying rows.

Steps:
1. Open the viewer settings (gear icon on the Scatter plot title bar) and, in
   the **Marker** section, set **Jitter Size** to 20 and **Jitter Size Y** to
   15.
2. Shift-drag a rectangle across the middle of the plot area and verify the
   selected-row count rose above zero; record the count.
3. Change **Jitter Size** to 30 and wait for the plot to re-render.
4. Verify the selected-row count is UNCHANGED from the recorded value
   (GROK-19147 guard).
5. Ctrl+Shift-drag a rectangle over part of the selected region and verify the
   selected-row count still drops — deselection keeps working after the jitter
   change.
6. Revert: click the empty plot background to clear the selection and set both
   **Jitter Size** and **Jitter Size Y** back to 0.

Expected:
- A Shift-drag with non-zero jitter selects rows (count above zero)
- Changing Jitter Size leaves the selected-row count unchanged (GROK-19147)
- A Ctrl+Shift-drag after the jitter change still lowers the count

### Scenario 3: Zoom, pan and reset the viewport

The viewer's default **Zoom and Filter** mode filters the table as the plot is
zoomed. This scenario isolates the viewport, so it switches that mode off
first and restores it at the end; zoom-driven filtering is covered separately.

Steps:
1. In the property panel's **Data** section, set **Zoom and Filter** to
   `no action`.
2. Read the viewer's viewport rectangle — record it as the full-data baseline.
3. Alt-drag a rectangle over a sub-region of the plot (about three tenths to
   six tenths in both directions) — a zoom to that area.
4. Verify the viewport shrank: both its width and its height are below the
   baseline.
5. Drag the plot with the plain left button from its center toward the upper
   left — a pan.
6. Verify the viewport translated: its x moved while its width stayed the
   same.
7. Right-click the plot area and choose **Reset View**.
8. Verify the viewport returned to the recorded baseline rectangle
   (round-trip).
9. Scroll the mouse wheel one step forward over the plot center — a wheel
   zoom.
10. Verify the viewport shrank again, then reset it to the baseline — either
   with the **Reset View** command or by double-clicking an empty spot of the
   plot area, well inside it and away from any marker. A double-click ON a
   marker does not reset the view; it is the point double-click instead.
11. Drag the **minimum handle** of the horizontal range slider along the
   bottom edge of the plot a short distance to the right.
12. Verify the viewport moved accordingly — its x rose and its width fell.
13. Revert: reset the view with **Reset View** and set **Zoom and
   Filter** back to `filter by zoom`.

Expected:
- The Alt-drag zoom shrinks the viewport (width and height below the baseline)
- The plain drag pans it (x moves, width constant)
- The wheel step shrinks it further
- The Reset View command restores the baseline full-data rectangle, and so
  does a double-click on a marker-free spot of the plot area
- Dragging the X range slider's minimum handle raises the viewport x and lowers its width, and the reset restores the baseline

### Scenario 4: Keyboard selection and view shortcuts

The keyboard shortcuts act on the focused plot, so each of them is pressed
after clicking the plot to give it focus. Select-all deliberately covers only
the rows that pass the current filter.

Steps:
1. Open the **Filters** panel from the Toolbox and narrow a categorical
   column (SEX) to a single value; wait for the table to settle and record the
   resulting filtered-row count.
2. Click the plot area once to give it focus.
3. Press Ctrl+A.
4. Verify the selected-row count equals the recorded FILTERED row count — not
   the table's total row count.
5. Press Ctrl+Shift+A and verify the selected-row count returned to zero.
6. Shift-drag a rectangle over the cloud, verify the selected-row count rose
   above zero, then press Escape and verify it returned to zero.
7. Scroll the mouse wheel one step forward over the plot center, then press H.
8. Verify the viewport returned to its pre-zoom rectangle — H resets the view.
9. Press L and verify the viewer's **Lasso Tool** setting reads on.
10. Shift-drag a closed polygon around a group of markers and verify the
   selected-row count rose above zero.
11. Revert: press L again to turn the **Lasso Tool** off, click the empty plot
   background to clear the selection, and reset the filter through the
   **Reset filter** button in the Filters panel; verify the table is
   unfiltered again.

Expected:
- Ctrl+A with a filter applied selects exactly the filtered rows, not the whole table
- Ctrl+Shift+A returns the selected-row count to zero
- Escape after a Shift-drag returns the selected-row count to zero
- H after a wheel zoom restores the pre-zoom viewport
- L turns the Lasso Tool on, a lasso Shift-drag selects rows, and a second L turns it off

## Automation notes

- Narrowing a categorical filter in these steps is driven through the Filter
  Panel's filter-group API rather than by clicking the card's canvas. The guard
  needs a DETERMINISTIC surviving category set, and the card's checkbox list is
  canvas-drawn: a coordinate click can toggle *a* category but cannot choose
  *which* one. Where a guard only needs "exactly one category left, whichever it
  is", the real coordinate click is used instead — see the labels-tooltip
  scenario, which does exactly that. The filter narrowing here is setup for the
  graded signal, never the signal itself.

Setup: the viewer handle is
`const sp = grok.shell.tv.viewers.find(v => v.type === 'Scatter plot');` the
viewer is added via the Toolbox icon `[name="icon-scatter-plot"]` (a synthetic
`.click()` works). Resolve the viewer root as
`[...document.querySelectorAll('[name="viewer-Scatter-plot"]')].find(e => !e.closest('.d4-dialog'))`.
X and Y are set through the on-viewer selectors
(`[name="div-column-combobox-x"]` / `-y`, lowercase role suffixes) with a
synthetic `mousedown` on `.d4-column-selector-column` plus real typing and
Enter; `sp.props.xColumnName` / `yColumnName` is the fallback. Never key a step
on the table name — `readCsv` does not name the table after the file. The two
graded signals are `grok.shell.tv.dataFrame.selection.trueCount` and
`sp.viewport` (a `{x, y, width, height}` rect); `sp.props.viewport` is always
`null` and must never be used.

Scenarios 1-3 — canvas gestures: every pointer gesture on this viewer drives
successfully from real pointer input on `canvas[name="canvas"]` — Shift-drag,
Ctrl+Shift-drag, background click, Alt-drag zoom, plain-drag pan, wheel zoom
and the range-slider handle drag. The double-click reset works the same way,
in both delivery modes, but it is position-sensitive: the viewer resets the
view only when the pointer is NOT over a data point, and over a marker the
gesture becomes the point double-click and leaves the viewport alone. A
double-click aimed at the centre of a dense cloud therefore looks like a
no-op. Aim it well inside the plot area but off the markers, hover the target
first and confirm `df.mouseOverRowIdx` reads `-1`, and grade the reset by
reading `sp.viewport` back. Where a spec needs a reset it can simply rely on
rather than the gesture under test, use the `[name="div-Reset-View"]`
context-menu leaf or a real `H`. Use bubbling events with
`buttons: 1` on mousedown/mousemove and `0` on mouseup, with six to eight
intermediate `mousemove`s about 40 ms apart; `PointerEvent`s work equally.
Recon reference counts on demog (5850 rows, a 430x748 viewer): first
Shift-drag 2968, Ctrl+Shift-drag down to 1690, second Shift-drag band 1692,
background click 0, lasso polygon 3237 — the asserts are directional
(above zero / rose / fell / equals zero), never these exact numbers, because
they depend on the viewer size. Aim a gesture at a known row with
`sp.worldToScreen(x, y)` plus the canvas bounding-rect origin (the helper
returns canvas-local coordinates). A background click only clears the
selection while `resetSelectionOnBackgroundClick` is on (the default); a click
on a marker does not clear it.

Scenario 2: the jitter rows are `[name="prop-jitter-size"]` and
`[name="prop-jitter-size-y"]` — numeric slider rows, driven by focusing the
`property-grid-slider-textbox`, selecting all, typing the value and pressing
Enter (range 0..50). Live recon confirmed the guard holds: at jitter 20/15 a
Shift-drag selected 3549 rows and raising Jitter Size to 30 left
`df.selection.trueCount` at 3549.

Scenario 3: `[name="prop-zoom-and-filter"]` is a choice row — click the
`[name="prop-view-zoom-and-filter"]` cell, set the revealed select's `.value`
and dispatch `change`; the live choices are
`no action | filter by zoom | zoom by filter | pack and zoom by filter` (note
the fourth value's exact wording). Under `no action` the viewport moves while
`df.filter.trueCount` stays at the full count, which keeps this scenario
independent of the debounced zoom-to-filter path. The range sliders are
`svg[name="x-slider"]` / `y-slider` with `[name="min-handle"]`,
`[name="max-handle"]` and `[name="pan-handle"]` children; they are
`visibility: hidden` until hover but drag synthetically while hidden (a 30 px
drag of the X minimum handle moved the viewport from x 30.2 / width 161.4 to
x 50 / width 141.5 in recon). Reset through the Reset View leaf or a real `H`;
never through `sp.zoom(x1, y1, x2, y2)` — setting the viewport from JS would
replace the gesture under test rather than exercise it.

Scenario 4 — keyboard: the shortcuts need REAL key events AND canvas focus —
a synthetic `KeyboardEvent` is ignored, and a real (trusted) click on the
canvas must precede the keys so `document.activeElement` is the canvas. Use
`page.mouse.click` on the canvas center followed by `page.keyboard.press`.
Ctrl+A selects the FILTERED rows only (recon: 38 of 100 rows on a filtered
table), so the assert compares against `df.filter.trueCount` read at that
moment, not against `df.rowCount`. `L` toggles `sp.props.lassoTool` and `R`
would toggle the regression line (out of scope here). The Filters panel is
opened through the Toolbox `[name="div-section--Filters"]`, which may be
present but invisible — guard on `offsetParent` and fall back to
`tv.getFiltersGroup()`; the panel builds slowly, so poll for its `.d4-filter`
card count before driving it. The reset control is
`[name="viewer-Filters"] [name="icon-arrow-rotate-left"]` (no confirmation
dialog), and `[name="span-filtered"]` in the status bar is ABSENT while the
table is unfiltered — check presence before reading it.

Console guards: filter the known infra noise of the shared dev server
(WebSocket reconnect errors, the Claude-runtime Docker timeout, the WebGPU
`powerPreference` notice) before comparing error counts.

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/demog.csv"]
}
