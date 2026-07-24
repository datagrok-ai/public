---
feature: pcplot
realizes_atlas:
  - pcplot.cp.normalize-and-density
  - pcplot-normalize-density-overlay
realizes:
  - viewers.pc-plot
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: github-1546
    status: fixed
realized_as:
  - pcplot-normalize-density-spec.ts
expected_results:
  - anchor: "Scenario 1 — switch vertical scale global then back"
    expectation: "Switching the vertical scale to a shared global scale and back to
      Normalized raises no page or console error and the PC Plot root element is
      still in the DOM"
  - anchor: "Scenario 1 — switch vertical scale global then back"
    expectation: "The settle-gated canvas ink counts in the Normalized and Global
      states differ by more than 500 pixels (either direction), and the round-trip
      back to Normalized returns the count to within 2000 pixels of the original"
  - anchor: "enable density, cycle circles/box/violin styles"
    expectation: "Before density is enabled, pc.props.densityStyle reads back its
      default value 'circles'"
  - anchor: "Scenario 2 — enable density, cycle circles/box/violin styles, drive per-part toggles"
    expectation: "Enabling density, cycling the style through circles, box plot and
      violin, and driving the per-part box/violin toggles (Show Median,
      Interquartile Range, Mean Cross, upper/lower whisker dash, Show Circles) and
      the Bins change raises no page or console error and the PC Plot root element
      is still in the DOM"
  - anchor: "Scenario 2 — enable density, cycle circles/box/violin styles, drive per-part toggles"
    expectation: "With all polylines hidden and the selection empty, enabling the
      density overlay in the Box Plot style adds ink over the sparse axes-only
      floor: the settle-gated canvas ink count with density on exceeds the
      hidden-lines floor by more than 1000 pixels"
  - anchor: "Scenario 3 — density recalculates on normalization double-toggle and AGE log scale (github-1546)"
    expectation: "github-1546 regression guard: with density on, double-toggling
      normalization and then switching AGE to a log scale raises no page or
      console error, a follow-up page.evaluate resolves (the overlay
      recalculated and did not freeze), and the PC Plot root element is still in
      the DOM"
  - anchor: "Scenario 3 — density recalculates on normalization double-toggle and AGE log scale (github-1546)"
    expectation: "After the normalization double-toggle with density still on,
      measured with all polylines hidden, the settle-gated canvas ink count
      stays more than 1000 pixels above the hidden-lines floor recorded in
      Scenario 2 — the density overlay is still painted, not stale-empty"
---

# PC Plot — Normalization and Density Overlay

## Purpose

Verifies the PC Plot's vertical scale (per-column normalization versus a shared
global scale) and its density overlay (circles, box plot, violin, and the
per-part box and violin components), including a regression guard for
github-1546, where the density overlay kept stale geometry after normalization
or log-scale changes. The plot draws entirely to canvas, so each change is
checked in two honest ways: it completes without page or console errors with
the viewer still on the page, and it visibly repaints the canvas by a clear
margin. The scenario deliberately does not judge whether the resulting picture
"looks right".

## Setup

1. Close all open views.
2. Open the demog dataset `System:DemoFiles/demog.csv` and wait for the table view.
3. Add a PC Plot viewer to the table view (Toolbox or via ribbon).
4. Assign AGE, HEIGHT and WEIGHT as the three axes and wait for the re-render.
5. Note the current page-error and console-error counts so each scenario can
   confirm that no new errors appeared.

## Scenarios

### Scenario 1: Vertical scale switched to a shared global scale

The vertical scale is driven by the **Normalize Each Column** checkbox
(Context Panel > Value section). Checked (the default), every axis spans its own
column's min-max; unchecked, all axes share one global scale so values are
comparable across columns (on demog, AGE visibly compresses because WEIGHT sets
the shared range). Switch to the non-default state first, then revert (round-trip).

Steps:
1. Note the current page-error and console-error counts, and measure the ink
   drawn on the plot canvas in the default Normalized state.
2. In the Context Panel > Value section, uncheck **Normalize Each Column**
   (shared global scale), wait for the plot to redraw, and measure the canvas
   ink in the Global state.
3. Re-check **Normalize Each Column** (per-column normalization), wait for the
   redraw, and measure the canvas ink again (round-trip).
4. Check that the Normalized and Global measurements clearly differ in either
   direction (the shared scale redistributes the polylines) and that the
   round-trip measurement lands close to the original Normalized one (exact
   margins in Automation notes).
5. Check that no new page or console error appeared and that the PC Plot
   viewer is still present on the page.

Expected:
- Switching the vertical scale to a shared global scale and back to Normalized raises no page or console error and the PC Plot root element is still in the DOM.
- The settle-gated canvas ink counts in the Normalized and Global states differ by more than 500 pixels (either direction), and the round-trip back to Normalized returns the count to within 2000 pixels of the original.

### Scenario 2: Density overlay enabled and its style cycled

Steps:
1. Note the current page-error and console-error counts, and check that the
   density style shows its default value, Circles, before any change. Then
   hide all polylines: clear the selection and turn **Show All Lines** off
   (right-click the plot > Selection > Show All Lines), and measure the ink
   on the now-sparse canvas — axes and labels only. This axes-only floor is
   reused by Scenario 3.
2. Enable the density overlay (**Show Density** in the viewer's Context
   Panel), set its style to Box Plot (right-click the plot > Style >
   box plot), wait for the redraw, and measure the canvas ink with density
   on. Check that it clearly exceeds the axes-only floor — the box-plot
   shapes are the dominant ink on the sparse canvas. Then turn
   **Show All Lines** back on (round-trip).
3. Set the density style to Circles (right-click the plot > Style > circles).
4. Set the density style to Box Plot (right-click the plot > Style >
   box plot).
5. Set the density style to Violin (right-click the plot > Style >
   violin plot).
6. Switch the style back to Box Plot, then in the Context Panel's Box Plot
   section toggle each per-part component: **Show Median**, **Interquartile
   Range**, **Mean Cross**, the upper and lower whisker dashes, **Show
   Circles**, and change **Bins** from 200 to 100. Check that no new page or
   console error appeared and that the PC Plot viewer is still present on the
   page.
7. Disable the density overlay.

Expected:
- Before any density toggle, the density style is at its default (Circles).
- With all polylines hidden and the selection empty, enabling the density overlay in the Box Plot style adds ink over the sparse axes-only floor: the settle-gated canvas ink count with density on exceeds the hidden-lines floor by more than 1000 pixels.
- Enabling density, cycling the style through circles, box plot and violin, and driving the per-part box/violin toggles (Show Median, Interquartile Range, Mean Cross, upper/lower whisker dash, Show Circles) and the Bins change raises no page or console error and the PC Plot root element is still in the DOM.

### Scenario 3: Density recalculates on normalization changes and log scale (github-1546)

github-1546 was a stale-overlay defect: density shapes kept the geometry of the
previous scale after normalization or log-scale was toggled. The observable guard
is that the recalculation completes without throwing or hanging the page.

Steps:
1. Note the current page-error and console-error counts.
2. Enable the density overlay (**Show Density** in the viewer's Context
   Panel) and set its style to Box Plot (right-click the plot > Style >
   box plot).
3. Uncheck **Normalize Each Column** (Context Panel > Value) — shared global
   scale — and wait for the redraw.
4. Re-check it (per-column normalization) and wait for the redraw.
5. Repeat the toggle once more (uncheck, then re-check) — this is the
   double-toggle path github-1546 left stale.
6. With density still on, hide all polylines again (clear the selection and
   turn **Show All Lines** off) and measure the canvas ink: check that it
   stays clearly above the axes-only floor from Scenario 2 — the overlay is
   still painted after the double-toggle, not stale-empty. Then turn
   **Show All Lines** back on.
7. With density still on, switch the AGE axis to a logarithmic scale
   (Context Panel > Value > Log Columns, add AGE).
8. Run a lightweight scripted check that the page still responds — it completing
   proves the page is not frozen.
9. Check that no new page or console error appeared, that the responsiveness
   check completed, and that the PC Plot viewer element is still present on
   the page.
10. Revert: set AGE back to a linear scale (remove it from Log Columns) and
    disable the density overlay.

Expected:
- github-1546 regression guard: with density on, double-toggling normalization and then switching AGE to a log scale raises no page or console error, the responsiveness check completes (the overlay recalculated and did not freeze), and the PC Plot viewer element is still present on the page.
- After the normalization double-toggle with density still on, measured with all polylines hidden, the settle-gated canvas ink count stays more than 1000 pixels above the hidden-lines floor recorded in Scenario 2 — the density overlay is still painted, not stale-empty.

## Automation notes

- Property drives: the axes are assigned via `pc.props.columnNames = ['AGE',
  'HEIGHT', 'WEIGHT']`; the vertical scale is `pc.props.normalizeEachColumn`
  (false = Global, true = Normalized); the density overlay is
  `pc.props.showDensity` with style `pc.props.densityStyle` (`'circles'`,
  `'box plot'`, `'violin plot'`); the per-part toggles are `showMedian`,
  `showInterquartileRange`, `showMeanCross`, `showUpperDash`, `showLowerDash`,
  `showCircles`, and `bins`; the AGE log scale is
  `pc.props.logColumnsColumnNames = ['AGE']` (revert with `[]`). Each set is
  followed by a 300-400 ms wait for the re-render.
- The Scenario 3 responsiveness check is a follow-up `page.evaluate` returning a
  simple value; "viewer element still present on the page" reads as
  `document.body.contains(pc.root)`.
- Every "measure the canvas ink" step is a settle-gated canvas ink count:
  `v.countCanvasPixels(page, 'PC Plot').total`
  (opaque non-white pixels on the viewer's canvas) re-read every 300 ms until
  two consecutive counts differ by fewer than 200 pixels, up to 5 iterations —
  so a delta between two measured states is the setter's effect, not a render
  tail in flight.
- Every measured ink value is written to `console.log` on green runs so the
  fixed thresholds (500-pixel Scenario 1 state delta, 2000-pixel round-trip
  tolerance, 1000-pixel density-over-floor and anti-stale margins) can be
  audited against live numbers. The thresholds are conservative wide margins,
  deliberately far below the deltas a real repaint produces.
- Density ink is measured on a sparse-canvas basis: `pc.props.showAllLines =
  false` with an empty selection (`df.selection.setAll(false)`) leaves only
  axes and labels painted, so the density shapes are the dominant ink. On a
  full canvas the overlay paints over the dense demog polylines and the pixel
  delta reads zero. `showAllLines` is restored to `true` after each
  measurement.
- Both pixel-delta measurements use the Box Plot style: the circles overlay
  hugs the already-painted axis line and adds no pixel-distinguishable ink
  even on the sparse canvas, so the circles style is exercised only by the
  no-error style-cycling floor, not by a delta assert.
- Scenario 2 records the hidden-lines ink floor before enabling the overlay
  and Scenario 3 reuses it for the anti-stale check, so the two scenarios must
  run in order within one session.
- Each per-part box/violin toggle in Scenario 2 redraws to canvas only, so the
  toggles share the scenario's single no-error + DOM-presence check rather
  than a per-toggle assertion.
- The steps' "clearly differ" / "lands close" / "clearly exceeds" margins are
  the fixed thresholds above: more than 500 pixels between the Normalized and
  Global states, within 2000 pixels on the Scenario 1 round-trip, and more
  than 1000 pixels over the hidden-lines floor for the density-on and
  anti-stale checks.
