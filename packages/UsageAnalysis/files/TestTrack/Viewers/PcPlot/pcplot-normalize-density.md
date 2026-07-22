---
feature: pcplot
realizes_atlas:
  - pcplot.cp.normalize-and-density
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
  - anchor: "Scenario 1 Step 4"
    expectation: "Switching the vertical scale to a shared global scale and back to
      Normalized raises no page or console error and the PC Plot root element is
      still in the DOM"
  - anchor: "Scenario 2 Step 6"
    expectation: "Enabling density, cycling the style through circles, box plot and
      violin, and driving the per-part box/violin toggles (Show Median,
      Interquartile Range, Mean Cross, upper/lower whisker dash, Show Circles) and
      the Bins change raises no page or console error and the PC Plot root element
      is still in the DOM"
  - anchor: "Scenario 3 Step 8"
    expectation: "github-1546 regression guard: with density on, double-toggling
      normalization and then switching AGE to a log scale raises no page or
      console error, a follow-up page.evaluate resolves (the overlay
      recalculated and did not freeze), and the PC Plot root element is still in
      the DOM"
---

# PC Plot — Normalization and Density Overlay

## Purpose

The vertical-scale and density-overlay surfaces of the PC Plot render entirely to
canvas. There is no DOM signal that says "the box plot now looks right". What can
be checked honestly is that each real setter leaves the plot working: no page error
is raised, the page still responds, and the viewer element is still present in the
DOM. That is the assertion this scenario makes, three times, for the three distinct
flows named in the atlas critical path.

## Setup

1. Close all open views.
2. Open the demog dataset `System:DemoFiles/demog.csv` and wait for the table view.
3. Add a PC Plot viewer to the table view (Toolbox or via ribbon).
4. Assign AGE, HEIGHT and WEIGHT as the three axes via
   `pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT']` and wait for the re-render.
5. Collect a baseline of page-error and console-error counts so each scenario can
   assert "no new errors since baseline".

## Scenarios

### Scenario 1: Vertical scale switched to a shared global scale

The default is per-column normalization. Switch to the non-default state first,
then revert (round-trip).

Steps:
1. Record the current page-error and console-error baseline counts.
2. Set the vertical scale to the shared global scale
   (`pc.props.normalizeEachColumn = false`) and wait 400 ms for the re-render.
3. Revert to per-column normalization
   (`pc.props.normalizeEachColumn = true`) and wait 400 ms.
4. Assert that no new page or console error appeared since the baseline and that
   the PC Plot root element (`pc.root`) is still present in the DOM.

Expected:
- Switching the vertical scale to a shared global scale and back to Normalized raises no page or console error and the PC Plot root element is still in the DOM.

### Scenario 2: Density overlay enabled and its style cycled

Steps:
1. Record the current page-error and console-error baseline counts.
2. Enable the density overlay (`pc.props.showDensity = true`), wait 300 ms.
3. Set the density style to Circles (`pc.props.densityStyle = 'circles'`), wait 300 ms.
4. Set the density style to Box Plot (`pc.props.densityStyle = 'box plot'`), wait 300 ms.
5. Set the density style to Violin (`pc.props.densityStyle = 'violin plot'`), wait 300 ms.
6. Switch the style back to Box Plot, then drive each per-part component toggle:
   Show Median (`showMedian`), Interquartile Range (`showInterquartileRange`),
   Mean Cross (`showMeanCross`), upper/lower whisker dash (`showUpperDash`,
   `showLowerDash`), Show Circles (`showCircles`) and the Bins change (`bins`
   200 → 100). Then assert that no new page or console error appeared since the
   baseline and that the PC Plot root element is still present in the DOM.
7. Disable the density overlay (`pc.props.showDensity = false`) and wait 300 ms.

Each per-part toggle redraws to canvas only, so they share the single no-error +
DOM-presence floor of this scenario rather than a per-toggle assertion.

Expected:
- Enabling density, cycling the style through circles, box plot and violin, and driving the per-part box/violin toggles (Show Median, Interquartile Range, Mean Cross, upper/lower whisker dash, Show Circles) and the Bins change raises no page or console error and the PC Plot root element is still in the DOM.

### Scenario 3: Density recalculates on normalization changes and log scale (github-1546)

github-1546 was a stale-overlay defect: density shapes kept the geometry of the
previous scale after normalization or log-scale was toggled. The observable guard
is that the recalculation completes without throwing or hanging the page.

Steps:
1. Record the current page-error and console-error baseline counts.
2. Enable the density overlay and set its style to Box Plot, wait 300 ms.
3. Switch the vertical scale to Global (`pc.props.normalizeEachColumn = false`),
   wait 400 ms.
4. Switch it back to Normalized (`pc.props.normalizeEachColumn = true`), wait 400 ms.
5. Repeat the toggle once more (Global then Normalized) — this is the double-toggle
   path github-1546 left stale.
6. With density still on, switch the AGE axis to a logarithmic scale
   (`pc.props.logColumnsColumnNames = ['AGE']`), wait 400 ms.
7. Run a follow-up `page.evaluate` that returns a simple value (e.g. `1 + 1`).
   It resolving proves the page is not frozen.
8. Assert that no new page or console error appeared since the baseline, that the
   follow-up evaluate resolved, and that the PC Plot root element is still present
   in the DOM.
9. Revert: set AGE back to a linear scale (`pc.props.logColumnsColumnNames = []`),
   disable the density overlay (`pc.props.showDensity = false`), wait 300 ms.

Expected:
- github-1546 regression guard: with density on, double-toggling normalization and then switching AGE to a log scale raises no page or console error, a follow-up page.evaluate resolves (the overlay recalculated and did not freeze), and the PC Plot root element is still in the DOM.
