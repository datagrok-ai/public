---
feature: pcplot
realizes_atlas:
  - pcplot.cp.setup-columns-color-filter
realizes:
  - viewers.pc-plot
priority: p0
target_layer: playwright
coverage_type: smoke
realized_as:
  - pcplot-setup-color-filter-spec.ts
related_bugs:
  - id: GROK-18000
    status: fixed
  - id: GROK-17754
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: "The viewer RENDERS one axis slider per selected column — the DOM
      axis-slider names equal AGE, HEIGHT, WEIGHT in order"
  - anchor: "Scenario 2 Step 4"
    expectation: "df.filter.trueCount drops below the full row count after dragging
      the AGE axis range-slider (its max-handle DOM element) to a narrower window"
  - anchor: "Scenario 2 Step 7"
    expectation: "df.filter.trueCount restores to the full row count after Reset
      View (actuated via the canvas context menu)"
  - anchor: "Scenario 3 Step 6"
    expectation: "No browser console errors after setting the secondary options in
      sequence — adding then removing a column (GROK-18000), with the rendered
      axis-slider set growing 3 → 4 (STARTED, a valid DateTime axis) on add and
      returning to 3 on remove; setting the color column to HEIGHT and switching
      the coloring type to categorical, numerical, then none (GROK-17754),
      toggling Show Filtered Out Lines, and setting AGE to a logarithmic scale;
      the plot keeps updating without a manual refresh"
  - anchor: "Scenario 4 Step 4"
    expectation: "The legend element inside the viewer lists exactly the RACE
      column's categories (Asian, Black, Caucasian, Other) once RACE is set as
      the color column, is absent while no color column is set, and disappears
      again when the color column is cleared"
---

# PC Plot — Setup, Column Selection, Color, In-Chart Range Filter, Log Scale

## Setup

1. Close all open views (`grok.shell.closeAll()` or close all tabs).
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load (full row count, read dynamically).
3. Add a PC Plot viewer to the current table view via the toolbar (Add viewer > PC Plot).
4. Acquire a JS handle to the viewer:
   `const v = grok.shell.tv.viewers.find(v => v.type === 'pc_plot');`

## Scenarios

### Scenario 1: Column setup and axis count assertion

Steps:
1. In the Context Panel > Value > Column Names, clear all selected columns and select AGE, HEIGHT, and WEIGHT (the three numeric columns from demog).
2. Wait for the PC Plot to render (the canvas repaints; at minimum wait for one animation frame via `page.waitForTimeout(300)` or an `evaluate` that resolves without error).
3. Via JS API, read the viewer's column list:
   `const cols = v.getOptions().look.columnNames; // array of active column names`
4. Assert the column list has exactly 3 entries matching AGE, HEIGHT, WEIGHT.

Expected:
- The viewer's columnNames list length is 3 (one axis per selected column: AGE, HEIGHT, WEIGHT)

**Actuation note:** the spec sets the columns by assigning `pc.props.columnNames`
rather than through the Context Panel > Value > Column Names control — that
Select-columns list is canvas-rendered and not scriptable headless. The read-back
is still the RENDERED DOM (one `axis-slider-<col>` element per
column), so a broken re-render fails instead of echoing the prop.

### Scenario 2: In-chart range-filter with Reset View round-trip (PRIMARY SIGNAL)

Steps:
1. Record the initial full row count:
   `const fullCount = grok.shell.tv.dataFrame.filter.trueCount;`
2. On the AGE axis of the PC Plot, drag the top range-slider handle downward to roughly the mid-range — this narrows the AGE window and activates an in-chart filter.
   (`page.mouse.move` + `page.mouse.down` + `page.mouse.move` + `page.mouse.up` on the slider handle DOM element inside the viewer canvas area; alternatively drive via JS `v.setOptions({look: {filter: ...}})` if a stable selector is not available — see Notes.)
3. Wait for the filter to propagate (`page.waitForTimeout(300)` then re-read).
4. Read the filtered count:
   `const filteredCount = grok.shell.tv.dataFrame.filter.trueCount;`
   Assert `filteredCount < fullCount` — the range-slider filter reduced the active row set.
5. Double-click whitespace on the PC Plot canvas to trigger Reset View (fires the RESET_VIEW event).
6. Wait for the filter reset to propagate (`page.waitForTimeout(300)`).
7. Read the restored count:
   `const restoredCount = grok.shell.tv.dataFrame.filter.trueCount;`
   Assert `restoredCount === fullCount` — the filter is fully reset to the original row count.

Expected:
- df.filter.trueCount drops below the full row count after dragging the AGE axis range-slider to a narrower window
- df.filter.trueCount restores to the full row count after Reset View

**Actuation note:** the range-filter is driven by dragging the real
`[name="axis-slider-AGE"] [name="max-handle"]` DOM element (standard mouse events);
Reset View is actuated through the canvas context menu (`Reset View` item), not a
whitespace double-click — the context-menu path is the scriptable-headless
equivalent and fully restores the filter. The full row count is read dynamically
from `df.filter.trueCount`, not hard-coded.

### Scenario 3: Secondary settings — no-error floor (GROK-18000 + GROK-17754 + display toggles)

All of these are secondary settings whose rendered outcome is canvas-only, so they
share a SINGLE consolidated no-error floor rather than a separate assertion each.
Drive them in sequence against a single pageerror/console-error baseline; do not
multiply into near-identical asserts.

Steps:
1. Record the pageerror baseline (subscribe `page.on('pageerror', ...)` and capture the current console-error count).
2. GROK-18000 — change the column selection: add a fourth column (STARTED, a valid DateTime axis — the PC Plot renders `axis-slider-STARTED` and round-trips it in `columnNames`) then remove it, restoring AGE/HEIGHT/WEIGHT. The axes must update immediately with no manual refresh: the rendered axis-slider set grows from 3 to 4 (STARTED present) on add and returns to 3 on remove.
3. GROK-17754 — set the color column to HEIGHT, then switch the coloring type categorical → numerical → none. No leftover legend whitespace / no error.
4. Toggle Show Filtered Out Lines on, then off.
5. Switch the AGE column to a logarithmic scale, then back to linear.
6. Assert once: zero new pageerrors and no new console errors across steps 2–5 combined (the error count equals the baseline). The plot kept updating throughout without a manual refresh.

Expected:
- No browser console errors after exercising the secondary settings in sequence — add then remove a column (GROK-18000), set the color column and switch the coloring type categorical/numerical/none (GROK-17754), toggle Show Filtered Out Lines, and switch a column to logarithmic scale; the plot keeps updating without a manual refresh

### Scenario 4: Categorical coloring renders a legend

Coloring by a categorical column is the one part of the color surface that leaves a
mark in the DOM: the viewer renders a legend element listing the column's own
categories. That list comes from the data, so it is a real check that the color
mapping was applied — unlike reading the color property back.

Steps:
1. With no color column set, look for a legend element inside the viewer root
   (`[name="viewer-PC-Plot"] .d4-legend`) — there should be none.
2. Set the color column to RACE, a categorical column of the demog dataset.
3. Wait for the legend to render.
4. Read the legend's text and assert it lists exactly the four RACE categories:
   Asian, Black, Caucasian, Other.
5. Clear the color column and confirm the legend element is gone again.

Expected:
- The legend element inside the viewer lists exactly the RACE column's categories (Asian, Black, Caucasian, Other) once RACE is set as the color column, is absent while no color column is set, and disappears again when the color column is cleared
