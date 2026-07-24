---
feature: pcplot
realizes_atlas:
  - pcplot.cp.setup-columns-color-filter
  - pcplot.cp.layout-project-persistence
  - pcplot-color-column-legend-coding
  - pcplot-range-filter-cross-viewer
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
  - anchor: "Column setup — select AGE, HEIGHT, WEIGHT"
    expectation: "The viewer RENDERS one axis slider per selected column — the DOM
      axis-slider names equal AGE, HEIGHT, WEIGHT in order"
  - anchor: "In-chart range-filter drop"
    expectation: "df.filter.trueCount drops below the full row count after dragging
      the AGE axis range-slider (its max-handle DOM element) to a narrower window"
  - anchor: "Reset View restore"
    expectation: "df.filter.trueCount restores to the full row count after Reset
      View (actuated via the canvas context menu)"
  - anchor: "GROK-18000 — add then remove a column"
    expectation: "No browser console errors while adding then removing a column
      (GROK-18000): the rendered axis-slider set grows 3 → 4 (STARTED, a valid
      DateTime axis) on add and returns to 3 on remove, with no manual refresh"
  - anchor: "GROK-17754 — color by HEIGHT"
    expectation: "No browser console errors after setting the color column to
      HEIGHT and switching the coloring type to categorical, numerical, then
      none (GROK-17754)"
  - anchor: "Show Filtered Out Lines toggle"
    expectation: "No browser console errors after toggling Show Filtered Out
      Lines on and off; the plot keeps updating without a manual refresh"
  - anchor: "Per-column logarithmic scale for AGE"
    expectation: "No browser console errors after setting AGE to a logarithmic
      scale and back to linear"
  - anchor: "Categorical coloring renders a legend"
    expectation: "The legend element inside the viewer lists exactly the RACE
      column's categories (Asian, Black, Caucasian, Other) once RACE is set as
      the color column, is absent while no color column is set, and disappears
      again when the color column is cleared"
  - anchor: "Numeric coloring gradient drive"
    expectation: "With AGE as the color column, Invert Color Scheme in the
      right-click Color Scheme group flips the invertColorScheme state and a
      second click restores it (menu-to-state round-trip); Edit... in the same
      group opens the 'Color-coding: AGE' dialog, dismissed via its CLOSE
      button. No browser console errors after driving the remaining gradient
      options in sequence — switch the color axis to logarithmic, clamp Color
      Min / Color Max, then reset everything back to no coloring (those options
      repaint the canvas only, so they stay a no-error floor)"
  - anchor: "Legend position cycle and visibility round-trip"
    expectation: "With RACE as the color column, cycling Legend Position through
      Left, Right, Top, Bottom keeps the legend rendered with the same labels;
      Legend Visibility = Never removes the legend element from the DOM and
      Auto restores it with the SAME labels (round-trip); clearing the color
      column removes the legend again"
  - anchor: "Grid conditional color coding surfaces its bins"
    expectation: "Setting conditional color coding on the HEIGHT grid column
      surfaces its bins ('20-150', '150-250') in the plot's DOM legend; switching
      the column to linear color coding drops the DOM legend (a linear/numeric
      gradient is canvas-drawn and has no DOM legend)"
  - anchor: "Layout round-trip"
    expectation: "After saving the layout of the configured viewer (axes AGE,
      HEIGHT, WEIGHT; color column RACE; title 'PC Persistence Probe'), adding
      a Scatter plot, and re-applying the saved layout, the view's viewer set
      equals the SAVED set (a PC Plot is present AND the later-added Scatter
      plot is absent) AND the restored PC Plot carries columnNames
      AGE/HEIGHT/WEIGHT, colorColumnName RACE, and the probe title; the probe
      layout is deleted afterwards even on failure"
  - anchor: "Project save / Close All / reopen"
    expectation: "After saving the view as a project through the ribbon Save
      button, Close All, and reopening the saved project, a PC Plot viewer is
      present in the reopened view AND its columnNames equal AGE/HEIGHT/WEIGHT,
      colorColumnName equals RACE, and the title equals 'PC Persistence Probe'
      (a cross-session round-trip); the probe project is deleted afterwards
      even on failure"
---

# PC Plot — Setup, Column Selection, Color, In-Chart Range Filter, Log Scale

## Purpose

Verifies the PC Plot's day-to-day configuration surface: choosing the axis
columns, filtering rows by dragging an axis range slider (and restoring them
with Reset View), color-coding the lines and the legend behavior that comes
with it, per-column logarithmic scale, and the survival of a configured plot
across a saved layout and a saved project. Settings whose only outcome is how
the canvas is painted are checked as a group for error-free operation rather
than judged by appearance.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a PC Plot viewer to the current table view via the toolbar (Add viewer > PC Plot).

## Scenarios

### Scenario 1: Column setup and axis count assertion

Steps:
1. In the Context Panel > Value > Column Names, clear all selected columns and select AGE, HEIGHT, and WEIGHT (the three numeric columns from demog).
2. Wait for the PC Plot to render (the canvas repaints).
3. Read back the RENDERED axes: the viewer shows one axis slider per active
   column.
4. Verify the rendered axis-slider names equal AGE, HEIGHT, WEIGHT in order.

Expected:
- The viewer renders one axis slider per selected column — the DOM axis-slider names equal AGE, HEIGHT, WEIGHT in order

### Scenario 2: In-chart range-filter with Reset View round-trip

Steps:
1. Record the initial full row count (all rows pass the filter).
2. On the AGE axis of the PC Plot, drag the top range-slider handle downward to roughly the mid-range — this narrows the AGE window and activates an in-chart filter.
3. Wait for the filter to propagate.
4. Read the filtered row count and verify it is lower than the full count — the range-slider filter reduced the active row set.
5. Double-click whitespace on the PC Plot canvas to trigger Reset View.
6. Wait for the filter reset to propagate.
7. Read the restored row count and verify it equals the full count — the filter is fully reset to the original row count.

Expected:
- The filtered row count drops below the full row count after dragging the AGE axis range-slider to a narrower window
- The filtered row count restores to the full row count after Reset View

### Scenario 3: Secondary settings — no-error floor (GROK-18000 + GROK-17754 + display toggles)

These secondary settings only change how the canvas is painted, so they are
exercised in sequence and checked together: the whole sequence must complete
without errors and the plot must keep updating on its own.

Steps:
1. Note the current page-error and console-error counts.
2. GROK-18000 — change the column selection (Context Panel > Value > Column Names): add a fourth column (STARTED, a valid date-time axis) then remove it, restoring AGE/HEIGHT/WEIGHT. The axes must update immediately with no manual refresh: the rendered axis-slider set grows from 3 to 4 (STARTED present) on add and returns to 3 on remove.
3. GROK-17754 — set the color column (Context Panel > Color > Color Column) to HEIGHT, then switch the column's coloring type categorical → numerical → none. No leftover legend whitespace / no error.
4. Toggle **Show Filtered Out Lines** (Context Panel > Data) on, then off.
5. Switch the AGE column to a logarithmic scale (Context Panel > Value > Log Columns), then back to linear.
6. Check once: no new page or console error appeared across steps 2–5 combined. The plot kept updating throughout without a manual refresh.

Expected:
- No browser console errors after exercising the secondary settings in sequence — add then remove a column (GROK-18000), set the color column and switch the coloring type categorical/numerical/none (GROK-17754), toggle Show Filtered Out Lines, and switch a column to logarithmic scale; the plot keeps updating without a manual refresh

### Scenario 4: Categorical coloring renders a legend

Coloring by a categorical column is the one part of the color surface that leaves a
mark in the DOM: the viewer renders a legend element listing the column's own
categories. That list comes from the data, so it is a real check that the color
mapping was applied — unlike reading the color property back.

Steps:
1. With no color column set, look for a legend element inside the viewer —
   there should be none.
2. Set the color column (Context Panel > Color > Color Column) to RACE, a categorical column of the demog dataset.
3. Wait for the legend to render.
4. Read the legend's text and assert it lists exactly the four RACE categories:
   Asian, Black, Caucasian, Other.
5. Clear the color column and confirm the legend element is gone again.

Expected:
- The legend element inside the viewer lists exactly the RACE column's categories (Asian, Black, Caucasian, Other) once RACE is set as the color column, is absent while no color column is set, and disappears again when the color column is cleared

### Scenario 5: Numeric coloring gradient drive with the context-menu Color Scheme group

With a numeric color column the right-click menu carries a **Color Scheme**
group — the menu-reachable part of the gradient surface: its **Invert Color
Scheme** item and its **Edit...** dialog. The remaining gradient options (log
axis, min/max clamps) only repaint the canvas and stay a no-error floor.

Steps:
1. Note the current page-error and console-error counts.
2. Open Context Panel > Color > Color Column and set it to AGE — lines are color-coded by age.
3. Right-click the plot > **Color Scheme** > **Invert Color Scheme** — the gradient direction reverses.
4. Right-click > **Color Scheme** > **Invert Color Scheme** again — the gradient is back to the original direction.
5. Right-click > **Color Scheme** > **Edit...** — the "Color-coding: AGE" dialog opens; close it with its CLOSE button.
6. Change **Color Axis Type** (Context Panel > Color) to logarithmic — the color mapping changes.
7. Set custom **Color Min** (30) and **Color Max** (60) values — the color range narrows.
8. Reset: clear Color Min / Color Max, set Color Axis Type back to linear, clear the color column.
9. Check once: no new page or console error appeared across steps 2–8.

Expected:
- Invert Color Scheme in the context-menu Color Scheme group flips the invertColorScheme state, and a second click restores it (menu-to-state round-trip)
- Edit... in the same group opens the "Color-coding: AGE" dialog, dismissed via its CLOSE button
- No browser console errors after driving the remaining gradient options in sequence (log axis, min/max clamp, full reset)

### Scenario 6: Legend position cycle and visibility round-trip

Steps:
1. Set the color column (Context Panel > Color > Color Column) to RACE and wait for the legend to render; record its labels.
2. Change **Legend Position** (Context Panel) to Left, Right, Top, Bottom — the legend stays rendered with the same labels.
3. Set **Legend Visibility** to Never — the legend disappears.
4. Set **Legend Visibility** to Auto — the legend reappears with the SAME labels (round-trip).
5. Reset Legend Position to Auto and clear the color column — the legend is gone again.

Expected:
- Cycling Legend Position keeps the legend rendered with the same labels; Legend Visibility Never removes the legend element and Auto restores it with the same labels; clearing the color column removes the legend

### Scenario 7: Color coding from the grid column

Conditional color coding renders a DOM legend listing its bins, while a
linear/numeric gradient has no DOM legend — that contrast is the readable signal
that the plot picked up the grid column's color-coding change.

Steps:
1. Set the color column (Context Panel > Color > Color Column) to HEIGHT.
2. In the grid, set conditional color coding on the HEIGHT column with the bins
   '20-150' and '150-250' — the PC Plot legend surfaces the conditional bins.
3. In the grid, switch HEIGHT to linear color coding — the DOM legend disappears (the linear gradient is canvas-drawn).
4. Clean up: reset the HEIGHT column coloring to linear and clear the color column.

Expected:
- Conditional grid color coding surfaces its bins ('20-150', '150-250') in the plot's DOM legend; switching to a linear scheme drops the DOM legend; after cleanup no legend remains

### Scenario 8: Layout and project persistence of the configured viewer

This scenario checks that the settings configured above survive a layout
re-apply and a project reopen. It re-configures a known state first so the
round-trips have concrete values to restore.

Steps:
1. Configure the known state: axes AGE, HEIGHT, WEIGHT; color column RACE; title "PC Persistence Probe".
2. Layout round-trip: save the view layout, add a Scatter plot viewer, then re-apply the saved layout.
3. Verify the view's viewer set equals the SAVED set: a PC Plot viewer is present AND the later-added Scatter plot is absent.
4. Verify the restored PC Plot kept its configuration: axis columns AGE, HEIGHT, WEIGHT; color column RACE; title "PC Persistence Probe".
5. Delete the probe layout.
6. Project round-trip: save the view as a project through the ribbon Save button, dismiss the Share dialog, Close All, then reopen the saved project.
7. Verify a PC Plot viewer is present in the reopened view AND its configuration (axis columns, color column, title) equals the persisted values — a cross-session round-trip.
8. Delete the probe project.
9. Clean up: clear the color column and the title.

Expected:
- Re-applying the saved layout restores the SAVED viewer set (PC Plot present, Scatter plot absent) with the configured columnNames/colorColumnName/title
- Reopening the saved project restores a PC Plot with the same columnNames/colorColumnName/title
- The probe layout and project are deleted even when a verification fails — they never leak

## Automation notes

Setup: the viewer handle is
`const v = grok.shell.tv.viewers.find(v => v.type === 'PC Plot');` views are
closed via `grok.shell.closeAll()`. Row counts are read dynamically from
`grok.shell.tv.dataFrame.filter.trueCount`, never hard-coded. Render waits are
`page.waitForTimeout(300)` or an `evaluate` that resolves without error.

Scenario 1: the spec sets the columns by assigning `pc.props.columnNames`
rather than through the Context Panel > Value > Column Names control — that
Select-columns list is canvas-rendered and not scriptable headless. The read-back
is still the RENDERED DOM (one `axis-slider-<col>` element per column, queried as
`[name="viewer-PC-Plot"] [name^="axis-slider-"]`), so a broken re-render fails
instead of echoing the prop.

Scenario 2: the range-filter is driven by dragging the real
`[name="axis-slider-AGE"] [name="max-handle"]` DOM element (standard mouse events);
Reset View is actuated through the canvas context menu (`Reset View` item), not a
whitespace double-click — the context-menu path is the scriptable-headless
equivalent and fully restores the filter. The full row count is read dynamically
from `df.filter.trueCount`, not hard-coded.

Scenario 3: the error baseline subscribes `page.on('pageerror', ...)` and
captures the current console-error count. The GROK-17754 coloring-type switch
stays a JS API drive (`df.col('HEIGHT').meta.colors.setCategorical()` /
`setLinear()`): the right-click Color Scheme group offers no
categorical/numerical type-switch item — only Edit... (a dialog) and, for a
numeric column, Invert Color Scheme — so there is no menu equivalent to
replace it with.

Scenario 5: the invert and Edit... drives go through the right-click Color
Scheme group. That group appears only while a color column is set, and the
label "Color Scheme" also names the scheme picker inside it, so the spec
resolves the group by its menu-group element and scopes the child clicks to
it. The scheme picker itself is a canvas-drawn control with no headless
handle, so picking a scheme swatch is not driven. Color Axis Type and Color
Min / Color Max remain prop drives: the Color Scheme group carries no items
for them (the Linear/Logarithmic radio entries elsewhere in the right-click
menu belong to the auto-generated property submenu, not the curated Color
Scheme group). The Edit... dialog is dismissed via its
`button[name="button-CLOSE"]` button.

Scenario 4: the legend element is queried as
`[name="viewer-PC-Plot"] .d4-legend`.

Scenario 7: conditional color coding is set via the JS API —
`df.col('HEIGHT').meta.colors.setConditional({'20-150': ..., '150-250': ...})`;
linear coloring via `setLinear()`.

Scenario 8: the layout is saved and re-applied via the JS API
(`tv.saveLayout()` / `grok.dapi.layouts.save` / `tv.loadLayout`) — the View >
Layout menu path has no headless handles; the round-trip end-state is the same.
The project is saved through the real ribbon Save button
(`[name="button-Save"]`) because only the UI Save captures the view layout; a
"Share <project>" dialog pops up after a successful save and is dismissed via
its CANCEL button. The probe project name carries a `Date.now()` suffix so
concurrent runs never collide. The probe layout and project are deleted in
`finally` teardowns (`grok.dapi.layouts.delete` / `grok.dapi.projects.delete`),
so they are removed even when an assertion fails.
