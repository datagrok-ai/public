---
feature: matrixplot
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - matrix-plot-spec.ts
related_bugs:
  - id: GROK-20439
    status: open
expected_results:
  - anchor: "Default State"
    expectation: "The Matrix plot viewer is present in the view and the default X
      and Y column sets each contain only numerical and datetime demog columns,
      at most 10 per set (a product read of the auto-picked sets)"
  - anchor: "Back Color"
    expectation: "Setting backColor to red reads back the set value and resetting it
      to white reads back white, with no console or page error; the painted
      background is deliberately not measured (open GROK-20439)"
  - anchor: "Title and Description"
    expectation: "With Show Title on and Title 'My Matrix', the title text renders
      on the viewer; the Description text renders inside the viewer element;
      turning Show Title off and clearing the Title removes the rendered title"
  - anchor: "Inner Viewer Look"
    expectation: "Switching Cell Plot Type between Density plot and Scatter plot
      repaints an off-diagonal cell (settle-gated pixel delta, exact restore);
      setting the inner marker size to 10 and the inner bin count to 20
      round-trips the look property with no console or page error and the 16
      inner cells still rendering — the inner size change does not repaint the
      cells on this build, so it is a prop round-trip plus a no-error/liveness
      floor rather than a render claim"
---

# Matrix plot tests

## Purpose

Verifies the Matrix plot's secondary settings surface: the auto-picked
default column sets, the background color property, the title and
description, and the look settings of the inner cell plots. The focused
scenarios in this folder own the primary workflows (column configuration,
cell plot type, axes and layout, scrolling, row source and filtering, cell
inspection); nothing here duplicates them. Settings whose only outcome is
how the canvas is painted are checked for error-free operation and a visible
repaint rather than judged by appearance.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Matrix plot

## Default State

1. Verify the Matrix plot viewer is present in the view.
2. Read the default **X** column set — it contains only numerical and
   datetime demog columns, at most 10.
3. Read the default **Y** column set — the same rule holds.

## Back Color

1. Set **Back Color** to red via the JS API.
2. Read the property back — it holds the red value; no console or page error
   appeared.
3. Reset **Back Color** to white via the JS API and read back white.

## Title and Description

1. Open the viewer settings and set **Show Title** to true.
2. Set **Title** to "My Matrix" — the title text renders on the viewer.
3. Set **Description** to "Test description" — the description text renders
   inside the viewer element.
4. Set **Show Title** to false and clear the **Title** — the rendered title
   disappears (on this build the titlebar always shows the title text, so
   clearing the Title is what removes it).
5. Clear the **Description**.

## Inner Viewer Look

1. Measure a settle-gated pixel reading of an off-diagonal cell in the default
   Density plot state.
2. Set **Cell Plot Type** to Scatter plot and verify the off-diagonal cell's
   reading changed — the reliable inner-cell render signal.
3. Set the inner viewer's marker size to 10 via the JS API on the look object;
   verify the property round-trips, no console or page error appears, and the
   16 inner cells still render. The size change does not repaint the cells on
   this build, so it is a prop round-trip plus a no-error/liveness floor, not a
   render claim (recorded gap).
4. Set **Cell Plot Type** back to Density plot and verify the off-diagonal
   cell's reading changed again (exact restore of the Density reading).
5. Set the inner viewer's bin count to 20 via the JS API on the look object;
   verify the property round-trips, no console or page error appears, and the
   16 inner cells still render (same no-repaint floor as the marker size).

## Automation notes

- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Matrix plot')`; on demog the
  auto-picked sets are AGE, HEIGHT, WEIGHT, STARTED (numeric + datetime),
  read as `mp.props.xColumnNames` / `yColumnNames` and checked against the
  dataframe's column types, never against a hard-coded list.
- Back Color is a property-echo check by design: the painted background is
  NOT measured because of the open GROK-20439 (back color not applied to the
  drawing); the step therefore asserts the prop round-trip plus a no-error
  floor. When GROK-20439 is fixed, promote this to a pixel assert.
- The title renders in the viewer's title bar and the description inside the
  viewer element — both are DOM text reads scoped to
  `[name="viewer-Matrix-plot"]` and its panel.
- The reliable inner-cell render signal is the **Cell Plot Type** switch: an
  off-diagonal cell's settle-gated non-white pixel count over
  `[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer` (row-major
  index 1) moves between Density plot (~1609) and Scatter plot (~3409) and
  restores exactly on the way back — re-read until two consecutive readings
  agree. The inner look object (`mp.props.innerViewerLook` — marker size for
  the scatter state, bin count for the density state) is not editable through
  a settings-panel row; live recon on this build found that assigning
  `markerSize` / `binCount` round-trips the property and keeps the 16 cells
  rendering but does NOT repaint them. So the inner-size steps assert a prop
  round-trip plus a no-error/liveness floor with the no-repaint gap recorded —
  never a vacuous prop echo presented as a render check.

---
{
  "order": 16,
  "datasets": ["System:DemoFiles/demog.csv"]
}
