---
feature: correlationplot
realizes_atlas:
  - correlationplot.cp.cell-interactions
  - correlationplot.int.ignore-double-click
realizes:
  - viewers.correlation-plot
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20125
    status: fixed
  - id: GROK-19482
    status: fixed
  - id: GROK-19054
    status: fixed
  - id: GROK-19053
    status: fixed
realized_as:
  - correlation-plot-cell-interactions-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      The d4-correlation-plot-corr-cell-click event fires with column1.name and
      column2.name equal to the clicked pair (e.g. HEIGHT and AGE), and value
      equal to the runtime-computed Pearson coefficient for that pair within
      tolerance 1e-3.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      The context panel shows a columnsCorrelation semantic object whose header
      reads "<col1> vs <col2>" for the clicked pair, and the panel contains an
      expandable Scatter plot pane.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      After double-clicking a correlation cell, exactly one new Scatter Plot
      viewer is added to the view — the viewer count increases by one — with
      xColumnName and yColumnName equal to the cell's column pair.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      After setting ignoreDoubleClick to true and double-clicking a correlation
      cell, the viewer count is unchanged (no new Scatter Plot is added),
      confirming the double-click is suppressed (GROK interaction guard).
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      Hovering a correlation cell shows exactly one .d4-tooltip element in the
      DOM (GROK-19482 duplicate-tooltip guard), whose text contains "<Type> R:
      <value>" where <value> matches the runtime-computed coefficient formatted
      to 3 decimal places within tolerance 1e-3.
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      The .d4-tooltip canvas pixel snapshot for a second cell differs from the
      snapshot for the first cell (GROK-20125 per-cell scatter guard) — the
      embedded scatter plot renders the correct column pair, not a static image.
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      Hovering the pinned row-header name column shows the column-statistics
      tooltip (not a correlation R value tooltip).
  - anchor: "Scenario 4 Step 2"
    expectation: >-
      The right-click context menu on a correlation cell contains Show Pearson R
      as a toggle item, a Tooltip submenu with Edit..., Visible, and
      Properties... as items, and a Columns submenu with X Columns and Y Columns
      entries.
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      After clicking Tooltip > Edit..., a modal dialog titled "Edit Tooltip"
      opens and can be closed without errors (GROK-19054 dialog save guard).
  - anchor: "Scenario 4 Step 6"
    expectation: >-
      After clicking Tooltip > Properties..., no dialog opens; instead, the
      embedded tooltip scatter plot's property panel appears in the context
      panel (the item pushes to the context panel, not a modal).
  - anchor: "Scenario 5 Step 3"
    expectation: >-
      After right-clicking a correlation cell and selecting Open as table, a new
      table view opens whose columns contain the column-pair names and
      correlation values; the values match the corresponding matrix cell values
      within tolerance 1e-3 (GROK-19053 table-content guard).
---

# Correlation Plot — Cell Interactions

## Setup

1. Close all open tables and viewers via the Close All action in the main toolbar.
2. Open the demog dataset (System:DemoFiles/demog.csv).
3. Add a Correlation Plot viewer to the current Table View using the Toolbox
   Viewers section (select the Correlation Plot icon), or via the **Add Viewer**
   menu.
4. Locate the Correlation Plot viewer in the layout and record its bounding
   rectangle for cell-coordinate arithmetic.
5. Calibrate the correlation cell grid geometry: click once on an off-diagonal
   cell and read the column pair reported by the viewer (using the refdoc
   §Cell Geometry formula with pinnedW=104, headerH=60, cellW=40, rowH=20);
   adjust pinnedW or headerH if the reported pair is off by a column or row.
6. Record the current viewer count as the baseline for later double-click
   assertions.

## Scenarios

### Scenario 1: Single-click fires event and updates context panel

Steps:
1. Choose an off-diagonal cell for the HEIGHT (X) vs AGE (Y) pair using the
   calibrated coordinates from Setup Step 5.
2. Record the runtime-computed Pearson coefficient for the HEIGHT/AGE pair
   over the full demog row set as the expected value (tolerance 1e-3).
3. Click the center of the HEIGHT vs AGE cell and confirm the viewer fires
   the cell-click event with column1.name = HEIGHT, column2.name = AGE, and
   value equal to the runtime Pearson coefficient within tolerance 1e-3.
4. Inspect the context panel immediately after the click and confirm it shows
   a columnsCorrelation semantic object for the HEIGHT/AGE pair with an
   expandable Scatter plot pane.

Expected:
- The cell-click event fires with column1.name = HEIGHT, column2.name = AGE,
  and value within 1e-3 of the runtime Pearson coefficient (Step 3).
- The context panel displays a columnsCorrelation object headed "<col1> vs <col2>"
  with a Scatter plot pane present (Step 4).

### Scenario 2: Double-click opens scatter plot; ignoreDoubleClick suppresses it

Steps:
1. Choose an off-diagonal cell for the WEIGHT (X) vs AGE (Y) pair using
   calibrated coordinates.
2. Record the current viewer count in the Table View as the pre-double-click baseline.
3. Double-click the center of the WEIGHT vs AGE cell.
4. Count the viewers in the Table View and confirm the count increased by exactly
   one compared to the baseline, and that the new viewer is a Scatter Plot with
   xColumnName equal to WEIGHT and yColumnName equal to AGE.
5. Close the newly opened Scatter Plot viewer (click its close icon scoped to the
   Scatter Plot container).
6. Confirm the viewer count returns to the baseline.
7. In the Context Panel, set the **Ignore Double Click** property to true.
8. Double-click the same WEIGHT vs AGE cell center again.
9. Confirm the viewer count is unchanged from the post-close baseline (no new
   Scatter Plot opened).
10. In the Context Panel, set the **Ignore Double Click** property back to false.

Expected:
- After double-clicking the cell, exactly one new Scatter Plot viewer appears
  with xColumnName = WEIGHT and yColumnName = AGE (Step 3-4).
- After enabling Ignore Double Click and double-clicking again, the viewer count
  is unchanged (Step 7-9), confirming suppression.

### Scenario 3: Hover tooltip — deduplication and per-cell scatter

Steps:
1. Choose two off-diagonal cells with distinct column pairs:
   Cell A = HEIGHT (X) vs AGE (Y); Cell B = WEIGHT (X) vs STARTED (Y).
2. Move the mouse to the center of Cell A and wait for the tooltip to appear.
3. Confirm exactly one tooltip element is visible in the interface (GROK-19482
   duplicate-tooltip guard). Confirm the tooltip text includes "Pearson R:"
   followed by a value within 1e-3 of the runtime Pearson coefficient for
   HEIGHT vs AGE.
4. Capture a canvas pixel snapshot of the embedded scatter plot inside the
   tooltip for Cell A (settle-gated).
5. Move the mouse to the center of Cell B and wait for the tooltip to update.
6. Capture a canvas pixel snapshot of the embedded scatter plot inside the
   tooltip for Cell B.
7. Assert that the Cell A and Cell B canvas snapshots differ (GROK-20125
   per-cell scatter guard — the embedded scatter changes with the cell pair).
8. Move the mouse to the pinned row-header name column area (left of the first
   value column) and confirm that the tooltip shown is a column-statistics
   tooltip, not a correlation R value tooltip.

Expected:
- Exactly one tooltip visible while hovering Cell A; its text contains
  "Pearson R: <value>" within 1e-3 of the runtime coefficient (Step 3).
- The canvas snapshots for Cell A and Cell B differ (Step 7).
- Hovering the row-header name shows a column-statistics tooltip (Step 8).

### Scenario 4: Right-click context menu items and Tooltip subgroup actions

Steps:
1. Right-click the center of an off-diagonal correlation cell.
2. Confirm the context menu opens. Verify that the following items are present
   in the menu: the **Show Pearson R** toggle item, the **Columns** submenu
   entry containing X Columns and Y Columns, and the **Tooltip** submenu entry.
3. Expand the **Tooltip** submenu (hover or click the Tooltip group item).
4. Confirm the Tooltip submenu contains at minimum: **Edit...**, **Visible**,
   and **Properties...** items.
5. Click **Tooltip > Edit...** and wait for the modal dialog to appear. Confirm
   the dialog title is "Edit Tooltip". Close the dialog via its CANCEL button
   and confirm the dialog closes without errors (GROK-19054 dialog save guard).
6. Re-open the context menu on the same cell and expand the Tooltip submenu again.
7. Click **Tooltip > Properties...** and confirm NO new dialog opens; instead
   confirm the context panel updates to show the embedded tooltip scatter plot's
   property panel (look for Row Source or X column scatter props in the context
   panel).
8. Close the context menu if still open (press Escape or click elsewhere).

Expected:
- The context menu contains Show Pearson R, a Columns submenu with X/Y Columns,
  and a Tooltip submenu with Edit..., Visible, Properties... (Step 2-4).
- Clicking Tooltip > Edit... opens the "Edit Tooltip" modal which closes without
  errors (Step 5).
- Clicking Tooltip > Properties... pushes the tooltip scatter property panel to
  the context panel without opening a dialog (Step 7).

### Scenario 5: Open as table — column pairs and values match matrix

Steps:
1. Record the runtime Pearson coefficient for the HEIGHT vs AGE pair (by reading
   from the viewer's backing DataFrame) as the expected value (tolerance 1e-3).
2. Right-click the center of the HEIGHT vs AGE correlation cell and locate the
   **Open as table** menu item (exact text "Open as table", lowercase).
3. Click **Open as table** and wait for a new table view to open.
4. Inspect the new table: confirm it contains a row whose column-pair fields
   identify HEIGHT and AGE, and whose correlation value matches the expected
   Pearson coefficient within tolerance 1e-3 (GROK-19053 table-content guard).
5. Close the table view produced by Open as table.

Expected:
- A new table view opens after clicking Open as table (Step 3).
- The table contains a row for the HEIGHT/AGE pair with a correlation value
  within 1e-3 of the runtime Pearson coefficient (Step 4).

## Automation notes

- target_layer rationale: all five scenarios require real Playwright page.mouse
  gestures (click, dblclick, move, contextmenu) on the ColumnGrid overlay canvas;
  cell coordinates computed from the viewer bounding box plus refdoc §Cell
  Geometry arithmetic (pinnedW, headerH, cellW, rowH); the d4-correlation-plot-
  corr-cell-click event and the context panel push are the primary assertable
  signals. The Playwright layer is the only one with live mouse input + in-page
  JS API access in a single test.
- Cell coordinate calibration (Setup Step 5): fire one probe click and validate
  the reported pair from the event args; if off by a cell, adjust pinnedW or
  headerH. Do NOT read gridCols or the inner ColumnGrid from JS — it is not
  exposed. Calibrate from event feedback only.
- Real vs synthetic input: specs use real Playwright page.mouse events (CDP
  path). The refdoc §Cell Geometry notes that synthetic dispatchEvent also fired
  the click handler in the 2026-07-31 build, but DO NOT rely on it — the spec
  must use real input so it is isTrusted-consistent across builds.
- Setup Step 2 actuation: open demog.csv via readDataframe and add it to the
  workspace via grok.shell.addTableView.
- Setup Step 4 actuation: locate the Correlation Plot root element via
  [name="viewer-Correlation-plot"] and record the viewer bounding rectangle.
- Setup Step 5 actuation: fire a real Playwright page.mouse.click at the
  computed off-diagonal center and read the column pair from the resulting
  d4-correlation-plot-corr-cell-click event.
- Scenario 1 Step 3 actuation: perform a real Playwright page.mouse.click at
  the cell center and capture the d4-correlation-plot-corr-cell-click event
  arguments (column1.name, column2.name, value).
- Scenario 2 Step 3 and Step 8 actuation: use a real Playwright
  page.mouse.dblclick at the computed cell center.
- Scenario 3 Step 2 actuation: use page.mouse.move to the cell center and wait
  for the .d4-tooltip element to appear in the DOM.
- Pearson reference computation: use the exposed getCorrelation(c1, c2) or
  compute Stats.corr(df.getCol(c1), df.getCol(c2)) at runtime; never hardcode
  the expected value; tolerance 1e-3 (widen only if live recon shows stable WASM
  divergence).
- GROK-19482 guard (Scenario 3 Step 3): assert
  document.querySelectorAll('.d4-tooltip').length === 1 while a tooltip is
  visible; move the mouse away first to dismiss any stale tooltip before hovering
  the new cell.
- GROK-20125 guard (Scenario 3 Step 7): capture canvas snapshots with
  snapshotCanvasColors / diffCanvasColors from ../../helpers/viewers.ts (or
  equivalent settle-gated canvas diff); a non-zero diff confirms the embedded
  scatter changed between cells.
- GROK-19054 guard (Scenario 4 Step 5): Tooltip > Edit... opens the modal titled
  "Edit Tooltip" (DOM confirmed 2026-07-31); click CANCEL to close; assert no
  exception in the console around the save path. The test does NOT need to
  submit the dialog to exercise the guard — opening and closing cleanly is the
  regression signal.
- GROK-19053 guard (Scenario 5): "Open as table" exact menu text is lowercase
  "Open as table" per the refdoc §Context Menu DOM observation; find by
  .d4-menu-item-label text match.
- Context menu navigation: menu items have no name= attributes; locate by
  iterating .d4-menu-item-label elements and matching text; to get a group's
  own label, use the refdoc §Context Menu DOM note (group label comes AFTER
  nested children, so use closest('.d4-menu-item') identity check).
- Scenario 4 Step 2 actuation: right-click via a real Playwright page.mouse.click
  with button 'right' (or context menu event on the overlay canvas); find menu
  items by iterating .d4-menu-item-label elements and matching text.
- Tooltip > Properties... vs Edit... distinction: Properties... pushes the
  tooltip scatter's props to the context panel (no modal); Edit... opens the
  "Edit Tooltip" modal. Assert the correct outcome for each per Scenario 4
  Steps 5 and 7.
- Viewer count assertion (Scenario 2): read
  grok.shell.tv.viewers.filter(v => v.type === 'Scatter plot').length before
  and after the double-click.
- Tooltip read channel: the single reusable .d4-tooltip element hosts the
  correlation tooltip content inline — its textContent carries the '<Type> R:
  <value>' line (probe-verified). An EMPTY .d4-tooltip element means the hover
  missed the cell (the element's idle state), not that the text lives on a
  hidden layer.
- Geometry staleness: re-read the viewer root's bounding rect immediately
  before every hover/click group that follows a viewer add/remove — dock
  relayout shifts the root, so cached Setup-time coordinates miss the 40x20
  cells. After hovering, poll for the '<Type> R:' line in the .d4-tooltip
  textContent instead of using fixed waits.
- Emission discipline for expected_results_coverage: every realized_by MUST be
  a verbatim, contiguous quote copied from the emitted spec (a single line is
  enough — e.g. the exact expect(...) line such as "expect(after).toBe(before)"
  with its real variable names); never paraphrase, never join multiple
  statements with ';', never use '(...)' or '…' ellipses — the mechanical
  verbatim matcher rejects them and Gate E fails with E-EXPECT-COVERAGE-01.
  For multi-line assertions quote only the final expect line.
- GrokML console allowlist: 'Package GrokML is not available...' is a benign
  console warning — add it to the console noise allowlist so it does not fail
  the clean-console assert.
- Column rename prohibition: no step in this scenario renames any column.
