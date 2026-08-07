---
feature: pivottable
realizes_atlas:
  - pivottable.cp.inner-grid-look-viewers
  - pivottable.int.viewer-columns-per-pivot-category
realizes:
  - viewers.pivot-table
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-16299
    status: fixed
  - id: GROK-15004
    status: fixed
realized_as:
  - pivottable-inner-grid-look-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      After Linear color coding is applied via the inner grid's column-header
      context menu, the value column's gridLook entry reports colorCodingType
      "Linear" and isColorCoded true, with no console error.
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      After Grid > Hide, the value column's gridLook entry reports visible false
      (the GROK-16299 guard — the column is hidden).
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      Clicking the select-all checkbox (the DOM checkbox next to the search
      field) in the Order or Hide Columns dialog restores the hidden column:
      after closing the dialog the value column's gridLook entry reports
      visible true again, with no console error.
  - anchor: "Scenario 4 Step 3"
    expectation: >-
      With one pivot column (SEVERITY), the viewer picker adds one viewer column
      — the Aggregate-row chip count grows by one — with no console error
      (GROK-15004 guard).
  - anchor: "Scenario 4 Step 6"
    expectation: >-
      With two pivot columns (SEVERITY and SEX), the viewer picker runs with no
      console error; the per-category viewer-column count is canvas-only.
  - anchor: "Scenario 6 Step 9"
    expectation: >-
      After re-applying the saved layout, the pivot is a member of the current
      view, its docked inline title reads "Pivot Overview", and the value
      column's gridLook colorCodingType is "Linear" again (title and coloring
      both restored from the layout).
  - anchor: "Scenario 6 Step 10"
    expectation: >-
      Driving the ribbon Save opens the "Save project" dialog with no console
      error (the project-persistence entry point is exercised).
---

# Pivot Table — Inner grid look and viewer columns

## Setup

1. Log in to Datagrok.
2. Close all open views.
3. Open the demog dataset via the File Browser (Files > App Files > Demo Files > demog.csv)
   and wait for the table view to appear.
4. Add the Pivot Table viewer from the Toolbox (Viewers section).
5. Wait for the tag-editor header with three rows labeled Group by, Aggregate, and Pivot.
6. Configure a cross-tab on demog: add DIS_POP as the Group by column using the + icon next
   to the Group by row; confirm the Aggregate row shows avg(AGE) and the Pivot row shows SEVERITY
   from the default auto-configuration.

## Scenarios

### Scenario 1: Apply color coding to an inner-grid aggregate column

Steps:
1. Locate the avg(AGE) column header inside the pivot's inner grid.
2. Record the current console error count.
3. Right-click the avg(AGE) column header in the inner grid to open its column-header context menu.
4. In the context menu, navigate to the color-coding section and select a color scheme
   appropriate for a numeric column — for example, Linear from the Color coding submenu.
   Wait for the color to be applied to the column cells.
5. Verify that at least one data cell in the avg(AGE) column shows a background color that
   differs from the default (uncolored) baseline.
6. Read the current console error count and verify it equals the count from Step 2.

### Scenario 3: Column visibility — Hide and restore via inner-grid column-header menu (GROK-16299)

Steps:
1. Verify the avg(AGE) column is visible in the pivot's inner grid.
2. Record the current console error count.
3. Click any avg(AGE) cell to make it the current column (Hide targets the current column),
   then right-click the avg(AGE) column header in the inner grid, open the "Grid" group and
   click Hide. Wait for the column to disappear from the inner grid.
4. Verify that the avg(AGE) column is no longer visible in the inner grid.
5. Read the current console error count and verify it equals the count from Step 2.
6. Right-click any remaining column header in the inner grid and select Order or Hide Columns
   from the context menu. In the dialog, click the select-all checkbox next to the search field
   (it shows all columns). Close the dialog and wait for the avg(AGE) column to reappear.

### Scenario 4: Viewer columns per pivot category — one column vs. none (GROK-15004)

Steps:
1. Ensure the pivot shows exactly one pivot column (SEVERITY) and at least one aggregate column
   (avg(AGE)). Re-configure if needed by removing excess pivot columns.
2. Record the current console error count.
3. Locate the viewer-type selector near the Aggregate row header. Click it and pick any available
   viewer type (e.g. Bar chart or Sparkline). Wait for the viewer columns to appear in the inner
   grid.
4. Verify one new chip appears in the Aggregate row; the per-category viewer columns are
   visible on the canvas only.
5. Read the current console error count and verify it equals the count from Step 2.
6. Add a second pivot column: click the + icon next to the Pivot row and add SEX. Wait for the
   pivot to recompute with two pivot columns (SEVERITY and SEX).
7. Observe the viewer columns in the inner grid and count them.

### Scenario 6: Persistence tail on spgi-100 — inline title, color coding, layout and project

Steps:
1. Close all open views.
2. Open the spgi-100 dataset via the File Browser (Files > App Data > Chem > tests >
   spgi-100.csv) and wait for the table view to appear.
3. Add the Pivot Table viewer. Configure a cross-tab: add Primary Series Names as the Group by
   column and add Chemist as a second Group by; confirm the Aggregate row shows at least one
   aggregate column.
4. Apply color coding to the first aggregate column via its inner-grid column-header context menu
   (same path as Scenario 1; confirm the non-default background color is visible on the cells).
5. Edit the pivot viewer's title inline: click on the title text in the viewer's header area and
   type "Pivot Overview"; confirm the title updates to show "Pivot Overview".
6. Save the current layout via the ribbon Save button in the toolbar, naming it
   "pivot-look-test". Wait for confirmation.
7. Change the Group by columns (remove one via its × chip) and verify the title still reads
   "Pivot Overview".
8. Re-apply the saved layout "pivot-look-test" from the gallery. Wait for the pivot to update.
9. Read the pivot viewer's inline title and verify the background color of a data cell in the
   colored aggregate column is still non-default.
10. Click ribbon Save, confirm the Save project dialog opens, then Cancel.
11. Delete the layout "pivot-look-test" via the layout management panel (cleanup).

## Automation notes

- target_layer rationale: all color-coding, column visibility, and viewer-column steps are
  driven through the real inner-grid header context menu — a Playwright-accessible DOM menu;
  the assertions read the durable look channel `pv.getOptions(true).look.gridLook.columns[]`.
- The key-column reorder check lives in pivot-table-ui.md as a manual scenario — the re-sort
  is display-only with no durable read-back channel, so no spec step drives an inner-grid
  header drag. Scenarios 2 and 5 are intentionally absent; numbering is kept stable for the
  anchors. The restore path is drivable headless: gridLook.columns[].visible flips back to
  true.
- Setup Step 1 uses `loginToDatagrok`; Setup Step 2 and Scenario 6 Step 1 use `closeAllViews`;
  Setup Step 4 and Scenario 6 Step 3 use `addViewerByIcon`.
- Color-coding menu path: confirmed live — for a numeric column the submenu offers Off /
  Conditional / Linear / Linked ("Categorical" is not offered).
- GROK-16299 Hide path: confirmed live — the Hide action sits under Grid > Hide in the inner
  grid's column-header context menu.
- Viewer-columns signal (Scenario 4): the durable read is
  `pv.getOptions(true).look.gridLook.columns[]` plus the Aggregate-row chip-count delta (one
  new chip per viewer pick); the per-category viewer columns render on the canvas only.
  `grok.shell.tableByName('demog').col('SEVERITY').categories.length` gives the number of
  distinct SEVERITY values.
- Inline title (Scenario 6): the title renders in the docked panel's title-bar text
  (`.panel-titlebar-tabhost .panel-titlebar-text`), confirmed live; the verification reads that
  rendered text back after the layout round-trip.
- Layout save and re-apply use the JS-API layout channel; the project-persistence entry point
  is exercised via the real ribbon Save button (the Save project dialog opens, then Cancel).
- Console-error delta checks use a `page.on('console')` collector that records console messages
  of type error (the harmless cloned-iframe noise is ignored); the count is read before and
  after the guarded action, and a positive delta is a FAIL.
- I6 closure (pivottable.int.viewer-columns-per-pivot-category): Scenario 4, Steps 3-7 — one
  pivot column (viewer columns = category count) and two pivot columns (zero viewer columns).
- Won't-fix guard (GROK-16595 / GROK-17179): a Row Source SWITCH resets the color coding —
  this is accepted behaviour and is never asserted green; no step in this scenario switches Row
  Source and then reads back the color.
