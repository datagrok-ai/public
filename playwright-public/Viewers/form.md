# Form tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Form

## Row navigation

1. Verify the form shows row 0 values by default
2. Click `[name="icon-chevron-right"]` in the Form toolbar — form updates to row 1
3. Click `[name="icon-chevron-left"]` — form returns to row 0
4. In the grid, click row 5 — form updates to show row 5 values

## Keyboard navigation

1. Click on the Form viewer to give it focus
2. Press Right arrow — form advances to the next row
3. Press Down arrow — form advances to the next row again
4. Press Left arrow — form goes to the previous row
5. Press Up arrow — form goes to the previous row again
6. Press Space — current row selection toggles

## Row selection

1. Click `[name="icon-square"]` in the Form toolbar — current row (row 0) becomes selected
2. Click `[name="icon-square"]` again — row 0 is deselected
3. In the grid, select row 3 — the square icon in the Form toolbar should appear checked
4. Navigate the form to a different row — the square icon updates to reflect whether the new row is selected

## Sync mode

1. Right-click the Form viewer to open the context menu
2. Open the **Track Row** submenu — verify **Current** is checked by default
3. Select **Mouse Over** — hover over different rows in the grid, form updates to the hovered row
4. Select **None** — form stops updating when current row or mouse-over changes
5. Select **Current** — form resumes tracking the current row

## Edit mode

1. Click `[name="icon-edit"]` in the Form toolbar — the icon becomes active (checked state)
2. Clear the AGE field and type a new numeric value
3. Verify the corresponding cell in the grid shows the updated value
4. Click `[name="icon-edit"]` again — the icon returns to unchecked state
5. Attempt to type in the AGE field — input should not be accepted

## Column selector

1. Click `[name="icon-list"]` in the Form toolbar — a column selection dialog opens
2. Uncheck all columns, then check only AGE and SEX
3. Click OK — the form shows only AGE and SEX fields
4. Click `[name="icon-list"]` again and additionally check RACE
5. Click OK — the form shows AGE, SEX, and RACE

## Build form from column list

1. Click `[name="icon-list"]` in the Form toolbar — a column selection dialog opens
2. Toggle several checkboxes (check AGE, HEIGHT, WEIGHT; uncheck the rest)
3. Click OK — form shows only the selected fields
4. Save the layout
5. Close the Form viewer
6. Apply the saved layout
7. Verify the Form viewer is restored with AGE, HEIGHT, and WEIGHT fields

## Toolbar visibility properties

1. Open the Context Panel (F4) with the Form viewer selected
2. Set **Show Navigation** to false — the entire toolbar disappears
3. Set **Show Navigation** back to true — toolbar reappears
4. Set **Show Prev Row Arrow** to false — `[name="icon-chevron-left"]` disappears from toolbar
5. Set **Show Next Row Arrow** to false — `[name="icon-chevron-right"]` disappears
6. Set **Show Row Selector** to false — `[name="icon-square"]` disappears
7. Set **Show Field Editor** to false — `[name="icon-edit"]` disappears
8. Set **Show Design Editor** to false — `[name="icon-object-ungroup"]` disappears
9. Set **Show Column Selector** to false — `[name="icon-list"]` disappears
10. Set **Show Save File** to false — `[name="icon-download"]` disappears
11. Set **Show Open File** to false — `[name="icon-folder-open"]` disappears
12. Restore all visibility properties to true

## Filtered data navigation

1. Open the filter panel
2. Set a filter on SEX to show only "M" rows
3. Click `[name="icon-chevron-right"]` — form advances only to rows where SEX = M
4. Click `[name="icon-chevron-left"]` — form goes back, skipping non-M rows
5. Remove the filter — navigation covers all rows again

## Context menu

1. Right-click the Form viewer
2. Verify **Edit Form...** item is present in the context menu
3. Verify **Select Columns...** item is present
4. Verify **Track Row** submenu is present
5. Open **Track Row** submenu — verify three options: None, Current, Mouse Over

## Column changes reaction

1. In the grid, rename the AGE column to AGE_NEW
2. Verify the form field label updates to AGE_NEW
3. Delete a column that is currently shown in the form (e.g., HEIGHT)
4. Verify the form removes the HEIGHT field without error

## Design mode (field dragging)

1. Click `[name="icon-object-ungroup"]` in the toolbar — icon becomes checked, fields become draggable
2. Drag the AGE field to a new position on the form canvas
3. Verify AGE stays in the new position
4. Click on a field while in design mode — its properties appear in the Context Panel
5. Click `[name="icon-object-ungroup"]` again to exit design mode — fields are no longer draggable

## Color coding

1. In the grid, apply color coding to the AGE column
2. Verify the AGE field in the form picks up the color coding from the grid

## Table switching

1. Open SPGI
2. Add a Form viewer
3. Open viewer settings (F4) → Data → Table — switch to demog
4. Verify the form fields rebind to demog columns
5. Switch Table back to SPGI — form fields rebind to SPGI columns
6. Close All

---
{
  "order": 13,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
