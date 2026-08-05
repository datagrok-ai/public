---
feature: form
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - form-spec.ts
related_bugs: []
---

# Form viewer (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Form** from **Toolbox > Viewers**

The viewer is plain DOM, so every check below reads the text on screen and
compares it with the table's own values.

## Add the viewer

1. Click the **Form** icon in **Toolbox > Viewers**
2. A field appears for every column of the table

## The form shows the current row

1. The values of *USUBJID*, *SEX*, *RACE* and *DIS_POP* match the current row

## Navigation

1. Click the **next row** arrow — the current row moves on by one and the fields follow it
2. Click the **previous row** arrow — both go back

## Following the grid

1. Make another row current in the grid — the fields show that row

## Row selector

1. Click the row selector on the form — exactly that row becomes selected
2. Click it again — the selection clears

## Toolbar visibility

1. Uncheck **Show Navigation** in **Context Panel > Misc** — the toolbar disappears
2. Check it again — the toolbar comes back
3. Uncheck **Show Next Row Arrow** — only the next arrow goes, the previous one stays

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`form-spec.ts`. Kept verbatim so no scenario is lost.

### Keyboard navigation

> Manual

1. Click on the Form viewer to give it focus
2. Press Right arrow — form advances to the next row
3. Press Down arrow — form advances to the next row again
4. Press Left arrow — form goes to the previous row
5. Press Up arrow — form goes to the previous row again
6. Press Space — current row selection toggles

### Sync mode

> Manual

1. Right-click the Form viewer to open the context menu
2. Open the **Track Row** submenu — verify **Current** is checked by default
3. Select **Mouse Over** — hover over different rows in the grid, form updates to the hovered row
4. Select **None** — form stops updating when current row or mouse-over changes
5. Select **Current** — form resumes tracking the current row

### Edit mode

> Manual

1. Click `[name="icon-edit"]` in the Form toolbar — the icon becomes active (checked state)
2. Clear the AGE field and type a new numeric value
3. Verify the corresponding cell in the grid shows the updated value
4. Click `[name="icon-edit"]` again — the icon returns to unchecked state
5. Attempt to type in the AGE field — input should not be accepted

### Column selector

> Manual

1. Click `[name="icon-list"]` in the Form toolbar — a column selection dialog opens
2. Uncheck all columns, then check only AGE and SEX
3. Click OK — the form shows only AGE and SEX fields
4. Click `[name="icon-list"]` again and additionally check RACE
5. Click OK — the form shows AGE, SEX, and RACE

### Build form from column list

> Manual

1. Click `[name="icon-list"]` in the Form toolbar — a column selection dialog opens
2. Toggle several checkboxes (check AGE, HEIGHT, WEIGHT; uncheck the rest)
3. Click OK — form shows only the selected fields
4. Save the layout
5. Close the Form viewer
6. Apply the saved layout
7. Verify the Form viewer is restored with AGE, HEIGHT, and WEIGHT fields

### Filtered data navigation

> Manual

1. Open the filter panel
2. Set a filter on SEX to show only "M" rows
3. Click `[name="icon-chevron-right"]` — form advances only to rows where SEX = M
4. Click `[name="icon-chevron-left"]` — form goes back, skipping non-M rows
5. Remove the filter — navigation covers all rows again

### Context menu

> Manual

1. Right-click the Form viewer
2. Verify **Edit Form...** item is present in the context menu
3. Verify **Select Columns...** item is present
4. Verify **Track Row** submenu is present
5. Open **Track Row** submenu — verify three options: None, Current, Mouse Over

### Column changes reaction

> Manual

1. In the grid, rename the AGE column to AGE_NEW
2. Verify the form field label updates to AGE_NEW
3. Delete a column that is currently shown in the form (e.g., HEIGHT)
4. Verify the form removes the HEIGHT field without error

### Design mode (field dragging)

> Manual

1. Click `[name="icon-object-ungroup"]` in the toolbar — icon becomes checked, fields become draggable
2. Drag the AGE field to a new position on the form canvas
3. Verify AGE stays in the new position
4. Click on a field while in design mode — its properties appear in the Context Panel
5. Click `[name="icon-object-ungroup"]` again to exit design mode — fields are no longer draggable

### Color coding

> Manual

1. In the grid, apply color coding to the AGE column
2. Verify the AGE field in the form picks up the color coding from the grid

### Table switching

> Manual

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
