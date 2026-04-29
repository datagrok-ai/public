# Annotation regions

Visual overlays on Scatter Plot, Line Chart, Density Plot, Box Plot, Histogram, and Bar Chart viewers
that highlight, select, and analyze specific data areas. Two-dimensional viewers (Scatter, Line,
Density) support rectangle / lasso / formula regions. One-dimensional viewers (Box, Histogram, Bar)
support a single axis-locked rectangle that the platform persists as a value-axis formula region —
lasso is disabled there because there is only one axis to drag along.

All scenarios start with:
1. Close all
2. Open `System:DemoFiles/demog.csv`
3. Add a Scatter Plot by clicking the Scatter Plot icon in the Toolbox > Viewers section

## Region creation

1. Open viewer properties (gear icon) → Tools — ensure **Lasso Tool** is unchecked
2. Right-click the scatter plot canvas → **Tools** → **Draw Annotation Region** — drag to draw a rectangle
3. Expected: rectangle region appears; Formula Lines dialog opens with the region selected; dialog preview uses the same X/Y axes as the source viewer
4. In the dialog editor change: Title, Description, Fill color, Outline color, Outline width, Opacity — expected: all edits are reflected in the preview and on the source viewer immediately
5. Click OK
6. Open viewer properties → Tools → set **Lasso Tool** = true
7. Right-click → **Tools** → **Draw Annotation Region** — drag a free-form polygon
8. Expected: polygon region created; Formula Lines dialog opens with the new region selected
9. Right-click → **Tools** → **Formula Lines...**; click **ADD NEW** → **Region - Formula Lines**
10. Expected: formula region created with two default boundary formulas referencing the viewer's current X/Y columns; appears in the dialog grid and preview

## Visibility controls

1. Ensure the Scatter Plot has at least one viewer region and one dataframe region
2. Right-click the canvas → **Tools** → **Show Viewer Annotation Regions** — expected: only the viewer region hides; dataframe region stays visible
3. Right-click → **Tools** → **Show Dataframe Annotation Regions** — expected: only the dataframe region hides
4. Right-click → **Tools** → **Show Annotation Regions** — expected: both viewer and dataframe regions become visible; global toggle overrides the two individual ones

## Region editing

1. Open Line Chart on demog and add one formula region via right-click → **Tools** → **Formula Lines...** → **ADD NEW** → **Region - Formula Lines**
2. Right-click on that region → **Edit** — expected: Formula Lines dialog opens with the region pre-selected
3. With the dialog open on an existing region, change each of: Title, Title color, Description, Fill color, Outline color, Outline width, Opacity
4. Expected: each edit reflected in the preview immediately; after OK, changes persist on the source viewer
5. Open viewer properties → find **Annotation Font** (Misc / Style section) — change font family and size
6. Expected: region title labels on all regions update immediately; plain formula line labels are unaffected

## Formula Lines dialog

1. Open the Formula Lines dialog on a Scatter Plot with at least one region
2. Add one of each: rectangle region, lasso region, formula region
3. Expected: preview accurately renders each region — correct boundary shape, fill color, opacity; formula regions show correct curves
4. Inspect the grid on the left side — expected: lists all regions with columns Title, Formula, Show

## Line Chart specific

1. Add a Line Chart; in viewer properties set **Y Columns** to more than one numeric column (multi-axis mode)
2. Right-click → **Tools** → try **Draw Annotation Region** — expected: item is absent from the menu
3. In viewer properties set **Y Columns** to a single numeric column
4. Create one rectangle region, one lasso region, one formula region (same steps as Region creation above)
5. Expected: all three regions created successfully; hover highlighting and click selection behave the same as on Scatter Plot
 behave identically to regions on a Scatter Plot (hover highlighting, click selection, editing).

## Density Plot

Density Plot is two-dimensional, so it behaves like Scatter Plot for region creation and visibility.

1. Close all viewers; add a Density Plot via Toolbox → Viewers → Density Plot
2. In viewer properties set X Column = `WEIGHT`, Y Column = `HEIGHT`; ensure **Lasso Tool** is unchecked
3. Right-click the canvas → **Tools** → **Draw Annotation Region**; drag a rectangle in the central area
4. Expected: a rectangle region is added; Formula Lines dialog opens with that region selected
5. Click OK to dismiss
6. In viewer properties set **Lasso Tool** = true
7. Right-click → **Tools** → **Draw Annotation Region**; drag a freeform polygon (release mouse to close)
8. Expected: a polygon region is added; Formula Lines dialog opens
9. Click OK to dismiss
10. Right-click → **Tools** → **Formula Lines...** → click **ADD NEW** → **Region - Formula Lines**
11. Expected: a formula region is added with default boundary formulas referencing `${WEIGHT}` and `${HEIGHT}`
12. Right-click → **Tools** → **Show Annotation Regions** → expected: regions hide; click again → regions show

## Box Plot

Box Plot is one-dimensional (Y axis = value column). Drawing produces a horizontal value-band that the
platform persists as a value-axis formula region (`${valueColumn} = lo`, `${valueColumn} = hi`). Lasso
is disabled.

1. Close all viewers; add a Box Plot via Toolbox → Viewers → Box Plot
2. In viewer properties: Category Column 1 = `RACE`, Value Column = `AGE`
3. Right-click the canvas → **Tools**; expected items in order: **Show Annotation Regions**, **Draw Annotation Region**, **Formula Lines...**
4. Right-click → **Tools** → **Draw Annotation Region**; drag any vertical extent within the chart
5. Expected:
   - selection rectangle locks to full chart width (only Y follows the cursor)
   - upon mouse release, a formula region is created and the Formula Lines dialog opens with it selected
   - the new region has `formula1 = '${AGE} = <lo>'`, `formula2 = '${AGE} = <hi>'`, `header = '${AGE} in [<lo>, <hi>]'`
6. Click OK
7. Confirm Lasso is unavailable: open viewer properties → Tools — expected: **Lasso Tool** is not present
8. Right-click on the **Y axis** (value axis, left side) → expected: an **Annotations** group appears with **Add Line**, **Add Band**, **Add Region** items
9. Click **Annotations** → **Add Line**; expected: a horizontal formula line is created at the Y column's median value (Q2); Formula Lines dialog opens on it
10. Click OK; right-click on Y axis → **Annotations** → **Add Band**; expected: a band between Q1 and Q3 is created; dialog opens
11. Click OK; right-click on Y axis → **Annotations** → **Add Region**; expected: a formula region between Q1 and Q3 is created; dialog opens
12. Click OK; right-click on the **X axis** (category axis, bottom) → expected: NO **Annotations** group appears (X is categorical)

## Histogram

Histogram is one-dimensional (X axis = value column). Drawing produces a vertical value-band persisted
as a value-axis formula region. Lasso is disabled.

1. Close all viewers; add a Histogram via Toolbox → Viewers → Histogram
2. In viewer properties: Value Column = `AGE`
3. Right-click the canvas → **Tools**; expected items in order: **Show Annotation Regions**, **Draw Annotation Region**, **Formula Lines...**
4. Right-click → **Tools** → **Draw Annotation Region**; drag any horizontal extent within the chart
5. Expected:
   - selection rectangle locks to full chart height (only X follows the cursor)
   - upon mouse release, a formula region is created and the Formula Lines dialog opens with it selected
   - the new region has `formula1 = '${AGE} = <lo>'`, `formula2 = '${AGE} = <hi>'`
6. Click OK
7. Right-click on the **X axis** (bottom) → expected: an **Annotations** group with **Add Line**, **Add Band**, **Add Region**
8. Click **Annotations** → **Add Line**; expected: a vertical formula line at Q2; dialog opens; click OK
9. Right-click on the **Y axis** (left, the count axis) → expected: NO **Annotations** group (Y axis is the derived count axis)
10. Right-click → **Tools** → **Show Annotation Regions** → toggles visibility

## Bar Chart

Bar Chart is one-dimensional with an orientation-dependent value axis: vertical bars place the value
on Y, horizontal bars on X. Drawing locks the categorical axis (full chart-box span) and tracks the
value axis only. Lasso is disabled. The persisted region's `header` includes the aggregation type
(e.g. `${avg(AGE)} in [...]`) because bars plot aggregated values.

1. Close all viewers; add a Bar Chart via Toolbox → Viewers → Bar Chart
2. In viewer properties: Split Column = `RACE`, Value Column = `AGE`, Value Aggr Type = `avg`, Orientation = `Vertical`
3. Right-click the canvas → **Tools**; expected items in order: **Show Annotation Regions**, **Draw Annotation Region**, **Formula Lines...**
4. Right-click → **Tools** → **Draw Annotation Region**; drag any vertical extent
5. Expected:
   - selection rectangle locks to full chart width (only Y follows the cursor); cursor is `ns-resize`
   - upon mouse release, a formula region is created with `formula1 = '${AGE} = <lo>'`, `formula2 = '${AGE} = <hi>'`, `header = '${avg(AGE)} in [<lo>, <hi>]'`
6. Click OK
7. Set Orientation = `Horizontal`; right-click → **Tools** → **Draw Annotation Region**; drag any horizontal extent
8. Expected:
   - selection rectangle now locks to full chart height (only X follows the cursor); cursor is `ew-resize`
   - on release, another formula region is created with the same shape (now mapped through the X axis)
9. Set Orientation back to `Vertical`; right-click on the **Y axis** (value axis) → expected: an **Annotations** group with **Add Line**, **Add Band**, **Add Region**
10. Click **Annotations** → **Add Line**; expected: a horizontal formula line at the aggregated column's Q2; dialog opens; click OK
11. Right-click on the **X axis** (categorical axis) → expected: NO **Annotations** group

## Axis context menu — Annotations group

Verifies the shared `addAxisAnnotationsMenu` helper produces a correctly populated Annotations submenu
(Add Line, Add Band, Add Region) on every viewer that has a numerical axis.

1. Close all viewers; add a Histogram, Value Column = `AGE`
2. Right-click the **X axis** → **Annotations** → **Add Line**
3. Expected: the look's `formulaLines` JSON now contains an entry with `formula = '${AGE} = <q2>'` and `orientation = 'Vertical'`; Formula Lines dialog opens at that entry
4. Click OK
5. Right-click X axis → **Annotations** → **Add Band**
6. Expected: a new entry is added with `type = 'band'`, `formula = '${AGE} in (<q1>, <q3>)'`, `orientation = 'Vertical'`; dialog opens
7. Click OK
8. Right-click X axis → **Annotations** → **Add Region**
9. Expected: the look's `annotationRegions` JSON gets a new entry with `type = 'formula'`, `formula1 = '${AGE} = <q1>'`, `formula2 = '${AGE} = <q3>'`, `header = '${AGE} in [<q1>, <q3>]'`; dialog opens
10. Click OK

## PowerPack absent

Verifies the user-visible warning surfaces when PowerPack is not installed (i.e. the
`PowerPack:formulaLinesDialog` function is unavailable).

1. With PowerPack uninstalled or not yet initialized: open a Bar Chart → right-click → **Tools** → **Formula Lines...**
2. Expected: a warning balloon appears with text "Install PowerPack to edit formula lines"
3. Right-click on the value axis → **Annotations** → **Add Line**
4. Expected: a formula line is still added to the look's JSON, AND a warning balloon "Install PowerPack to edit formula lines" appears (the dialog can't open without PowerPack, but the line is created so the look-side state is consistent)

---
{
  "order": 29,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
