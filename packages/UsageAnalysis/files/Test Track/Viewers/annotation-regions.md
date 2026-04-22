# Annotation regions

Visual overlays on Scatter Plot and Line Chart viewers that highlight, select, and analyze specific data areas.

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

---
{
  "order": 29,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
