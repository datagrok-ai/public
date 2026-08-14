---
feature: 3d-scatter-plot
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [viewers.scatter-plot3d]
realized_as:
  - 3d-scatter-plot-spec.ts
related_bugs: []
---

# 3D Scatter Plot tests — Playwright

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **3D Scatter Plot** from **Toolbox > Viewers**

Everything below is driven through the UI: the Viewers toolbox, the on-viewer column selectors, the Context Panel property grid, the viewer context menu, and real mouse and keyboard input on the plot. The JS API is used only to read back what the UI does not show as text (current row, selected-row count, filtered-row
count).

## Add the viewer

1. Click the **3D Scatter Plot** icon in **Toolbox > Viewers**
2. The on-viewer selectors read **X: AGE**, **Y: HEIGHT**, **Z: WEIGHT** — lfthe axes are assigned automatically

## Axis column assignment

1. Set **X** to *WEIGHT* and **Z** to *AGE* using the on-viewer selectors
2. The selectors show the new columns and the scene repaints
3. Set them back to *AGE* and *WEIGHT*

## Color coding — categorical

1. Set **Color** to *SEX* with the on-viewer selector
2. A legend appears listing exactly *F* and *M*, and the scene repaints

## Color coding — numerical

1. Set **Color** to *AGE*
2. The legend switches away from the categorical items and the scene repaints

## Marker type

1. Go to **Context Panel > Marker**
2. Set **Marker Type** to *box*, then *sphere*, then *cylinder*
3. Each shape paints a visibly different scene

## Marker opacity

1. Set **Marker Opacity** to *25* — the scene repaints
2. Set it back to *100*

## Axes

1. Go to **Context Panel > Axes**
2. Uncheck **Show Axes** — the scene repaints without axes
3. Check it again
4. Set **X Axis Type** to *logarithmic*, then back to *linear* — the scene
   repaints each time

## Camera

1. Drag across the plot — the scene rotates
2. Right-click the plot and choose **Reset View** — the menu closes and the
   camera returns
3. Scroll the wheel up over the plot — the scene zooms in
4. Scroll the wheel down — the scene zooms back out

## Clicking points

1. Click a point — the corresponding row becomes current in the table
2. **Shift+click** a point — it is added to the selection

## Filtered out points

1. Filter the table to `SEX = F`
2. Go to **Context Panel > Misc** and enable **Show Filtered Out Points** —
   the filtered-away rows come back as ghost markers and the scene repaints
3. Disable it and reset the filter

## Mouse-over row group highlight

1. Add a **Bar Chart** to the view
2. Hover a bar — the matching points highlight in the 3D plot (the scene repaints)
3. Uncheck **Show Mouse Over Row Group** in **Context Panel > Misc** and confirm
   the setting is off
4. Re-check it and close the Bar Chart

## Legend

1. Go to **Context Panel > Legend**, set **Legend Visibility** to *Always*
2. Set **Legend Position** to *Left* — the legend moves and the viewer repaints
3. Set it back to *Auto*

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`3d-scatter-plot-spec.ts`. Kept verbatim so no scenario is lost.

### Axis types

> Manual

1. Open properties, set X Axis Type to logarithmic
2. Set Y Axis Type to logarithmic
3. Set Z Axis Type to logarithmic
4. Set X Axis Type back to linear
5. Set Y Axis Type back to linear
6. Set Z Axis Type back to linear

### Size coding

> Manual

1. Open properties, set Size to WEIGHT
2. Set Size to AGE
3. Clear Size

### Labels

> Manual

1. Open properties, set Label to SEX
2. Clear Label

### Marker opacity and rotation

> Manual

1. Open properties, set Marker Opacity to 20
2. Set Marker Opacity to 100
3. Set Marker Opacity back to 69
4. Enable Marker Random Rotation
5. Disable Marker Random Rotation

### Axes visibility and grid lines

> Manual

1. Open properties, disable Show Axes
2. Enable Show Axes
3. Disable Show Vertical Grid Lines
4. Disable Show Horizontal Grid Lines
5. Enable Show Vertical Grid Lines
6. Enable Show Horizontal Grid Lines

### Background and colors

> Manual

Note: color properties set via JS API (UI color picker is not automatable).
1. Set Back Color to black
2. Set Axis Line Color to white
3. Restore Back Color to white and Axis Line Color to default

### Dynamic camera movement

> Manual

1. Open properties, enable Dynamic Camera Movement
2. Disable Dynamic Camera Movement

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
