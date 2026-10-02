---
feature: trellisplot
target_layer: manual-only
coverage_type: smoke
---

# Trellis plot tests (manual checklist)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Trellis plot


## Floating viewer after applying layout

1. Close all and open demog
2. Use Ctrl+'-' to zoom out the screen view
3. Add a trellis plot
4. Undock the viewer and drag it to the bottom of the screen
5. Save the layout
6. Use Ctrl+'+' to zoom in to the original screen size
7. Apply the saved layout — verify the viewer is properly positioned and accessible on screen

## Viewer basics

1. Undock the viewer and move it around the screen
2. Dock the viewer in a different location
3. Open the inner viewer type selector at the top of the trellis plot and look through the
   list — the drop-down renders in full, every entry is readable, and after picking a type
   the control panel stays legible (the functional ladder over all inner types is automated).
   With the Curves package deployed the list ends with **Curves**; without it the list is one
   entry shorter, which is expected and not a defect
4. Hover and select elements inside the viewer
5. Open/Close the Context Panel
6. Display the help window

## Context menu (detailed)

1. Right-click a trellis cell — six to eight items are visible, among them **General**,
   **Tooltip**, **Properties...** and (for a Scatter plot or Box plot inner viewer) a group
   header named after the inner type. **General** and **Tooltip** are groups: they open on
   HOVER, not on click, and their submenu slides out to the right. Move the pointer straight
   along the row to the header — wandering off the row on the way closes the submenu again
2. Hover **General**, then click through its items (Clone, Full Screen, Save to Gallery, Save as
   PNG, Embed...) and check that each does what it says — leave **Close** for last, it removes
   the viewer
3. Hover **Tooltip** and walk its items the same way
4. Open the Context Panel
5. Open **Properties...** and verify that property changes are consistent between the Context
   Panel and the context menu

## Create a trellis from another viewer

This section starts from a scatter plot, not from a trellis plot.

1. Close all, open demog, and add a **Scatter plot** with X = AGE, Y = HEIGHT, coloured by SEX
2. Right-click the scatter plot — **Use in Trellis** is not one of the visible top-level items;
   it lives inside the **General** group
3. Hover **General**, moving the pointer straight along the row — the submenu opens to the right
   within about a second, and the rest of the menu stays where it was
4. Click **General | Use in Trellis** — a new Trellis plot appears with Scatter plot as its inner
   viewer type, carrying over the source viewer's X, Y and colour settings; the original scatter
   plot stays open

## Category scrolling by hand

1. Set X to DIS_POP and Y to RACE, then press the **-** icon on each axis (or shrink the
   viewer) until fewer categories are shown than the columns hold — the scroll sliders at the
   bottom and left edges of the chart area now have somewhere to travel
2. Drag the horizontal scroll slider — the visible X categories shift as you drag and the
   labels stay in step with the cells
3. Drag the vertical scroll slider — the visible Y categories shift the same way
4. Hover over the chart area and turn the mouse wheel — the grid scrolls vertically through
   the Y categories
5. Drag both sliders back to their starting ends — the original category window is restored
   and no cell is left blank or half-drawn

## Curves inner viewer (and table switching)

1. Open **Files > Demo > curves.csv**, then go back to the demog view
2. In the Context Panel > **Data**, set **Table** to curves
3. Select **Curves** as the inner viewer — that is the label the picker shows for the multi-curve
   viewer. The entry comes from the Curves package rather than from the platform itself, so on a
   server where that package is not deployed it is simply absent and this whole section is
   skipped
4. Set the X and Y axes on the inner viewer — every cell redraws with the chosen axes and the
   curves are readable
5. Change the number of categories shown on each axis with the **+** / **-** icons — the grid
   grows and shrinks accordingly and the curves stay correctly scaled
6. Move and resize the zoom slider for the X axis, then for the Y axis — the curves rescale in
   every cell and the axis labels follow the new range
7. Click the **Gear** icon and check the inner viewer properties

## Inner viewer color coding

1. Set the inner viewer type to **Pie chart**
2. Click the **Gear** icon on the trellis title bar and switch to the inner viewer's own tab of
   the Context Panel
3. Set the inner viewer's **Marker Color** to RACE
4. Look at a populated cell — the colour coding is applied, one category keeps the same colour
   across all cells, and no cell is left showing the pre-change colours
5. Switch the inner viewer to **Box plot** and repeat — the colour change lands there too

## Cleanup

1. Delete whichever of these this run saved (a partial run may have created only some): the
   layout from the "Floating viewer" section and the gallery entry from **Save to Gallery**
   in the context-menu walk (View > Layout > Open gallery / Browse > Platform > Layouts).
2. Discard the PNG file downloaded by **Save as PNG**.
3. Close all.

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/curves.csv"]
}
