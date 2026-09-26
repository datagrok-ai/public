# Annotation regions: drawing and editing regions on viewers

These scenarios check how a user creates annotation regions on every viewer that supports them, and
how an existing region is changed. Scatter plot, density plot and line chart have two axes, so a
region can be a rectangle, a lasso polygon or a pair of formulas. Histogram, box plot and bar chart
have one value axis. On these viewers the drawn rectangle is locked to the value axis and saved as a
two-formula region, and lasso is not offered. Every new region opens the Formula Lines dialog
(PowerPack) right after it is created. A region is edited from its own context menu, and the changes
stay on the viewer after OK. The region title font is a viewer setting. Viewer regions and dataframe
regions can be hidden separately or together.

## Setup

1. Close all views.
2. Open the demog-1000 dataset (`System:DemoFiles/demog-1000.csv`) and wait for the grid to load.

Each scenario adds its own viewer from **Toolbox > Viewers** and configures it through the controls on
the plot, as listed in the `<configuration>` column:

| viewer | configuration |
|---|---|
| Scatter plot | X selector = WEIGHT, Y selector = HEIGHT |
| Density plot | X selector = WEIGHT, Y selector = HEIGHT |
| Line chart | right-click the HEIGHT chart and pick **HEIGHT > Hide other charts**; X selector = AGE (the Y axis then shows avg(HEIGHT)) |
| Histogram | value selector = AGE |
| Box plot | value selector = AGE, category selector = RACE |
| Bar chart (vertical) | split selector = RACE, value selector = AGE, aggregation next to the value selector = avg; right-click the plot, **Orientation > Vertical** |
| Bar chart (horizontal) | the same, but **Orientation > Horizontal** |

## Scenario 1: The Tools menu offers the region items

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click an empty spot of the plot and point at **Tools**.
3. The Tools submenu lists **Show Annotation Regions**, **Draw Annotation Region** and **Formula Lines...**, in this order.
4. The Tools submenu `<lasso in Tools>` a **Lasso Tool** item.
5. Close the menu with Escape.
6. No errors appear in the console.

| viewer | configuration | lasso in Tools |
|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | does not have (Lasso Tool is at the top level of the menu) |
| Density plot | X = WEIGHT, Y = HEIGHT | has |
| Line chart | HEIGHT chart only, X = AGE | has |
| Histogram | value = AGE | does not have |
| Box plot | value = AGE, category = RACE | does not have |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | does not have |
| Bar chart (horizontal) | split = RACE, value = avg(AGE), horizontal | does not have |

## Scenario 2: Draw a rectangle region on a two-axis viewer and edit it

1. Add a `<viewer>` and apply `<configuration>`.
2. The viewer holds no annotation regions.
3. Right-click an empty spot of the plot and pick **Tools > Draw Annotation Region**.
4. A balloon says "Click and drag to draw region. Press ESC to exit drawing mode.", and the viewer is in region drawing mode.
5. Drag from about one third to about two thirds of the plot across, and from the top edge of the plot down to about one quarter of its height.
6. The viewer leaves drawing mode and shows one viewer region, drawn on the plot.
7. The **Formula Lines** dialog opens on the **Viewer** tab with the new region selected. The **Area** section shows **X column** = `<x column>`, **Y column** = `<y column>` and four corner **Points**. The dialog preview uses the same X and Y axes as the viewer.
8. In the **Description** section, type `Band A` into **Title** and click **OK**. The dialog closes, and the viewer shows one region titled "Band A".
9. Point at a spot inside the region that has no marker or bin under it. The tooltip shows "Band A".
10. Right-click that spot. The context menu lists **Edit...**, **Show Annotation Regions** and **Title Font**.
11. Pick **Edit...**. The **Formula Lines** dialog opens with "Band A" selected, and **Title** shows `Band A`.
12. Set **Title** to `Band B`, **Description** to `Edited band`, **Region Color** to `#ff8800` and **Outline Color** to `#003366`. Set the outline **Width** to `3`, drag the **Opacity** slider to about 60, and set the title **Color** in the **Description** section to `#ff0000`.
13. Each change except **Description** shows in the dialog preview right away. The description is not drawn in the preview; it appears in the region's tooltip (step 15).
14. Click **OK**. The viewer shows one region titled "Band B", filled orange with a thick dark-blue outline.
15. Point at the same spot again. The tooltip shows "Band B" and "Edited band".
16. Right-click an empty spot of the plot, pick **Tools > Formula Lines...**, delete the region with its trash button in the list, and click **OK**.
17. The viewer holds no annotation regions again.
18. No errors appear in the console.

| viewer | configuration | x column | y column |
|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | WEIGHT | HEIGHT |
| Density plot | X = WEIGHT, Y = HEIGHT | WEIGHT | HEIGHT |
| Line chart | HEIGHT chart only, X = AGE | AGE | avg(HEIGHT) |

## Scenario 3: Draw an axis-locked region on a one-axis viewer and edit it

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click an empty spot of the plot and pick **Tools > Draw Annotation Region**.
3. The viewer is in region drawing mode. The pointer looks the same as on a two-axis viewer: it does not hint at the axis the region is locked to; the lock shows up while dragging (step 5).
4. Drag a short diagonal stroke from lower left to upper right in the middle of the plot. On a bar chart, draw it beyond the end of the shortest bar instead, so that part of the region has no bar under it.
5. While dragging, the selection rectangle spans `<locked extent>`, and only the `<free axis>` follows the pointer.
6. On release, the viewer shows one viewer region that spans `<locked extent>`.
7. The **Formula Lines** dialog opens with the region selected. In the **Formula** section, **Formula 1** is `<formula column> = ` followed by the value where the stroke started, and **Formula 2** is `<formula column> = ` followed by the value where it ended, so the two values mark the edges of the region, the smaller one first.
8. The **Title** field in the **Description** section is empty: a drawn region gets no title of its own.
9. Click **OK**.
10. Point at a spot inside the region that has no bar, bin or box under it. The tooltip says "(no title)" and the number of rows the region holds.
11. Right-click that spot and pick **Edit...**. The **Formula Lines** dialog opens with the region selected.
12. Set **Title** to `Band B`, **Description** to `Edited band`, **Region Color** to `#ff8800` and **Outline Color** to `#003366`. Set the outline **Width** to `3`, drag the **Opacity** slider to about 60, and set the title **Color** in the **Description** section to `#ff0000`.
13. Each change except **Description** shows in the dialog preview right away. The description is not drawn in the preview; it appears in the region's tooltip (step 15).
14. Click **OK**. The viewer shows one region titled "Band B", filled orange with a thick dark-blue outline.
15. Point at the same spot again. The tooltip shows "Band B" and "Edited band".
16. Delete the region through **Tools > Formula Lines...** (trash button, **OK**). The viewer holds no regions.
17. No errors appear in the console.

| viewer | configuration | locked extent | free axis | formula column |
|---|---|---|---|---|
| Histogram | value = AGE | the full height of the plot | horizontal position | `${AGE}` |
| Box plot | value = AGE, category = RACE | the full width of the plot | vertical position | `${AGE}` |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | the full width of the plot | vertical position | `${AGE}` |
| Bar chart (horizontal) | split = RACE, value = avg(AGE), horizontal | the full height of the plot | horizontal position | `${AGE}` |

## Scenario 4: Rectangle, lasso and formula regions together, and the dialog lists them all

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click an empty spot of the plot and pick **Tools > Draw Annotation Region**. Drag a rectangle in the upper left part of the plot, and click **OK** in the dialog that opens.
3. Right-click an empty spot of the plot and pick `<lasso item>`. The **Lasso Tool** property of the viewer is on.
4. Right-click an empty spot of the plot and pick **Tools > Draw Annotation Region**.
5. Press the mouse button in the lower right part of the plot, drag through four more points that outline a pentagon, and release the button.
6. The viewer shows a second viewer region. The **Formula Lines** dialog opens with it selected, and its **Points** list has more than four points.
7. Click **OK**.
8. Right-click an empty spot and pick `<lasso item>` again. The **Lasso Tool** property is off.
9. Right-click an empty spot of the plot and pick **Tools > Formula Lines...**. The dialog opens on the **Viewer** tab.
10. Click **ADD NEW** and pick **Region - Formula Lines**.
11. The new region appears in the list and in the preview. Its **Formula 1** is `<formula 1>` and its **Formula 2** is `<formula 2>`: both refer to the viewer's current axis columns.
12. Click **OK**. The viewer holds three viewer regions.
13. The formula region now covers the whole plot, so right-click anywhere inside the plot and pick **Edit...** to open the **Formula Lines** dialog again.
14. The list on the left has the columns **Title**, **Formula** and **Show**, and one row per region: the rectangle, the lasso polygon and the formula region.
15. Select each row in turn. The preview draws the selected region with its own shape, fill color and opacity, and the formula region follows its two formulas.
16. Delete all three rows with their trash buttons and click **OK**. The viewer holds no regions.
17. No errors appear in the console.

| viewer | configuration | lasso item | formula 1 | formula 2 |
|---|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | **Lasso Tool** | `${HEIGHT} = ${WEIGHT} + 168.5` | `${HEIGHT} = ${WEIGHT} - 168.5` |
| Density plot | X = WEIGHT, Y = HEIGHT | **Tools > Lasso Tool** | `${HEIGHT} = ${WEIGHT} + 168.5` | `${HEIGHT} = ${WEIGHT} - 168.5` |
| Line chart | HEIGHT chart only, X = AGE | **Tools > Lasso Tool** | `${avg(HEIGHT)} = ${AGE} + 168.8` | `${avg(HEIGHT)} = ${AGE} - 168.8` |

## Scenario 5: The title font changes region titles only

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT.
2. Right-click an empty spot and pick **Tools > Formula Lines...**. Add **Region - Formula Lines** with **Formula 1** = `${HEIGHT} = 180`, **Formula 2** = `${HEIGHT} = 200` and **Title** = `Tall`.
3. Click **ADD NEW**, pick **Line - Vertical**, type `Reference weight` into **Title**, and click **OK**.
4. The plot shows the region title "Tall" and the line label "Reference weight", both in the default 10 px font.
5. Point at an empty spot inside the "Tall" region and right-click it. In the context menu, click the plus button next to the **Title Font** size four times.
6. Close the menu. The **Annotation Font** property of the viewer is 14 px, and the "Tall" title is drawn larger.
7. The **Formula Font** property is still 10 px, and the "Reference weight" label is unchanged.
8. Right-click the same spot inside "Tall" and click the minus button next to the **Title Font** size four times. The **Annotation Font** property is 10 px again.
9. Delete the region and the line through **Tools > Formula Lines...** and click **OK**. The viewer holds no regions and no formula lines.
10. No errors appear in the console.

## Scenario 6: Hide viewer regions and dataframe regions separately and together

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT.
2. Right-click an empty spot and pick **Tools > Formula Lines...**. Switch to the **DataFrame** tab and add **Region - Formula Lines** with **Formula 1** = `${WEIGHT} = 80`, **Formula 2** = `${WEIGHT} = 100` and **Title** = `Medium weight`. Click **OK**.
3. Open **Tools > Formula Lines...** again. On the **Viewer** tab, add **Region - Formula Lines** with **Formula 1** = `${HEIGHT} = 185`, **Formula 2** = `${HEIGHT} = 205` and **Title** = `Tall`. Click **OK**.
4. The viewer holds one viewer region and one dataframe region, and shows both "Tall" and "Medium weight".
5. Right-click an empty spot and pick **Properties... > Annotations > Show Viewer Annotation Regions**. The "Tall" region disappears, and "Medium weight" stays.
6. Right-click an empty spot and pick **Properties... > Annotations > Show Dataframe Annotation Regions**. "Medium weight" disappears too, and no region is shown.
7. Both regions are still stored: the viewer holds one viewer region and one dataframe region.
8. Right-click an empty spot and pick **Tools > Show Annotation Regions**. Both **Show Viewer Annotation Regions** and **Show Dataframe Annotation Regions** are on, and both regions are shown.
9. Pick **Tools > Show Annotation Regions** again. Both properties are off, and no region is shown.
10. Pick **Tools > Show Annotation Regions** once more. Both regions are shown.
11. Delete both regions through **Tools > Formula Lines...** on the **Viewer** and **DataFrame** tabs, and click **OK**. The viewer holds no regions.
12. No errors appear in the console.

## Scenario 7: A multi-axis line chart offers no region drawing

1. Add a Line chart and set X to AGE with the X selector. Keep the three default Y columns AGE, HEIGHT and WEIGHT.
2. Right-click an empty spot of a chart and point at **Tools**. The submenu lists **Draw Annotation Region**.
3. Close the menu. Right-click the plot and pick **Controls > Multi Axis**. The line chart switches to multi-axis mode.
4. Right-click an empty spot of the plot and point at **Tools**. The submenu has neither **Draw Annotation Region** nor **Show Annotation Regions**.
5. Close the menu. Right-click the plot and pick **Controls > Multi Axis** again. The line chart shows separate charts, and **Tools > Draw Annotation Region** is back.
6. No errors appear in the console.
