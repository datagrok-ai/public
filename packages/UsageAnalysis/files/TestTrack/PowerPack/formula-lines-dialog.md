# Formula lines: the Formula Lines dialog

These scenarios check the Formula Lines dialog (PowerPack) as the place where a user adds and edits
formula lines and bands. A formula whose left side is not a single column is rejected, and the item
still opens normally the next time. Lines that share a formula but have different ranges are all
drawn. The Show checkboxes hide and bring back items, and color and style changes are kept. A datetime
axis accepts lines. The preview shows the viewer's own axes. A dataframe line appears on every viewer
with the same columns, and a line chart keeps its lines through a saved layout.

## Setup

1. Close all views.
2. Open the demog-1000 dataset (`System:DemoFiles/demog-1000.csv`) and wait for the grid to load.

Viewers are added from **Toolbox > Viewers**. Their columns are set with the selectors on the plot.
The dialog is opened by right-clicking an empty spot of the plot and picking **Tools > Formula
Lines...**. In a formula field, typing `$` opens the column list, and the column is picked from it.

## Scenario 1: A formula with a column on both sides is rejected and does not break the line

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT, and open the Formula Lines dialog.
2. Click **ADD NEW** and pick **Line**. The formula editor holds `${HEIGHT} = ${WEIGHT}`.
3. Change the formula to `${HEIGHT} - ${WEIGHT} = 1`. The formula editor turns invalid, **OK** is disabled, and the dialog stays open and editable.
4. Change the formula to `${HEIGHT} = ${WEIGHT} + 1`. The editor is valid again and **OK** is enabled.
5. Click **OK**. The viewer holds exactly one formula line, `${HEIGHT} = ${WEIGHT} + 1`.
6. Open the Formula Lines dialog again. The line is listed and selected, its formula editor shows `${HEIGHT} = ${WEIGHT} + 1`, and the editor is not invalid (GROK-17282).
7. Delete the line with its trash button and click **OK**. The viewer holds no formula lines.
8. No errors appear in the console.

## Scenario 2: Two lines with the same formula and different ranges are both drawn

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT, and open the Formula Lines dialog.
2. Click **ADD NEW**, pick **Line**, type `60` into the **Range** min field and `90` into the max field, and set **Title** to `Light`.
3. Click **ADD NEW**, pick **Line**, type `100` into the **Range** min field and `150` into the max field, and set **Title** to `Heavy`.
4. Click **OK**. The viewer holds two formula lines, both `${HEIGHT} = ${WEIGHT}`, and both are drawn: one for WEIGHT 60 to 90 and one for WEIGHT 100 to 150.
5. Delete both lines through the dialog and click **OK**. The viewer holds no formula lines.
6. No errors appear in the console.

## Scenario 3: Unchecked items disappear and come back when checked again

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT, and open the Formula Lines dialog.
2. Add two **Line** items (titles `Light` and `Heavy`) and one **Band - Horizontal** item (title `Band`), and click **OK**. The viewer draws three formula items.
3. Open the dialog again, uncheck **Show** in two of the three rows, and click **OK**. The viewer draws one formula item.
4. Open the dialog again. The two rows are still unchecked, and the third is checked.
5. Check **Show** in both unchecked rows and click **OK**. The viewer draws all three formula items again.
6. Delete all three items through the dialog and click **OK**. The viewer holds no formula lines.
7. No errors appear in the console.

## Scenario 4: The color and style of a line are kept

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT, and open the Formula Lines dialog.
2. Click **ADD NEW** and pick **Line - Horizontal**. The **Constant line** section shows **Column** = HEIGHT and **Value** = 168.5.
3. In the **Format** section, set **Color** to `#ff0000` and **Style** to `dashed`. The preview draws the line red and dashed.
4. Click **OK**. The viewer draws the line `${HEIGHT} = 168.5` red and dashed.
5. Open the dialog again. **Color** is `#ff0000` and **Style** is `dashed`.
6. Delete the line and click **OK**. The viewer holds no formula lines.
7. No errors appear in the console.

## Scenario 5: A line on a datetime axis

1. Add a Scatter plot, set X to STARTED (a datetime column) and Y to HEIGHT, and open the Formula Lines dialog.
2. Click **ADD NEW** and pick **Line - Vertical**.
3. The dialog accepts the datetime column: the **Constant line** section shows **Column** = STARTED, no error appears, and the preview draws a vertical line at the median start date, in 1991.
4. Click **OK**. The scatter plot draws the vertical line at the same date (#2487).
5. Delete the line through the dialog and click **OK**. Set X back to WEIGHT. The viewer holds no formula lines.
6. No errors appear in the console.

## Scenario 6: The dialog preview follows the viewer's current axes

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT, and open the Formula Lines dialog.
2. The preview shows X = WEIGHT and Y = HEIGHT.
3. Click **ADD NEW** and pick **Line - Horizontal**. The preview shows the new line. Click **OK**.
4. Set the scatter plot's X to AGE and Y to WEIGHT.
5. Open the Formula Lines dialog again. The preview shows X = AGE, the viewer's new X axis, and Y = HEIGHT, the column the existing line is built on, so the line is still drawn in the preview. The scatter plot itself, now plotting WEIGHT on Y, does not draw the line (#671).
6. Select the line in the list. The preview highlights the selected line.
7. Delete the line, click **OK**, and set X back to WEIGHT and Y back to HEIGHT.
8. No errors appear in the console.

## Scenario 7: A dataframe line is drawn wherever the axis carries its column

The line is built on WEIGHT, which has no empty values, so a line chart over it renders whatever its
X axis is.

1. Add a Scatter plot with X = AGE and Y = WEIGHT.
2. Add a Line chart, keep only the WEIGHT chart (right-click it and pick **WEIGHT > Hide other charts**), and set its X to SEX. Its Y axis now shows avg(WEIGHT), an aggregate over the two sexes.
3. Open the Formula Lines dialog on the scatter plot and switch to the **DataFrame** tab.
4. Click **ADD NEW** and pick **Line - Horizontal**. The **Constant line** section opens on **Column** = WEIGHT and **Value** = 77.4, the column's median. Set **Value** to `100` — away from the default, so the line's place proves the edit — and **Title** to `Reference weight`. Click **OK**.
5. The scatter plot draws the line "Reference weight" at WEIGHT 100.
6. The line chart, whose Y axis is avg(WEIGHT), does not draw it: the line is built on the raw WEIGHT column, and neither sex's average is anywhere near 100.
7. Set the line chart's X to USUBJID, so that every row has its own X value and the Y axis carries WEIGHT itself. The chart draws a line through the 1000 values (a polyline, no markers).
8. The line chart now draws the line "Reference weight" at 100 as well (#2747).
9. Open the Formula Lines dialog on the line chart. The preview shows the viewer's own axes, X = USUBJID and Y = WEIGHT. Pick AGE in the preview's X column selector: the preview switches to AGE, while the line chart keeps USUBJID. Close the dialog with **CANCEL**.
10. Open the Formula Lines dialog on the scatter plot. The preview shows X = AGE and Y = WEIGHT. Pick HEIGHT in the preview's X column selector: the preview switches, and the scatter plot keeps AGE. Close the dialog with **CANCEL**.
11. Set the line chart's X back to SEX. Delete the line on the **DataFrame** tab and click **OK**. Neither viewer draws a formula line.
12. No errors appear in the console.

## Scenario 8: Line chart formula lines come back with a saved layout

1. Add a Line chart, keep only the HEIGHT chart, and set its X to AGE.
2. Open the Formula Lines dialog. Add **Line - Horizontal** (a constant line at 168.8) and **Band - Horizontal** (a band from 164.9 to 171.8), and click **OK**.
3. The line chart draws two formula items: the line "avg(HEIGHT) = 168.8" and the band "avg(HEIGHT) in(164.9, 171.8)".
4. In the top menu, pick **View > Layout > Save to Gallery**. The status bar shows "Saving layout...", and the context panel says "Layout saved".
5. Close the line chart with **General > Close** from its context menu. Only the grid is left.
6. In the top menu, pick **View > Layout > Open Gallery**, type `demog-1000` into the layout filter, and click the layout just saved.
7. The line chart comes back and draws the same two formula items, the line and the band (GROK-19943).
8. Right-click the saved layout in the gallery, pick **Delete**, and confirm with **DELETE**. The layout is gone from the gallery.
9. No errors appear in the console.
