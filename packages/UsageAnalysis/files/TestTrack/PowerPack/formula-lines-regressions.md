# Formula lines: regression checks

These scenarios guard fixed defects around formula lines on viewers. A column used on both axes can
be renamed while a formula line is configured. Renaming any column used by a formula line rewrites
the line. A formula line survives a change of the axis columns. A band and a line survive a switch to
a logarithmic axis. Hovering markers next to a formula line raises no errors. Each defect is named
next to the check that guards it.

## Setup

1. Close all views.
2. Open the demog-1000 dataset (`System:DemoFiles/demog-1000.csv`) and wait for the grid to load.

Viewers are added from **Toolbox > Viewers**. Their columns are set with the selectors on the plot.
The Formula Lines dialog is opened by right-clicking an empty spot of the plot and picking **Tools >
Formula Lines...**. A column is renamed in the grid: right-click its header, pick **Column
Properties...**, type the new name into **New name**, and click **OK**.

## Scenario 1: Rename a column that is on both axes while a formula line is configured

1. Add a Scatter plot and set both X and Y to AGE.
2. Open the Formula Lines dialog, click **ADD NEW**, pick **Line**, and click **OK**. The viewer holds one formula line, `${AGE} = ${AGE}`.
3. Rename the AGE column to `AGE_RENAMED`.
4. The scatter plot's X column and Y column are both AGE_RENAMED (GROK-19334).
5. The formula line now reads `${AGE_RENAMED} = ${AGE_RENAMED}` and is still drawn.
6. No errors appear in the console during the rename.
7. Rename the column back to AGE. X, Y and the formula line refer to AGE again.
8. Delete the line through the dialog, click **OK**, and set Y back to HEIGHT.
9. No errors appear in the console.

## Scenario 2: Renaming a column rewrites the formula lines that use it

1. Add a `<viewer>` and apply `<configuration>`.
2. `<add the line>`. The viewer holds one formula line, `<formula>`.
3. Rename the `<column>` column to `<new name>`.
4. The viewer's formula line now reads `<renamed formula>`, and it is still drawn.
5. Rename the column back to `<column>`. The formula line reads `<formula>` again.
6. Delete the line through the Formula Lines dialog and click **OK**.
7. No errors appear in the console.

| viewer | configuration | add the line | formula | column | new name | renamed formula |
|---|---|---|---|---|---|---|
| Scatter plot | X = AGE, Y = HEIGHT | Open the dialog, **ADD NEW > Line**, **OK** | `${HEIGHT} = ${AGE}` | AGE | AGE_R | `${HEIGHT} = ${AGE_R}` |
| Line chart | HEIGHT chart only (**HEIGHT > Hide other charts**), X = AGE | Right-click the Y axis, go straight down to **Annotations**, then right to **Add Line**, **OK** | `${avg(HEIGHT)} = 168.8` | HEIGHT | HEIGHT_R | `${avg(HEIGHT_R)} = 168.8` |

## Scenario 3: A formula line survives a change of the axis columns

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT.
2. Open the Formula Lines dialog, click **ADD NEW**, pick **Line**, and click **OK**. The viewer holds one formula line, `${HEIGHT} = ${WEIGHT}`.
3. Set X to AGE and Y to WEIGHT with the selectors on the plot.
4. The viewer's stored formula lines are exactly the same as before the change: one line, `${HEIGHT} = ${WEIGHT}` (GROK-16214).
5. Set X back to WEIGHT and Y back to HEIGHT. The line is drawn again.
6. Delete the line through the dialog and click **OK**.
7. No errors appear in the console.

## Scenario 4: Formula items survive a switch to a logarithmic Y axis

1. Add a `<viewer>` and apply `<configuration>`.
2. `<add the item>`. The viewer draws one formula item.
3. Right-click the Y axis and pick **Y Axis Type > Logarithmic**. The Y axis is logarithmic, and the viewer still draws the formula item, placed on the logarithmic scale.
4. Right-click the Y axis and pick **Y Axis Type > Linear**. The Y axis is linear, and the item is still drawn.
5. No errors appear in the console during the switch and back (`<guards>`).
6. Delete the item through the Formula Lines dialog and click **OK**.
7. No errors appear in the console.

| viewer | configuration | add the item | guards |
|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | Open the dialog, **ADD NEW > Band - Horizontal**, **OK** | GROK-20458 |
| Line chart | HEIGHT chart only (**HEIGHT > Hide other charts**), X = AGE | Right-click the Y axis, go straight down to **Annotations**, then right to **Add Line**, **OK** | formula line on a logarithmic line chart axis |

## Scenario 5: Hovering markers next to a formula line raises no errors

1. Add a Scatter plot with X = WEIGHT and Y = HEIGHT.
2. Open the Formula Lines dialog, click **ADD NEW**, pick **Line**, and click **OK**. The viewer holds one formula line.
3. Move the pointer over eight different markers in turn, and pause on each until its tooltip appears.
4. Each pause shows the row tooltip of the marker under the pointer (HEIGHT, WEIGHT, USUBJID and so on).
5. No errors appear in the console during the sweep (github-2530).
6. Delete the line through the dialog and click **OK**. The viewer holds no formula lines.
7. No errors appear in the console.
