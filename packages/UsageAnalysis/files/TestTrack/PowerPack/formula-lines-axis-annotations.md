# Formula lines: the Annotations group on viewer axes

These scenarios check the Annotations group in the context menu of a numeric viewer axis. It holds
**Add Line**, **Add Band** and **Add Region**. Each item adds an item bound to the axis column,
places it at the column's quartiles, and opens the Formula Lines dialog (PowerPack) on the new item.
Categorical axes and derived axes, such as the count axis of a histogram, have no Annotations group.

## Setup

1. Close all views.
2. Open the demog-1000 dataset (`System:DemoFiles/demog-1000.csv`) and wait for the grid to load.

Viewers are added from **Toolbox > Viewers** and configured through the controls on the plot, with
the same configurations as in `annotation-regions.md`. On the scatter plot, density plot and
line chart, the X axis is right-clicked near its left end, away from the column selector under it.

## Scenario 1: Add Line places a constant line at the median of the axis column

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click `<axis>`. The context menu has an **Annotations** group.
3. Move the pointer straight down the menu to **Annotations**, then right along its row into the submenu, and click **Add Line**.
4. The **Formula Lines** dialog opens on the **Viewer** tab with the new line selected. Its **Constant line** section shows **Column** = `<column>` and **Value** = `<median>`.
5. Click **OK**. The viewer holds one formula line, `<formula>`, drawn as a `<orientation>` line across the plot.
6. Delete the line through **Tools > Formula Lines...** (trash button, **OK**). The viewer holds no formula lines.
7. No errors appear in the console.

| viewer | configuration | axis | column | median | formula | orientation |
|---|---|---|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | the X axis | WEIGHT | 77.4 | `${WEIGHT} = 77.4` | vertical |
| Scatter plot | X = WEIGHT, Y = HEIGHT | the Y axis | HEIGHT | 168.5 | `${HEIGHT} = 168.5` | horizontal |
| Density plot | X = WEIGHT, Y = HEIGHT | the X axis | WEIGHT | 77.4 | `${WEIGHT} = 77.4` | vertical |
| Density plot | X = WEIGHT, Y = HEIGHT | the Y axis | HEIGHT | 168.5 | `${HEIGHT} = 168.5` | horizontal |
| Line chart | HEIGHT chart only, X = AGE | the X axis | AGE | 45 | `${AGE} = 45.0` | vertical |
| Line chart | HEIGHT chart only, X = AGE | the Y axis | avg(HEIGHT) | 168.8 | `${avg(HEIGHT)} = 168.8` | horizontal |
| Histogram | value = AGE | the X axis (value axis) | AGE | 45 | `${AGE} = 45.0` | vertical |
| Box plot | value = AGE, category = RACE | the Y axis (value axis) | AGE | 45 | `${AGE} = 45.0` | horizontal |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | the value axis on the left | AGE | 46.2 | `${AGE} = 46.2` | horizontal |
| Bar chart (horizontal) | split = RACE, value = avg(AGE), horizontal | the value axis at the bottom | AGE | 46.2 | `${AGE} = 46.2` | vertical |

## Scenario 2: Add Band places a band between the quartiles of the axis column

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click `<axis>`, move the pointer straight down to **Annotations**, then right into the submenu, and click **Add Band**.
3. The **Formula Lines** dialog opens with the new band selected, and its **Band** section is shown.
4. Click **OK**. The viewer holds one formula band, `<formula>`, drawn as a `<orientation>` strip.
5. Delete the band through **Tools > Formula Lines...** (trash button, **OK**). The viewer holds no formula lines.
6. No errors appear in the console.

| viewer | configuration | axis | formula | orientation |
|---|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | the X axis | `${WEIGHT} in (64.2, 91.0)` | vertical |
| Scatter plot | X = WEIGHT, Y = HEIGHT | the Y axis | `${HEIGHT} in (160.9, 177.6)` | horizontal |
| Density plot | X = WEIGHT, Y = HEIGHT | the X axis | `${WEIGHT} in (64.2, 91.0)` | vertical |
| Density plot | X = WEIGHT, Y = HEIGHT | the Y axis | `${HEIGHT} in (160.9, 177.6)` | horizontal |
| Line chart | HEIGHT chart only, X = AGE | the X axis | `${AGE} in (36.0, 56.0)` | vertical |
| Line chart | HEIGHT chart only, X = AGE | the Y axis | `${avg(HEIGHT)} in (164.9, 171.8)` | horizontal |
| Histogram | value = AGE | the X axis (value axis) | `${AGE} in (36.0, 56.0)` | vertical |
| Box plot | value = AGE, category = RACE | the Y axis (value axis) | `${AGE} in (36.0, 56.0)` | horizontal |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | the value axis on the left | `${AGE} in (45.7, 47.1)` | horizontal |
| Bar chart (horizontal) | split = RACE, value = avg(AGE), horizontal | the value axis at the bottom | `${AGE} in (45.7, 47.1)` | vertical |

## Scenario 3: Add Region adds a region between the quartiles of the axis column

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click `<axis>`, move the pointer straight down to **Annotations**, then right into the submenu, and click **Add Region**.
3. The **Formula Lines** dialog opens with the new region selected. Its **Formula** section shows **Formula 1** = `<formula 1>` and **Formula 2** = `<formula 2>`. The **Title** field is empty: a region added from the axis menu gets no title of its own.
4. Click **OK**. The viewer holds one viewer region, which spans the whole plot across the other axis.
5. Delete the region through **Tools > Formula Lines...** (trash button, **OK**). The viewer holds no regions.
6. No errors appear in the console.

| viewer | configuration | axis | formula 1 | formula 2 |
|---|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | the X axis | `${WEIGHT} = 64.2` | `${WEIGHT} = 91.0` |
| Scatter plot | X = WEIGHT, Y = HEIGHT | the Y axis | `${HEIGHT} = 160.9` | `${HEIGHT} = 177.6` |
| Density plot | X = WEIGHT, Y = HEIGHT | the X axis | `${WEIGHT} = 64.2` | `${WEIGHT} = 91.0` |
| Density plot | X = WEIGHT, Y = HEIGHT | the Y axis | `${HEIGHT} = 160.9` | `${HEIGHT} = 177.6` |
| Line chart | HEIGHT chart only, X = AGE | the X axis | `${AGE} = 36.0` | `${AGE} = 56.0` |
| Line chart | HEIGHT chart only, X = AGE | the Y axis | `${avg(HEIGHT)} = 164.9` | `${avg(HEIGHT)} = 171.8` |
| Histogram | value = AGE | the X axis (value axis) | `${AGE} = 36.0` | `${AGE} = 56.0` |
| Box plot | value = AGE, category = RACE | the Y axis (value axis) | `${AGE} = 36.0` | `${AGE} = 56.0` |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | the value axis on the left | `${AGE} = 45.7` | `${AGE} = 47.1` |
| Bar chart (horizontal) | split = RACE, value = avg(AGE), horizontal | the value axis at the bottom | `${AGE} = 45.7` | `${AGE} = 47.1` |

## Scenario 4: Categorical and derived axes have no Annotations group

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click `<axis>` in its middle, near both of its ends, and at its upper and lower edge.
3. None of these context menus has an **Annotations** group or an **Add Line** item at its top level.
4. Close the menu with Escape.
5. No errors appear in the console.

| viewer | configuration | axis |
|---|---|---|
| Box plot | value = AGE, category = RACE | the X axis with the category names |
| Histogram | value = AGE | the Y axis with the counts |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | the axis with the category names |
| Bar chart (horizontal) | split = RACE, value = avg(AGE), horizontal | the axis with the category names |
