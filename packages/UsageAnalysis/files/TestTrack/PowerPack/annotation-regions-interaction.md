# Annotation regions: hovering and clicking regions

These scenarios check how a user works with annotation regions that are already on a viewer.
Pointing at a region shows its title, description and row count, and highlights its rows as the
mouse-over group in every viewer of the table. Pointing at an overlap of two regions describes the
rows the regions share. Clicking a region selects its rows, and Ctrl+click on a fully selected region
removes them again. A region reacts to the pointer only where no marker, bar or bin lies under it.

## Setup

1. Close all views.
2. Open the demog-1000 dataset (`System:DemoFiles/demog-1000.csv`) and wait for the grid to load.

Viewers are added from **Toolbox > Viewers** and configured through the controls on the plot, with
the same configurations as in `annotation-regions.md`. Regions are added with **Tools > Formula
Lines... > ADD NEW > Region - Formula Lines**. In a formula field, typing `$` opens the column list,
and the column is picked from that list.

## Scenario 1: Hover and click overlapping regions on a two-axis viewer

1. Add a `<viewer>` and apply `<configuration>`.
2. Add a region with **Formula 1** = `<A formula 1>`, **Formula 2** = `<A formula 2>` and **Title** = `<A>`, and click **OK**.
3. Add a region with **Formula 1** = `<B formula 1>`, **Formula 2** = `<B formula 2>` and **Title** = `<B>`, and click **OK**.
4. The viewer shows two regions, each with its title drawn next to it. If markers cover the regions densely, zoom in with the mouse wheel until empty spots show inside them.
5. Point at an empty spot inside `<A>` outside `<B>`. The viewer reports one hovered region, its rows are highlighted as the mouse-over group in the grid, and the tooltip says `<A>` and `<A rows> rows`.
6. Point at an empty spot inside the overlap of `<A>` and `<B>`. The viewer reports two hovered regions, and the tooltip lists both titles and `<AB rows> rows`, the rows that lie in both regions.
7. Move the pointer out of the viewer. No region is hovered any more.
8. Click the empty spot inside `<A>` outside `<B>`. `<A rows>` rows are selected in the grid and in the viewer.
9. Hold Ctrl and click the same spot. The rows of `<A>` are removed from the selection, and no rows are selected.
10. Click the empty spot inside the overlap. `<AB rows>` rows are selected: the rows that lie in both regions.
11. Pick **Select > None** in the top menu. No rows are selected. Delete both regions through **Tools > Formula Lines...** and click **OK**.
12. No errors appear in the console.

| viewer | configuration | A | A formula 1 | A formula 2 | A rows | B | B formula 1 | B formula 2 | AB rows |
|---|---|---|---|---|---|---|---|---|---|
| Scatter plot | X = WEIGHT, Y = HEIGHT | Tall | `${HEIGHT} = 180` | `${HEIGHT} = 200` | 159 | Heavy | `${WEIGHT} = 110` | `${WEIGHT} = 170` | 21 |
| Density plot | X = WEIGHT, Y = HEIGHT | Tall | `${HEIGHT} = 180` | `${HEIGHT} = 200` | 159 | Heavy | `${WEIGHT} = 110` | `${WEIGHT} = 170` | 21 |
| Line chart | HEIGHT chart only, X = AGE | Medium height | `${avg(HEIGHT)} = 165` | `${avg(HEIGHT)} = 175` | 822 | Older | `${AGE} = 60` | `${AGE} = 90` | 102 |

## Scenario 2: Hover and click a region on a one-axis viewer

1. Add a `<viewer>` and apply `<configuration>`.
2. Right-click an empty spot of the plot, pick **Tools > Formula Lines...**, click **ADD NEW** and pick `<region item>`. Keep the formulas `<formulas>`, type `Middle age` into **Title**, and click **OK**.
3. Point at a spot inside the region above the bars or between the boxes, where nothing else is under the pointer. The viewer reports one hovered region, and the tooltip says "Middle age" and `<rows> rows`.
4. Click that spot. `<rows>` rows are selected.
5. Hold Ctrl and click the same spot. No rows are selected.
6. Delete the region through **Tools > Formula Lines...** and click **OK**.
7. No errors appear in the console.

| viewer | configuration | region item | formulas | rows |
|---|---|---|---|---|
| Histogram | value = AGE | Region - Vertical Range | `${AGE} = 36.0`, `${AGE} = 56.0` | 517 |
| Box plot | value = AGE, category = RACE | Region - Horizontal Range | `${AGE} = 36.0`, `${AGE} = 56.0` | 517 |
| Bar chart (vertical) | split = RACE, value = avg(AGE), vertical | Region - Horizontal Range | `${AGE} = 45.7`, `${AGE} = 47.1` | 89 |
