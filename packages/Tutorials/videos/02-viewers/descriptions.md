# Viewers

Playlist: **Viewers**

## Plan

Every clip is one idea, filmed from a feature in `public/packages/Tutorials/bdd/features/howto/02-viewers/`.
Demo data opened off camera (demog unless noted); nothing is saved on the server.

1. **Add viewers to a table** - the toolbox's Viewers icons and the Add viewer gallery (tags, search); then
   arrange them: dock one beside or under the grid, full screen, close.
2. **Scatter plot** - X, Y, Color and Size from the selectors on the plot; zoom and pan; a regression line;
   then the gear icon and one or two settings in each section (markers, labels, axes).
3. **Histogram** - the value column, the bins slider, split and color by a category, its properties.
4. **Bar chart** - category, value and aggregation, stacking, sorting, orientation.
5. **Line chart** - a date or number on X, several Y columns, split by a category, aggregation.
6. **Viewers work together** - selection, the current row and hover shared by the grid and every viewer;
   filtering applies to all.
7. **Viewer properties in general** - the gear icon, the sections of the context panel, title and
   description, legend, style (nothing saved as a default).
8. **Pivot table** - group by, pivot and aggregate; ADD the result as a table (added: the next most used
   viewer that no other section covers).
9. **Tooltips** - what a tooltip shows and Tooltip > Edit... (added).
10. **Reuse a viewer** - Pick Up / Apply Style Settings, and General > Clone (added).
11. **Show only some rows in a viewer** - Data > Filter formula and Row Source (added).
12. **Use the legend** - click, Ctrl+click and the x on categories; recolor one (added).

## 01 - Add viewers to a table (1:14)

Title: Add viewers to a table in Datagrok

Description: There are two ways to add a viewer to an open table. The toolbox on the left has icons for
the common viewers, one click each; the Add viewer icon in the toolbar opens a gallery of every viewer,
which you can narrow with a tag or a search. The video then drags a viewer by its title bar to dock it
under the others, shows one full screen (Alt+F, Option+F on a Mac, or Esc to go back) and closes one.
It uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/01-add-viewers.feature`

## 02 - Configure a scatter plot (1:42)

Title: Scatter plot in Datagrok: columns, zoom, regression line, settings

Description: Pick X and Y in the selectors on the axes, and Color and Size in the ones that appear at the
top right of the plot. The mouse wheel zooms in (the table then keeps only the rows you see), dragging
moves around, and a double-click on an empty spot shows everything again. Right-click > Tools adds a
regression line. The gear icon on the title bar opens all the settings in the context panel; the video
changes the marker shape, labels each point with the car model, adds a histogram along the X axis and
shows the correlation next to the fit. It uses the cars demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/02-scatter-plot.feature`

## 03 - Configure a histogram (1:07)

Title: Histogram in Datagrok: bins, range filter, split by category

Description: Pick the column in the selector under the histogram, and set the number of bins with the
small slider that appears at the top when you point at the plot. Hovering a bar shows how many rows it
holds; the handles under the bars filter the whole table, and a double-click on that slider shows all rows
again. The split selector at the top right draws one line per category, and in the settings (the gear
icon) Split Stack turns them into stacked bars and Show Values writes the counts. It uses the demog demo
dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/03-histogram.feature`

## 04 - Configure a bar chart (1:21)

Title: Bar chart in Datagrok: category, aggregation, stacking, sorting

Description: The selector beside the bars picks the category; the one under them picks the column and
how it is aggregated (count, average, sum and so on). Stacking splits each bar by a second category and
needs a count or a sum. Right-click the chart to sort the bars (Order) or turn them upright
(Orientation), and in the settings (the gear icon) Relative Values shows each bar as 100%. It uses the
demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/04-bar-chart.feature`

## 05 - Configure a line chart (1:12)

Title: Line chart in Datagrok: dates on X, split, several Y columns

Description: X can be a number or a date; for a date, the list next to the X selector groups the values
by year, month, week or day, and the chart then combines each group's rows. The plus at the top left
splits the line by a category. In the settings (the gear icon) you add more Y columns, each in its own
chart, or turn on Multi Axis to draw them together; the Y Aggr Type list then picks how the rows are
combined (average, median, sum...), and Interpolation smooths the lines. It uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/05-line-chart.feature`

## 06 - Viewers work together (0:38)

Title: Datagrok viewers are linked: hover, current row, selection, filter

Description: The grid and every viewer of a table work on the same rows. Pointing at a histogram bin
highlights its rows on the scatter plot, a row clicked in the grid becomes the current row in every
viewer, and rows selected in one viewer (Shift+drag, or a click on a bar) show as selected in all the
others. A filter from the Filters panel applies to the grid and every viewer at once. It uses the demog
demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/06-viewers-work-together.feature`

## 07 - Change a viewer's settings (1:08)

Title: Viewer settings in Datagrok: title, legend, search, style

Description: Every viewer keeps its settings in the context panel on the right; the gear icon on the
viewer's title bar opens them, grouped in sections. The video gives a scatter plot a title and a
description, moves its legend to the bottom, then uses the search box at the top of the panel to find
the grid line switches and the background color. The changes apply to this viewer only. It uses the demog
demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/07-viewer-settings.feature`

## 08 - Summarize a table with a pivot table (0:50)

Title: Pivot table in Datagrok: group, pivot and aggregate

Description: The pivot table viewer has three rows of chips: Group by makes one row per value, Pivot makes
one column per value, and Aggregate says what to compute in each cell. The plus on each row adds a column
(or drag a column header from the grid onto a row), the x on a chip removes it, and a right-click on an
aggregate chip changes its function (average, max, count and so on). ADD puts the result into a new
table. It uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/08-pivot-table.feature`

## 09 - Choose what a tooltip shows (0:40)

Title: Choose the columns a Datagrok tooltip shows

Description: Pointing at a marker on a scatter plot shows that row's values in a tooltip. To change which
columns it lists, right-click the plot and pick Tooltip > Edit..., untick the columns you don't need and
click OK; the tooltip lists up to seven columns, top to bottom. The dialog edits the table's shared
tooltip ("Use tooltip: Table (shared)"), which other viewers of the table use as well. It uses the demog
demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/09-tooltips.feature`

## 10 - Copy a viewer or its look (0:34)

Title: Copy a Datagrok viewer, or copy its look to another viewer

Description: Once one viewer looks the way you want, you can reuse that look. Right-click it and pick
Pick Up / Apply > Pick Up, then right-click another viewer of the same kind and pick Apply Style Settings:
it takes the look (here hollow markers and the axis histograms) and keeps its own columns, while Apply
copies the columns too. General > Clone adds a copy of a viewer next to it. It uses the demog-1000 demo
dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/10-reuse-a-viewer.feature`

## 11 - Show only some rows in a viewer (0:39)

Title: Show only some rows in one Datagrok viewer

Description: A viewer can show fewer rows than its table without filtering the table. In the viewer's
settings (the gear icon), Data > Filter takes a formula such as ${AGE} > 44, and only this viewer
changes: the grid and any other viewer keep all rows. Row Source picks which rows the viewer starts from:
all, the filtered ones, or, as here, only the rows selected in the grid. It uses the demog-1000 demo
dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/11-viewer-rows.feature`

## 12 - Use the legend (0:39)

Title: Use a Datagrok viewer's legend: show categories, change colors

Description: The legend lists the categories behind a viewer's colors. Click one to show only that
category on this viewer (the table itself is not filtered), Ctrl+click (Cmd+click on a Mac) to add
another, and click the x next to one to take it back. Pointing at an item shows a palette icon that
changes the category's color. It uses the demog-1000 demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/02-viewers/12-legend.feature`

