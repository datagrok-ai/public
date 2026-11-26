# More viewers

Playlist: **More viewers**

Each clip is one viewer: what it is for, the column choices on the viewer itself, then one or two
settings from the sections of its properties that matter most. Adding viewers in general, the scatter
plot, histogram, bar chart, line chart and how viewers work together are in the Viewers playlist.

## Plan

1. Box plot: value and category, the statistics and the p-value, a second category, quartiles, dot colors, the violin style. Done (01).
2. Pie chart: category, a slice hovered and clicked, sorting, the donut, the counts. Done (02).
3. Trellis plot: small multiples split by two categories, the inner viewer type, one scale, click to filter, the inner chart's own settings. Done (03).
4. Tree map and sunburst: a hierarchy of two categories, color, click to select. Done as one clip (04).
5. Correlation plot: hover and double-click a cell, Spearman, the columns of one axis. Done (05).
6. PC plot (parallel coordinates): color, filtering on an axis, the axes' columns and scale. Done (06).
7. Statistics viewer (07); then, as extra clips, the density plot (11), the heat map (12) and the matrix plot (13).
8. Charts package viewers: the Sankey diagram (08) and the radar (09). Word cloud and timelines left out (less used; no demo
   table with event dates).
9. Form and Tile viewers: rows as cards. Done as one clip (10).
10. Extra: the network diagram (14) and the calendar (15).

## 01 - Compare groups with a box plot (1:22)

Title: Compare groups with a box plot in Datagrok

Description: A box plot shows how a number is spread in each group: the box holds the middle half of the
values, the line is the median, and the dots outside the whiskers are outliers. This video adds one from the
toolbox, picks the value and the category right on the plot, and reads the statistics and the p-value under it.
The settings icon then splits the groups by a second category, adds quartiles to the statistics, colors the dots
by age and switches to a violin plot. The data is the demog demo table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/01-box-plot.feature

## 02 - Show shares with a pie chart (1:00)

Title: Pie and donut charts in Datagrok

Description: A pie chart splits the rows by one category and shows each category's share. This video adds one
from the toolbox, picks the category on the chart, and hovers and clicks a slice to see and select its rows.
In the settings it sorts the slices largest first, turns the pie into a donut and shows the row count on each
slice. The data is the demog demo table (Files > Demo), split by SEVERITY.

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/02-pie-chart.feature

## 03 - Small multiples with a trellis plot (1:17)

Title: Small multiples with a trellis plot in Datagrok

Description: A trellis plot draws the same chart once for every category, so groups can be compared side by
side. This video adds one from the toolbox, splits it by SEX across and RACE down with the selectors on the plot,
and puts a scatter plot colored by age in every cell. In the settings it gives all cells one scale, makes a click
on a cell filter the table, and turns on a regression line from the inner chart's own tab. The data is the demog
demo table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/03-trellis-plot.feature

## 04 - Show a hierarchy with a tree map or a sunburst (1:36)

Title: Tree map and sunburst: categories inside categories in Datagrok

Description: A tree map and a sunburst both show categories nested inside other categories. This video adds a
tree map from the toolbox, drags it by its title bar to the right of the grid, picks its levels (SEX, then DIS_POP) and a color column right on the map, and clicks a
rectangle to select its rows. Then it adds a sunburst from the Add viewer gallery and sets the same levels in its
Hierarchy setting; a click on a segment selects those rows in both viewers. The data is the demog demo table
(Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/04-tree-map.feature

## 05 - Find related columns with a correlation plot (1:06)

Title: Find related columns with a correlation plot in Datagrok

Description: A correlation plot shows how strongly every pair of numeric columns moves together: red cells rise
together, blue ones move in opposite directions. This video adds one from the toolbox and hovers a cell to see
that pair as a small scatter plot. In the settings it switches from Pearson to Spearman and keeps only mileage and
price on one axis; a double-click on a cell then opens the pair as a full scatter plot. The data is the cars demo
table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/05-correlation-plot.feature

## 06 - Compare many columns with a PC plot (1:02)

Title: Parallel coordinates (PC plot) in Datagrok

Description: A PC plot, or parallel coordinates plot, draws every row as a line across one vertical axis per
column, so you can follow a group across many measurements at once. This video adds one from the toolbox, colors
the lines by species and drags the top handle of one axis to filter the table. In the settings it keeps the
filtered-out lines in grey, removes the row-number column from the axes and puts all four measurements on one
scale; Reset View in the plot's menu clears the filter. The data is the iris demo table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/06-pc-plot.feature

## 07 - Summarize every column with the Statistics viewer (0:51)

Title: Summarize every column with the Statistics viewer in Datagrok

Description: The Statistics viewer gives one row per column of your table: how many values, how many are
missing, how many unique, and min, max, average, median and standard deviation. This video adds it from the
toolbox, drags it by its title bar to the bottom so every statistic fits, uses the right-click menu of a row to add a sum and a small histogram split by sex, and then sets a
filter formula in its settings so the numbers cover only people over 60. The data is the demog demo table
(Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/07-statistics.feature

## 08 - Show flows with a Sankey diagram (0:43)

Title: Show flows with a Sankey diagram in Datagrok

Description: A Sankey diagram shows how an amount flows from one category to the next; the wider the band, the
bigger the flow. This video adds one from the Add viewer gallery to a table of UK energy flows, where it picks the
source, target and amount columns by itself. A click on a node selects the rows of its flows, and the Color
section of the settings colors each band by where it comes from. The data is energy_uk.csv from Files > Demo;
the Sankey comes with the Charts package.

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/08-sankey.feature

## 09 - Compare rows across many columns with a radar (1:10)

Title: Compare rows across many columns with a radar chart in Datagrok

Description: A radar chart draws a row as a shape over several axes, one per column, on top of a band that shows
where most rows fall. This video adds one from the Add viewer gallery, picks five columns for its axes in the
Values setting and turns on Show Values. Clicking a car in the grid then draws that car's shape, so a sports car
and a family car look clearly different. The data is the cars demo table (Files > Demo); the radar comes with the
Charts package.

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/09-radar.feature

## 10 - Show rows as cards with the Tile Viewer and the Form (1:05)

Title: Rows as cards: the Tile Viewer and the Form in Datagrok

Description: The Tile Viewer shows every row of a table as a card, and the Form shows one row at a time as a form.
This video adds the Tile Viewer from the Add viewer gallery, drags it by its title bar to the bottom, splits its
cards into lanes by severity in the settings, and drags a card to another lane, which changes that row's value in
the table. Then it hides the context panel and adds a Form:
clicking a card shows that row in the form, and the arrows at the top step through the rows. The data is the
demog demo table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/10-tile-viewer.feature

## 11 - See where points crowd with a density plot (0:55)

Title: See where points crowd with a density plot in Datagrok

Description: When a scatter plot has too many points to read, a density plot counts them in small bins and colors
each bin by its count. This video adds one from the toolbox, plots earthquake longitude against latitude so the
plate boundaries show up, and hovers a bin to see its range and how many quakes it holds. In the settings it
raises the number of bins and switches to a logarithmic color scale so the quieter areas show too. The data is the
earthquakes demo table (Files > Demo > geo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/11-density-plot.feature

## 12 - See a whole table at once with a heat map (1:00)

Title: See a whole table at once with a heat map in Datagrok

Description: A heat map draws every value of a table as a colored cell, high values red and low ones blue, so rows
that look alike stand out. This video adds one from the toolbox and double-clicks the price header to sort by it:
the big, fast cars line up at the top. In the settings it hides the yes/no columns through the Columns setting.
The data is the cars demo table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/12-heat-map.feature

## 13 - See every pair of columns with a matrix plot (1:12)

Title: See every pair of columns with a matrix plot in Datagrok

Description: A matrix plot draws a small plot for every pair of numeric columns, with a histogram of each column
on the diagonal, so all the relationships in a table can be scanned at once. This video adds one from the toolbox,
switches the cells to scatter plots, removes the row-number column from both axes, and then uses the Scatter
plot tab of the settings to color every cell by species. The data is the iris demo table (Files > Demo).

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/13-matrix-plot.feature

## 14 - Show links between things with a network diagram (1:00)

Title: Show links between things with a network diagram in Datagrok

Description: A network diagram draws each value as a node and each row as a link between two nodes. This video
adds one from the Add viewer gallery to a table of Game of Thrones relationships and picks the two columns each
link joins right on the viewer. In the settings it makes the nodes bigger for characters with more links and the
links wider for stronger ties, and at the end a click on a character selects the rows of their links. The data is
got-s1-edges.csv from Files > Demo.

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/14-network-diagram.feature

## 15 - See dates on a calendar (0:37)

Title: See dates on a calendar in Datagrok

Description: The calendar lays a date column out by day, month and weekday; the bigger the circle, the more rows
fall on that day. This video adds it from the toolbox, clicks a month to select its rows and a weekday header to
select every Sunday. In the settings it switches On Click to Filter, so a click on a month keeps only its rows in
the table. The data is the demog demo table (Files > Demo), by its STARTED date.

Feature: public/packages/Tutorials/bdd/features/howto/03-more-viewers/15-calendar.feature

