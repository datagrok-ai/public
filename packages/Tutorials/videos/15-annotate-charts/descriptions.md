# Annotate and explain charts

Playlist: **Annotate and explain charts**. Features: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/`
(feature name `How-to · Annotations`).

Setting charts up is in the Viewers (02) and More viewers (03) playlists; these clips add meaning to a chart that
is already there. The demo table is opened and the chart added off camera; nothing is saved on the server.

## Plan

| NN | Clip | One idea | Data |
|---|---|---|---|
| 01 | Draw reference lines with formula lines | A line at a value from an axis's menu, a curve from a formula, its title and color | demog-1000 |
| 02 | Mark a group of points with an annotation region | Draw a region, name it, hover and click it | iris |
| 03 | Label points with structures | Structures and IDs as labels, for the selected points only | spgi-100 |
| 04 | Color points by conditions | Conditional color coding of a column, its ranges, the chart follows | demog-1000 |
| 05 | Compare each group with a control | The box plot's p-value, then each group against a control group | demog |
| 06 | Shade a range of values on a chart | A region on a histogram's axis: the quartiles, then a typed range | demog-1000 |
| 07 | Add a note that explains a chart | A Markup viewer whose ${COLUMN} references follow the clicked point | demog-1000 |

Left out: the regression line and its statistics (in 02-viewers, Configure a scatter plot) and the tooltip editor
(02-viewers, Choose what a tooltip shows).

## 01 - Draw reference lines with formula lines (1:02)

Title: Reference lines and curves on a Datagrok chart (formula lines)

Description: Formula lines draw a value or a formula over a chart. Right-click an axis and pick Annotations > Add
Line for a line at one value; the Formula Lines dialog opens on it, and Add new > Line there takes any formula of X
and Y, here the weight at a body mass index of 25. The dialog also sets the line's title, color and style, and
Tools > Formula Lines... in the chart's menu lists and edits them later. It uses the demog-1000 demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/01-formula-lines.feature`

## 02 - Mark a group of points with an annotation region (0:44)

Title: Mark and name a group of points on a Datagrok scatter plot

Description: An annotation region marks part of a chart and gives it a name. Right-click the scatter plot, pick
Tools > Draw Annotation Region and drag a box around the points; the dialog that opens takes the region's title
and a note. Pointing at the region shows both, with the number of rows inside, and a click on it selects those
rows. It uses the iris demo dataset, where the box goes round the petals on which versicolor and virginica mix.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/02-annotation-regions.feature`

## 03 - Label points with structures (0:55)

Title: Show structures as labels on a Datagrok scatter plot

Description: A scatter plot can write any column's values next to its points, and a molecule column is drawn as
the structure itself. In the plot's settings (the gear icon), Labels > Label Columns picks the structure and ID
columns, and Label Content Size makes the drawings bigger. Show Labels For, also in the plot's right-click menu
under Labels, can keep labels for the selected points only: Shift+drag round a few points and their structures
appear. It uses the spgi-100 demo dataset of 100 molecules.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/03-label-points-with-structures.feature`

## 04 - Color points by conditions (1:05)

Title: Color a Datagrok chart by ranges of values (conditional colors)

Description: Conditional color coding gives each range of a number column its own color, and every chart colored
by that column uses the same colors. Right-click the column's header in the grid and pick Color Coding >
Conditional; Color Coding > Edit... then sets the ranges: Bin... says how many, a click on a color changes it,
and each range is typed as a rule such as <40, 40-65 or >65. The video colors a scatter plot by age with the
demog-1000 demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/04-color-by-conditions.feature`

## 05 - Compare each group with a control (0:43)

Title: Compare each group with a control group on a Datagrok box plot

Description: The p-value at the top of a box plot says whether the groups differ at all; the icon beside it
compares each group with one control group instead. Pick the control in the selector that appears, and each other
group is
marked *** if it differs from the control or ns if it does not; pointing at a mark shows its adjusted p-value.
Right-click the comparison and pick Add Control Comparisons Table to get the numbers as a table. The video
compares height by race in the demog demo dataset, with Caucasian as the control: Asian and Other differ, Black
does not.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/05-compare-groups-with-a-control.feature`

## 06 - Shade a range of values on a chart (0:55)

Title: Shade and name a range of values on a Datagrok chart

Description: A shaded region marks a range of one column on a chart. Right-click the X axis of a histogram and
pick Annotations > Add Region: the region starts from the first to the third quartile, the middle half of the
values, and the dialog that opens gives it a title. A second region gets its bounds typed as formulas, here ages 65
to 89, and pointing at a region shows its name and how many rows fall in it. Bar charts, line charts and box plots
offer the same item on their value axis. It uses the demog-1000 demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/06-shade-a-range-of-values.feature`

## 07 - Add a note that explains a chart (0:44)

Title: Add a text note next to your Datagrok charts (Markup viewer)

Description: The Markup viewer puts a note next to your charts. Add it from the Add viewer gallery, then right-click
it and pick Edit content... to type the text in Markdown. A column name written as ${COLUMN} shows that column's
value for the current row, so clicking a point on the scatter plot changes the note to describe that patient. It
uses the demog-1000 demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/15-annotate-charts/07-add-a-note.feature`
