# Dose-response curves

Playlist: **Dose-response curves**. Features: `public/packages/Tutorials/bdd/features/howto/16-curves/` (feature name
`How-to · Curves`). For a scientist who analyses assay data: curves drawn in grid cells, fitted in the browser (the
Curves package). Everything shown runs in the browser; Calculate MSR (Python) is left out.

## Plan

| NN | Clip | One idea | Data |
|---|---|---|---|
| 01 | Dose-response curves in a table | A file of fitted curves from Browse: a curve in every cell, a point's values, the curve in the context panel (Chart, Fit statistics), a larger window | Curves sample Assay_Curves.csv (43 curves) |
| 02 | Change how curves are fitted and drawn | The Options pane: the fit function for the whole column, then display options (confidence band, IC80 drop lines, statistics on the plot) | Assay_Curves |
| 03 | IC50 and other statistics as columns | The Fit pane's + adds IC50 and R² for every curve as new columns; scientific format; sort by R² for the poorest fits | Assay_Curves |
| 04 | Leave out an outlier | Click a point to leave it out of the fit: the curve, its R² and IC50 follow; click it again to bring it back | Assay_Curves |
| 05 | Compare curves in one chart | The Multi curve viewer overlays the curves of the selected rows (one compound's runs) | Assay_Curves |
| 06 | Open a GraphPad Prism file | A .pzfx file from Browse opens with its XY tables as curves | Curves sample dose-response-4pl.pzfx |
| 07 | Curves per compound in a trellis plot | A trellis plot with Curves in every cell, one cell per compound; a click filters the table to it | Assay_Curves |

Not planned: Data > Curves > Data to Curves (it needs a table of wells with compound, batch, assay, target, run,
concentration and readout columns, and no demo file has that shape).

## 01 - Dose-response curves in a table (0:48)

Title: Dose-response curves in a table in Datagrok

Description: Datagrok draws fitted dose-response curves right in the grid. This video opens the Curves package's
sample file (Files > App Data > Curves > demo > Assay_Curves.csv in Browse): every cell of the Fitted Curve column
shows the points and a four-parameter fit on a log concentration axis. Point at a measurement to read its concentration
and response; click a curve and the context panel shows it larger (Chart) with its IC50, Hill slope, R² and other
statistics (Fit); double-click it to open it in a window of its own.

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/01-curves-in-a-table.feature

## 02 - Change how curves are fitted and drawn (0:48)

Title: Change how dose-response curves are fitted and drawn in Datagrok

Description: Click a curve in the grid and open Options in the context panel. With Level set to Column, every change
applies to all curves in the column: Fit Function picks the model (here a straight line, then back to the
four-parameter dose-response fit), Confidence Interval draws a band around each fit, Droplines adds IC80 next to IC50,
and Statistics prints values such as IC50 and R² on every plot. Uses the Curves sample Assay_Curves.csv (Files > App
Data > Curves > demo).

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/02-fit-and-display-options.feature

## 03 - IC50 and other statistics as columns (0:51)

Title: Get IC50 and R² for every dose-response curve as columns

Description: Click a curve and open Fit in the context panel: it lists the curve's statistics (IC50, Hill slope, R²,
AUC, pIC50 and more). The + next to a value adds it as a new column with a value for every curve. Here IC50 and R²
are added, the IC50 column is switched to scientific notation (Format in the column header's menu), and sorting by
R² brings the poorest fits to the top. Uses the Curves sample Assay_Curves.csv (Files > App Data > Curves > demo).

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/03-statistics-as-columns.feature

## 04 - Leave out an outlier (0:38)

Title: Leave an outlier out of a dose-response fit in Datagrok

Description: A point that does not belong can be left out of a curve's fit with a click. In this video one point of a
curve sits well below its neighbours; clicking it marks it as an outlier (a red cross), the curve is fitted again
without it, and the R² and IC50 columns next to the curves update (R² from 0.93 to 0.98). Clicking the cross puts
the point back. Uses the Curves sample Assay_Curves.csv (Files > App Data > Curves > demo); the R² and IC50 columns
were added from the Fit pane, as shown in the previous video.

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/04-leave-out-an-outlier.feature

## 05 - Compare curves in one chart (0:38)

Title: Compare dose-response curves in one chart in Datagrok

Description: To compare curves, select their rows and add the Multi curve viewer (Add viewer in the toolbar, search
for "curve"). It overlays the curves of the selected rows, with a legend naming each one; here the four runs of one
compound, where one run sits to the left of the others. Pointing at another row adds its curve to the chart. Uses the
Curves sample Assay_Curves.csv (Files > App Data > Curves > demo).

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/05-compare-curves.feature

## 06 - Open a GraphPad Prism file (0:37)

Title: Open a GraphPad Prism (.pzfx) file in Datagrok

Description: Datagrok opens GraphPad Prism project files (.pzfx) directly. This video double-clicks a sample Prism file
in Browse (Files > App Data > Curves > dose-response-4pl.pzfx): each XY table becomes a row, and its data sets are
fitted as dose-response curves in one cell, with the values excluded in Prism shown as outliers. Double-click the cell
for a larger chart with a legend, and the context panel's Fit section lists R², IC50, slope and more for every
compound.

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/06-open-a-prism-file.feature

## 07 - Curves per compound in a trellis plot (1:00)

Title: Dose-response curves per compound in a trellis plot

Description: A trellis plot draws one small chart per group. This video adds one from the toolbox, picks Curves as
the chart in every cell and splits the cells by Compound ID (Reset Y columns, in the right-click menu of the split
down, leaves one cell per compound), so each compound's runs are overlaid in a chart of their own. In the settings,
Misc > On Click set to Filter makes a click on a compound's chart filter the table to its runs. Uses the Curves sample
Assay_Curves.csv (Files > App Data > Curves > demo).

Feature: public/packages/Tutorials/bdd/features/howto/16-curves/07-curves-per-compound.feature
