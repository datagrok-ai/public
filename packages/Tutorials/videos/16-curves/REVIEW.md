# Review: 16-curves (Dose-response curves)

Fresh-eyes review (2026-11-22/23, ~23:55-00:05) of the 7 delivered clips: feature, descriptions.md block, and contact
sheets at the key lines and results of each delivered take (every delivered video matched its latest take's length
before the review). Judged on: is the one idea clear to someone new to Datagrok, is each line true when it shows,
plain wording, jumps, nothing about how the video was made, and no clip under ~35 s. The m:ss in descriptions.md
match the delivered videos.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Dose-response curves in a table | OK. Assay_Curves.csv opened from Files > App Data > Curves > demo on camera (43 rows); the point's tooltip reads concentration 1.00E-7, percent_inhibition 48.813; the context panel's Chart and Fit (R² 0.977, IC50 1.63E-7, Hill 0.79 ...); the double-click's Edit chart window, closed. |
| 02 | Change how curves are fitted and drawn | OK. Level Column; linear turns every curve in the column into a straight line, 4pl brings the curves back; Confidence Interval draws the bands, IC80 droplines next to IC50, ic50 and R² printed on every plot ("R²: 0.98, IC50: 1.63E-7" on row 1); the click on row 2 shows the same options. |
| 03 | IC50 and other statistics as columns | refilmed (0:51). The line "A new column with the IC50 of every curve, next to the curves" sat over a column that reads 0.000000 in every row (molar IC50s under the default format), which looks like a failure until the scientific-format step; it now says "A new column with every curve's IC50, too small to read yet", which the next line resolves. The clip ended on the menu closing after Sort, so "the poorest fits come first" was never shown being read: a closing hover on the top row's R² (0.00, a flat curve through scattered points) now carries "R² near 0: the curve does not fit these points". The maker's note: the file's own AUC/Hill/Min/Max/IC50 columns are hidden off camera at the start (a starting state, not a jump). |
| 04 | Leave out an outlier | OK. Row 24: R² 0.93, IC50 4.96E-8; the low point (1.11E-9, 83.23) clicked becomes a red cross, R² 0.98 and IC50 4.25E-8; the second click restores 0.93 / 4.96E-8. Every number in the lines is on screen when it shows. The R² and IC50 columns are a starting state (the first line and the description say where they come from). |
| 05 | Compare curves in one chart | OK. Rows 34-37 (compound_id 20, runs 5, 4, 3, 2) selected by their row headers, MultiCurveViewer from Add viewer, 4 curves with a legend; run 3 sits left of the others (its IC50 in the file is -8.54 log M against -7.98 to -8.02 for the other three: lower, as the line says); the hover on row 38 adds compound_id 21's run in red. The viewer's title keeps "compound_id 20 ... run_id 2" (the maker's note: a product bug, the title cannot be changed). |
| 06 | Open a GraphPad Prism file | OK. dose-response-4pl.pzfx double-clicked in Browse: one row, four compounds' curves in one cell; the Edit chart window with the legend (Compound A-D, "D (with outliers)") and the excluded values as crosses; the Fit pane lists R², AUC, Top, Slope, IC50 ... per compound. |
| 07 | Curves per compound in a trellis plot | OK. Trellis plot from the toolbox, the toolbox closed on camera, Curves as the inner chart, X split by Compound ID, Reset Y columns, one cell per compound; Misc > On Click: Filter; the click on compound_id 20's cell leaves its 4 runs (status bar "Filtered: 4"). The On Click property's tooltip stays over the plot for ~3 s at the end (the maker's note, cosmetic). |

Verdicts: 6 OK (01, 02, 04, 05, 06, 07), 1 refilmed (03).

## Data claims

All true in the frames: 01's point values and fit statistics; 02's printed IC50 and R²; 03's IC50 column and the 0.00
R² on top after the sort; 04's 0.93 -> 0.98 and 4.96E-8 -> 4.25E-8 and back; 05's run 3 to the left (lower IC50,
checked in public/packages/Curves/files/demo/Assay_Curves.csv) and compound_id 21 in red; 06's crosses and per-compound
statistics; 07's 4 runs of compound_id 20.

## Titles and descriptions

Plain and specific, no hype words; each names where the file is in Browse. No card, line or description mentions how
a video was made (the playlist header's note on what runs in the browser is for maintainers); the shortest clips are
06 (0:37), 04 and 05 (0:38).

## Systematic

- A grid's "cell N of <column>" names the table row, not the row's place on screen: after a sort, the top row is
  whichever table row sorted first (03: row 7).
- A new double column of tiny values (molar IC50s) shows 0.000000 everywhere until its format is changed; a line over it
  must not promise readable values.

## For the lead

- Product: the Multi curve viewer's title stays on the first row's curve name and cannot be renamed (05, the maker's).
- Product (cosmetic): a property's description tooltip in the context panel outlives the pick (07).
- Product/design: a column added from the Fit pane's "+" gets the default number format, so IC50s in molar read
  0.000000 in every row (03).
