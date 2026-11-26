# Review: 05-data (Work with data)

Fresh-eyes review (2026-11-22, ~21:30-21:45): feature lines, descriptions.md, and a contact sheet at the key lines and
results of each clip. Every delivered video matched its latest guide take (same length), so the frames come from the
takes. The m:ss in the descriptions.md headings were the older renders: all refreshed to the delivered lengths.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Add a calculated column | OK. The formula editor's autocomplete pops up briefly while the formula is typed (harmless). |
| 02 | Aggregate rows | OK. |
| 03 | Join two tables | text fixed: the description said the dialog "tells you how many rows match" and that we set ProductID "on both sides"; it only says whether the key rows match, and we pick ProductID on the order-details side. The dialog opens with "None of key rows are matched" in red until ProductID is picked: the caption covers it. Result table named "result" (platform default). |
| 04 | Link tables | text fixed: Mac Cmd+click added. Product: once linked, the products view's top-right label reads "94 filtered rows" / "104 filtered rows" (the order-details count) although products is not filtered; while the dialog is open it reads "Dialog: Link Tables". |
| 05 | Column properties | OK. The column's tag list shows the internal tag `.semantic-detection-duration` (product bug 8 of the maker). |
| 06 | Edit data | OK. After Ctrl+Z the grid shifts left and clips the first letter of model for ~4 s (product bug 10 of the maker). |
| 07 | Find and replace | OK. With Filter matching rows on, the dialog covers the model column of the 4 filtered rows, so the viewer sees 4 rows but not their "toyota" names. The Variables pane that pops up after OK (product bug 1) is closed off camera (the scan's no-jumps warning); checked frame by frame at 27.6-28.5 s, it never shows, OK dissolves straight into the result. Left as is: closing it on camera would film a product bug into the tutorial and make the feature fail once the bug is fixed. A viewer repeating the steps today will see the pane open; the lead may prefer fixing the bug (shell.dart ~1227). |
| 08 | What is in the Data menu | OK. |
| 09 | Export data | OK. Nothing on screen confirms a download (headless browser); the line "the file goes to your downloads" covers it. |
| 10 | Split a column | text fixed: the menu caption "Data > Split Column" (a menu path, no meaning) is now "Split Column: one text column into several" (feature and take, re-rendered, re-delivered). The dialog covers part of its own preview window (product bug 9). |
| 11 | Unpivot columns into rows | OK. Slow middle: the four Merge picks take 25 s of a 64 s clip, each a full drop-down gesture. |
| 12 | Categorize numbers into ranges | OK. The Matches column of the dialog stays blank after Bin... (product bug 5). |
| 13 | Cluster data | OK. With 2 clusters, two versicolor points land in the setosa cluster (orange dots inside the blue group on the final plot); the line "setosa forms a group of its own" is still true for the group, but a sharp viewer sees the strays (k-means init, product bug 4). |
| 14 | Fill in missing values | OK. |
| 15 | Batch edit selected rows | OK. |
| 16 | Append tables | OK. Result named "result"; the dialog's labels greyed (product bug 3 of the round-2 maker). |
| 17 | Compare two versions of a table | OK. Top-right label reads "Dialog: Compare Tables" after the dialog closed (known). The old value's tooltip is small but readable at 1080p. |
| 18 | Split a column by a pattern | OK. New columns named r, r (2), r (3) (known product bug). |
| 19 | Anonymize data | OK. |
| 20 | Compare two columns | text fixed: the caption "Compare Columns: do two columns agree, row by row?" was a question; now "Compare Columns: where two columns agree, row by row" (feature and take, re-rendered, re-delivered). |

Verdicts: 16 OK, 4 text fixed (2 of them re-rendered from the same take with the new line), 0 refilmed, 0 need work.

## Titles and descriptions

Plain, accurate and specific; each names the menu path and the demo data. No hype words, no questions. Mac shortcuts
are given where a modifier key is part of the gesture (06, now 04).

## Jumps

None. Every clip's off-camera `Given` steps build the starting state before the first visible step; 07's mid-clip
`Given the variables pane is closed` hides a pane a product bug opens and never shows on screen. New views (join,
aggregation, unpivot, append, diff results) open on a filmed OK/ADD press.

## Systematic

- Several dialogs open over the grid and cover what the line talks about (07 filtered rows, 10 preview, 14 HEIGHT
  rows 400-411). A step to move a dialog aside, or a dialog that opens away from the current cell, would help.
- Result tables of Join/Unpivot/Append are named "result": fine for a tutorial, but a viewer may read it as a test
  name. Platform default.
- A caption can be changed after filming without a refilm: edit `caption` in the take's `steps.json` (and the
  feature), then `python3 public/libraries/bdd/tool/guide-render.py <take> --tutorial --quiet` and `deliver.sh`.

## Merges (2026-11-22, ~23:15-23:40)

The lead found clips under ~35 s too thin. Three merges and one extension, filmed again, reviewed with the
checklist (pass, audit 0 off, timeline, frames at every note and result) and delivered under the kept number:

| NN | Clip | Took in | Length | Notes |
|---|---|---|---|---|
| 06 | Edit cells and rows, find and replace, edit many rows at once | 07 (Find and replace), 15 (Batch edit) | 1:40 | cars throughout: a price edited, two rows added and removed, symbol removed, Ctrl+Z, Ctrl+Y; Ctrl+H on toyota1, All Columns, Filter matching rows (4), OK; rows 22-25 selected by their headers, Batch Edit puts 1 into turbo for the four Toyotas (row 26 keeps 0). First take undid the removal of turbo and kept it: a column Ctrl+Z brings back returns selected, tinted over its whole height for the rest of the clip; now the removal is redone (Ctrl+Y) on symbol. The Variables pane that Find and Replace opens is still closed off camera (never on screen). Title without a colon so the file name stays plain. |
| 10 | Split a column by a separator or a pattern | 18 (Split Column by RegExp) | 1:03 | demog's DEMOG split at a space, then a filmed click on the titanic tab (both tables opened off camera at the start) and Split Column by RegExp on name. The grid scrolls a little sideways after each split (USUBJID / pclass clipped at the left edge): the known product quirk. |
| 14 | Fill in missing values | (extended) | 0:47 | demog instead of demog-1000: 751 heights, one weight, one age missing (the Using list's nulls column shows exactly these). Impute's list (AGE, HEIGHT, WEIGHT), USUBJID taken out of Using (10 to 9), In-place unticked: AGE(copy), HEIGHT(copy), WEIGHT(copy) appear at the right, rows 2221-2223 show the gaps next to their filled values (HEIGHT 159.990 and 163.043, WEIGHT 79.80). The library's column-list toggle lights the whole Select columns dialog (known library gap). |

Data claims: every count in a line was computed from data/demo (cars.csv, demog.csv) and is asserted by the feature.

## Merge after the renumbering (2026-11-23, ~00:00)

| NN | Clip | Took in | Length | Notes |
|---|---|---|---|---|
| 03 | Combine tables by joining or stacking them | 14 (Append tables, 0:31) | 1:30 | The four files (plates first, then order-details, products) open off camera before the first step; every view change after it is a filmed tab click or an OK. With the plates open, Join Tables proposes them (the pair with the most shared columns, `findOptimalTables`), so the clip picks order-details and products in the dialog ("With more tables open, pick the two to join"): an honest extra step a person meets too. The old take's "Order lines refer to a product only by its ID" was read over the products view as the pointer went to the order-details tab; now a note after the tab click, over order-details. Append from plate-96-1: the dialog lists all five open tables, the two plates ticked, result (2) of 192 rows, rows 96/97 at the plate boundary. Kept at 1:30 by folding the plates' intro into the Append menu caption (the renderer puts that caption at the top, over part of the grid's header row, for ~3 s while the menu is open). Results are named "result" and "result (2)" (platform default). Data checked: 2155 order lines all matched, ProductName never empty; both plates 96 rows, the same 8 columns, one barcode each (1737398, 5899277). |

