# Work with data

Playlist: **Work with data**. Feature name of every clip: `How-to · Data`. Features in
`public/packages/Tutorials/bdd/features/howto/05-data/`.

## Plan

One idea per clip, most important first. Selecting and filtering are another section's.

1. Add a calculated column: Add New Column, a formula over columns and functions, the preview, editing the
   formula later in the context panel, the column recalculating when a value it reads changes (demog).
2. Aggregate rows: Data > Aggregate Rows: group by, pivot, the aggregations and their functions, ADD to get
   the result as a table (demog).
3. Combine tables by joining or stacking them: Data > Join Tables on two Northwind demo files that share
   ProductID (order-details, products): the two tables, key columns, the columns to take, the join type, the
   result table; then Data > Append Tables stacks two hts plates with the same columns. Took in the former 14
   (Append tables), under 35 seconds on its own.
4. Link tables: Data > Link Tables: selecting products filters their order lines, side by side.
5. Column properties: F2 or the column menu: rename, format and units tags, change type (cars).
6. Edit data: edit a cell, add rows, remove rows and a column, undo; find and replace (Ctrl+H, filter matching
   rows, replace); Batch Edit sets one value in the selected rows (cars). Took in the former 7 (find and replace)
   and 15 (batch edit), each under 35 seconds on its own.
8. What is in the Data menu: a short tour of its items (all of them run in the browser; Curves is left out).
9. Export data: download the table, or only the filtered rows through As CSV (options).
10. and on: the transforms, one clip each, browser-only: Split Column, Unpivot, Categorize, Batch Edit,
   then ML > Cluster Data (k-means, core) and ML > Impute Missing Values if its code path checks out.
   Calculate MSR (Python) and anything server-side stay out.
16. Second round, the rest of the Data menu, one clip each: Append Tables (two plates; later merged into 3,
   "Combine tables by joining or stacking them"), Compare Tables (two
   versions of cars), Split Column by RegExp (titanic names: a pattern, not a separator; later merged into 10,
   "Split a column by a separator or a pattern"), Anonymize (demog), Compare Columns (PowerPack; tennis players' country and birthplace). Curves > Data to
   Curves is left out: it wants dose-response data with batch, assay and run columns, not a general table.

## 01 - Add a calculated column (0:56)

Title: Add a calculated column in Datagrok

Description: How to add a column computed from a formula. The Add New Column icon on the table's toolbar
opens a dialog where you name the column, pick functions from the list on the right and refer to other
columns as ${NAME}; the preview shows the first values before you add it. Later you can change the formula
in the column's Formula pane in the context panel, and the column recalculates when a value it reads changes.
Uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/01-add-a-calculated-column.feature`

## 02 - Aggregate rows (0:52)

Title: Aggregate rows: group by, pivot and summarize in Datagrok

Description: How to summarize a table group by group. Data > Aggregate Rows adds a pivot table to the view:
Group by sets the rows of the result, Pivot spreads them into one column per value, and Aggregate says what to
compute; right-click an aggregate to change its function (here, the median). The ADD button turns the
result into a table of its own. Uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/02-aggregate-rows.feature`

## 03 - Combine tables by joining or stacking them (1:30)

Title: Combine tables in Datagrok: join on a key or stack the rows

Description: Two ways to put tables together, both in the Data menu. Join Tables adds columns from one table to
another by a key they share: pick the two tables and the key columns (here ProductID), choose the columns to bring
over and the join type, and OK opens the result as a new table, here every order line with its product's name.
Append Tables stacks tables that have the same columns one below the other: tick the tables, and two 96-well
plates become one table of 192 rows. Uses the Northwind order-details and products files (Files > Demo >
northwind) and the plate-96-1 and plate-96-2 files (Files > Demo > hts).

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/03-combine-tables.feature`

## 04 - Link tables (0:56)

Title: Link two tables so one filters the other in Datagrok

Description: How to make one table follow another without joining them. Data > Link Tables picks the two
tables, the key columns that match their rows, and the link type: what travels from the first table (the
current row, the selection, the filter) and what it does in the second. Here selecting products filters the
order lines next to them to those products; Ctrl+click (Cmd+click on a Mac) adds one more. Uses the Northwind products and
order-details demo files (Files > Demo > northwind).

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/04-link-tables.feature`

## 05 - Column properties (1:06)

Title: Rename a column, set its format and units, change its type

Description: How to edit a column's properties. Press F2 on any cell (or pick Column Properties in the column
header's menu) to rename the column and edit its tags: a format such as #,##0 for thousands separators, units,
or any key and value you want to keep with the column. The tags also show in the column's Details in the
context panel. Change Type in the header's menu converts the column, here from whole numbers to decimals. Uses
the cars demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/05-column-properties.feature`

## 06 - Edit cells and rows, find and replace, edit many rows at once (1:40)

Title: Edit a table in Datagrok: cells, rows, find and replace, batch edit

Description: How to change a table by hand. Double-click a cell to edit it, use Edit > Add Rows to insert empty
rows, select rows by their row headers and press Shift+Delete to remove them, and remove a column from its
header's menu; Ctrl+Z undoes the last change and Ctrl+Y redoes it (Cmd on a Mac). Ctrl+H opens Find and Replace,
where Filter matching rows shows what will change before OK replaces it. To set one value in many rows, select
them and open Data > Batch Edit. Uses the cars demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/06-edit-a-table.feature`

## 07 - What is in the Data menu (0:39)

Title: A tour of the Data menu in Datagrok

Description: A quick walk down the Data menu of a table view, one line per item: Aggregate Rows, Unpivot,
the commands that work on two tables (Join, Link, Append and Compare Tables, available once a second table is
open), Batch Edit, Split Column, Split Column by RegExp, Anonymize, Categorize and Compare Columns. Other
videos in this playlist show several of them step by step. Uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/07-the-data-menu.feature`

## 08 - Export data (0:35)

Title: Export a table, or only its filtered rows, from Datagrok

Description: How to download a table as a file. The Export icon on the table's toolbar lists the formats (CSV,
Excel, JSON, Parquet and more); As CSV saves the whole table. To save only some rows, filter them first and
pick As CSV (options): tick Filtered rows only (or Selected rows only, or only selected columns) and the file
holds just those. Uses the demog-1000 demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/08-export-data.feature`

## 09 - Split a column by a separator or a pattern (1:03)

Title: Split a text column by a separator or a pattern in Datagrok

Description: Two ways to split one text column into several. Click a cell of the column, then Data > Split Column:
choose the separator (here a space) and demog's DEMOG ("53 C F") becomes three columns: age, race and sex. When no
single separator does, Data > Split Column by RegExp takes a regular expression, and every group in parentheses
becomes a column; here (.+), (\w+)\. (.+) splits Titanic names such as "Allen, Miss. Elisabeth Walton" into
surname, title and given names. Both dialogs preview the new columns before you add them. Uses the demog and
titanic demo files (Files > Demo).

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/09-split-a-column.feature`

## 10 - Unpivot columns into rows (1:04)

Title: Unpivot columns into rows in Datagrok

Description: How to turn several measurement columns into rows. Data > Unpivot asks what to do with each
column: Copy keeps it as it is in every new row, Merge folds it into two new columns, one with the column's
name and one with its value. Here the four iris measurements become one row per flower and measurement, 600
rows from 150. Uses the iris demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/10-unpivot.feature`

## 11 - Categorize numbers into ranges (0:36)

Title: Bin a number column into ranges in Datagrok

Description: How to turn numbers into named ranges. Data > Categorize: pick the column, click Bin... and choose
how many ranges of equal width you want; the dialog counts the rows that match. SAVE COLUMN adds a new
column with each row's range, color-coded, ready to group, color or filter by. Uses the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/11-categorize-numbers.feature`

## 12 - Cluster data (0:40)

Title: Cluster rows with k-means in Datagrok

Description: How to group similar rows. ML > Cluster > Cluster Data runs k-means right in the browser: choose
the feature columns (leave out IDs and row numbers), the number of clusters, the distance and normalization.
Show Scatter Plot colors a plot by cluster while you adjust the settings, and OK adds a clusters column. Uses
the iris demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/12-cluster-data.feature`

## 13 - Fill in missing values (0:47)

Title: Fill in missing values with k-nearest neighbors in Datagrok

Description: How to fill the empty cells of a table from the rows most like them. ML > Impute Missing Values opens
k-NN Imputation: Impute lists the columns that have gaps (here height, weight and age), Using the columns that
decide which rows are alike (leave out IDs), and Neighbors how many similar rows each value comes from. RUN fills
the gaps in place, or, with In-place unticked as here, adds a filled copy of each column next to the original so
you can compare them. Uses the demog demo dataset, which misses 751 heights, one weight and one age.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/13-impute-missing-values.feature`

## 14 - Compare two versions of a table (0:42)

Title: Compare two versions of a table: added, removed, changed rows

Description: How to see what changed between two versions of a table. Data > Compare Tables takes the old and the
new table, a key column that tells which row is which, and the columns to compare (all shared ones to start
with). OK opens a diff table with only the rows that differ: green rows are new, red rows were removed, and an
orange cell holds a changed value; point at it to see the old one. Uses two versions of the cars demo table.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/14-compare-tables.feature`

## 15 - Anonymize data (0:37)

Title: Anonymize a table before you share it in Datagrok

Description: How to hide the real values of a table. Data > Anonymize lists the columns: tick the ones to
anonymize and press OK. Text values become codes (each subject ID turns into "USUBJID 0", "USUBJID 1" and so on),
and numbers are multiplied by a random factor; the Number randomization factor sets how far, about 25% either
way at 0.1. Remember the columns that repeat a value, such as DEMOG here, which carries the age. Uses the demog
demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/15-anonymize-data.feature`

## 16 - Compare two columns (0:46)

Title: Check where two columns agree, row by row, in Datagrok

Description: How to find the rows where two columns disagree. Data > Compare Columns asks for the table (or two
tables of the same length) and the two columns, then counts the rows that match and those that don't; Select
selects the mismatches, and Add Column keeps the result as a true/false column. Here the country tennis players
represent is compared with where they were born: some differences are real, others are messy data, such as a
country written twice. Uses the tennis_players demo file (Files > Demo).

Feature: `public/packages/Tutorials/bdd/features/howto/05-data/16-compare-columns.feature`
