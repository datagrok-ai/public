# From data to answers

Playlist: **From data to answers**

## Plan

Each clip is one realistic task from start to finish (1 to 2.5 minutes). The notes on screen say why each stage
is there; the single-idea playlists (Basics, Viewers, More viewers, Select and filter, Data, Dashboards and
sharing, Molecules) explain the controls themselves. After the starting state, every move is on camera: panels are
opened from the status bar or the sidebar, views are reached by clicks, and Browse is expanded to show where things
are.

1. From a file to a shared dashboard: demog.csv from Browse > Files > Demo, a look at a column in the Context
   Panel, three viewers from the Toolbox, the filter panel narrowed to one patient group, SAVE as a dashboard,
   share it, find it in Browse > Dashboards.
2. Explore a new dataset: cars, the status bar and the column statistics in the Context Panel, a histogram to spot
   the outliers and select them, a correlation plot to find the column that goes with price, a scatter plot of the
   pair, Selection to Filter. (No box plot: cars has no text column to group by.)
3. Compare groups: demog, the filter panel to choose who is compared, a box plot by sex with its p-value, the same
   box plot split by race and sex, then the arrangement saved as a layout. (A trellis plot of box plots was tried:
   its cells scale apart and the comparison did not read.)
4. Clean and reshape: demog, Select > Missing Values and Remove Selected Rows, a calculated column (BMI), Aggregate
   Rows into a summary table by race, a chart of the summary.
5. Combine two tables: Northwind order lines and products, Join Tables to get product names on the order lines,
   then a bar chart to answer "which products sell the most?"
6. Molecules from structure to shortlist: spgi, a substructure filter, a scatter plot of two properties colored by
   activity, select the best, extract them to a new table.

Everything a clip saves (dashboards, layouts, shares) is deleted at the end of its feature and swept at its start.

## 01 - From a file to a shared dashboard (2:13)

Title: From a CSV file to a shared dashboard in Datagrok

Description: One task from start to finish. We open demog.csv from Browse > Files > Demo, check a column's numbers in
the Context Panel, add a bar chart, a scatter plot and a histogram from the Toolbox, and use the filter panel to
keep one patient group. Then SAVE turns it into a dashboard, the Share dialog that follows gives a colleague access,
and Browse > Dashboards shows where to find it later. The Viewers, Select and filter, and Dashboards and sharing
playlists explain each of these steps in detail.

Feature: public/packages/Tutorials/bdd/features/howto/08-end-to-end/01-from-a-file-to-a-shared-dashboard.feature

## 02 - Explore a new dataset (1:12)

Title: A first look at a new dataset in Datagrok

Description: A first look at a table you have not seen before, using the cars demo dataset (8 of its columns). The
status bar and the Stats section of the Context Panel give the size and the numbers, a histogram of price shows three
cars far above the rest, and a Shift+drag over those bars selects them. A correlation plot shows that engine size goes
with price most, a scatter plot of the two puts the selected cars at the top right, and Select > Selection to Filter
shows only them.

Feature: public/packages/Tutorials/bdd/features/howto/08-end-to-end/02-explore-a-new-dataset.feature

## 03 - Compare groups (1:16)

Title: Compare two groups in Datagrok with a filter, a box plot and a layout

Description: One question from start to finish: whether men and women differ in weight, and whether that holds in
every race. We keep one patient group with the filter panel, compare weight by sex in a box plot and read its p-value, then
split the box plot by race and sex to check each race. Finally Toolbox > Layouts > Save keeps the arrangement under
a name, so it can be applied to any table with the same columns. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/08-end-to-end/03-compare-groups.feature

## 04 - Clean and reshape a table (1:38)

Title: Clean a table, add a column, summarize and chart it in Datagrok

Description: From raw rows to a small summary table and a chart. We find the rows with no height through Select >
Missing Values and remove them with Edit > Remove > Selected Rows, add a BMI column with a formula, then group by race
with Data > Aggregate Rows and turn the result into a table of its own. A bar chart of that table shows the average
BMI per race. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/08-end-to-end/04-clean-and-reshape.feature

## 05 - Combine two tables to answer a question (1:27)

Title: Join two tables in Datagrok, then answer a question with a chart

Description: The task is to find the products that sell the most. The quantities are in the order lines and the
product names in a second table, so we join them with Data > Join Tables on ProductID, bringing over only the name. A bar chart of the summed
quantity per product then gives the answer. Uses the Northwind order-details and products demo files from
Browse > Files > Demo > northwind.

Feature: public/packages/Tutorials/bdd/features/howto/08-end-to-end/05-combine-two-tables.feature

## 06 - From molecules to a shortlist (1:25)

Title: From 100 molecules to a shortlist of potent compounds in Datagrok

Description: A small medicinal chemistry task from start to finish, using the spgi demo dataset of 100 compounds.
Chem > Search > Substructure Search keeps the molecules with a pyrrolidine ring, a scatter plot of logP against
TPSA colored by activity shows them, and a Shift+drag over the first bars of an activity histogram selects the most
potent. Select > Extract Selected Rows then puts the shortlist into a table of its own, with structures and every
property. The Molecules playlist explains each Chem tool on its own.

Feature: public/packages/Tutorials/bdd/features/howto/08-end-to-end/06-molecules-to-a-shortlist.feature
