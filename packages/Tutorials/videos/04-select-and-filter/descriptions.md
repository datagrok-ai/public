# Select and filter

Playlist: **Select and filter**

Plan, one idea per clip, in this order:

1. **Select rows**: a click, Shift+click and Ctrl+click in the grid, a box drawn on a scatter plot, the Select menu,
   the Selected count in the status bar, Esc to clear. demog.
2. **The filter panel**: open it, a categorical and a numeric filter, every viewer follows, the Filtered count,
   one filter switched off, reset. demog.
3. **Categorical filter tricks**: keep one category, step through them with the arrow keys, search in the card,
   invert, the counts (a click on a count selects those rows). demog.
4. **Expression filter**: rules from the Column / Operation / Value form, OR and AND, a rule switched off. demog.
5. **Hierarchical filter**: levels picked in the card, a branch kept, a level opened, one child switched off. demog.
6. **Text search**: the text filter on a long-text column, several terms with OR and AND, then the table's Search
   box (Ctrl+F). beer.
7. **Substructure filter**: a structure drawn or pasted into the Structure card, the matching rows, Not contains.
   spgi (100 molecules).
8. **Selection and filter together**: Select > Missing Values, Selection to Filter, Filter to Column. demog.
9. **Filter from a viewer**: a bar chart that filters on click, an axis slider on a scatter plot that zooms and
   filters, Reset View. demog.

Added after the plan, in order of usefulness:

10. **Duplicate rows**: Select > Duplicates, Selection to Top, keep one copy, Edit > Remove > Selected Rows. iris.
11. **Filter to a list of IDs and extract the rows**: a card for the Id column, a pasted list in its search, then
    the status bar's Filtered menu, Extract Rows into a new table. spgi.
12. **Saved filters**: Save or Apply > Save..., reset, apply the preset again. demog.
14. **Select and filter from a column**: Shift+Enter on a cell, Shift+Enter in the Search box, the Column options
    popup's filter, Add filter to the panel. demog.

(13 and 15 were merged into 11 and 14: on their own they were under 30 seconds.)

## 01 - Select rows (0:41)

Title: How to select rows in Datagrok

Description: The ways to select rows in a table: a click makes a row current, Shift+click selects a range,
Ctrl+click (Cmd+click on a Mac) adds or removes one row, and Shift+drag on a scatter plot selects the points
inside the box. The Select menu has All, None, Invert and more, the status bar shows how many rows are selected,
and Esc clears the selection. Selected rows are highlighted in the grid and in every viewer. Uses the demog demo
dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/01-select-rows.feature

## 02 - Filter rows with the filter panel (0:46)

Title: Filter rows with the filter panel in Datagrok

Description: The filter icon in the toolbar opens the filter panel, with a card for each column. Click a category
to keep only its rows, or drag the ends of a numeric range to narrow it; the grid and every viewer show only the
rows that pass, and the status bar counts them. A card's checkbox turns that filter off, and the reset icon at
the top of the panel clears all filters. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/02-the-filter-panel.feature

## 03 - Categorical filter tricks (0:39)

Title: Categorical filter tricks in Datagrok

Description: A few things a categorical filter card does besides ticking boxes. Click a category name to keep only
that category, then use the arrow keys to step to the next one; a checkbox adds a category back. Invert all, in
the menu next to the card's name, swaps what is kept and what is hidden, the search icon finds categories by
name, and a click on a count selects those rows. Uses the DIS_POP column of the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/03-categorical-filter-tricks.feature

## 04 - Filter with an expression (1:03)

Title: Filter with an expression in Datagrok

Description: The expression filter builds rules such as AGE > 50 from a column, an operation and a value. Add it
from the filter panel's menu (Add Filter > Expression), then add one rule at a time with the plus button. Rules
combine with OR until you click OR to switch to AND, and each rule has its own checkbox to turn it off. Uses the
demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/04-expression-filter.feature

## 05 - Filter with a hierarchy (0:43)

Title: Hierarchical filter in Datagrok

Description: The hierarchical filter shows several columns as levels of one tree, here SEX and then RACE. Add it
from the filter panel's menu (Add Filter > Hierarchical) and pick the columns for each level with the tree icon
on the card. Click a branch to keep only its rows, open it with the triangle to see the next level and its
counts, and untick a branch to leave it out. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/05-hierarchical-filter.feature

## 06 - Search text in a table (0:38)

Title: Search text in a table in Datagrok

Description: A column of long text gets a search box in the filter panel: type a word, press Enter, and only
the rows that mention it stay, with the word highlighted in the grid. Several words combine with OR, and a
click switches to AND. For the whole table, press Ctrl+F on the grid and type a condition such as
Origin = Belgium in the Search box. Uses the beer demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/06-text-search.feature

## 07 - Filter molecules by a substructure (0:36)

Title: Filter molecules by a substructure in Datagrok

Description: A molecule column gets a structure card in the filter panel. Click Sketch, draw a fragment or type or
paste its SMILES, and only the molecules that contain it stay, with the matching atoms highlighted in the grid.
The search type under the card switches to Not contains, Exact, Similar and more. Uses the spgi demo dataset
(100 molecules).

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/07-substructure-filter.feature

## 08 - Turn a selection into a filter (0:42)

Title: Turn a selection into a filter in Datagrok

Description: Rows can move between the selection, the filter and a column. Select > Missing Values selects the
rows with an empty cell in the columns you tick, Select > Selection to Filter shows only the selected rows, and
Select > Filter to Column keeps the current filter as a true/false column you can use later. Uses the demog demo
dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/08-selection-to-filter.feature

## 09 - Filter from a viewer (0:46)

Title: Filter from a viewer in Datagrok

Description: Viewers can filter the table too. Right-click a bar chart and pick On Click > Filter: a click on a bar
keeps only its rows, a click on another bar switches to it, and a click on empty space lets go. On a scatter plot,
zooming in (here by dragging the end of the axis slider) filters the table to the points still in view, and
Reset View in the right-click menu brings every row back. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/09-filter-from-a-viewer.feature

## 10 - Find and remove duplicate rows (0:48)

Title: Find and remove duplicate rows in Datagrok

Description: Select > Duplicates selects the rows whose values repeat another row's in the columns you tick; the
dialog also shows how many unique values each column has. Select > Selection to Top sorts the grid so the
selected rows come first, so you can compare them. Ctrl+click (Cmd+click on a Mac) one copy to keep it, then
Edit > Remove > Selected Rows deletes the rest. Uses the iris demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/10-remove-duplicate-rows.feature

## 11 - Filter to a list of IDs and copy the rows into a new table (0:46)

Title: Filter to a list of IDs and copy the rows to a new table in Datagrok

Description: How to keep only the rows on a list of compound IDs, such as one a colleague sends you, and take them
out as a table of their own. The + at the top of the filter panel adds a filter for any column, here Id; open the
card's search with the magnifier and paste the list, one ID per line as it comes from a spreadsheet or an email.
Then click the Filtered count in the status bar: its menu can keep or delete the rows, and Extract Rows copies them
into a new table, leaving the original as it was. Uses the spgi demo dataset (100 molecules).

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/11-filter-to-a-list-of-ids-and-copy-the-rows.feature

## 12 - Save filters and apply them later (0:43)

Title: Save filters and apply them later in Datagrok

Description: Once the filter panel is set the way you need, open the panel's menu and pick Save or Apply > Save...
to keep its filters under a name. Reset clears them, and the name under Save or Apply puts them all back in one
click. Saved presets are kept on the server, so they are there the next time you open the table. Uses the demog
demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/12-save-and-reuse-filters.feature

## 13 - Select and filter from a column (0:45)

Title: Select and filter from a column in Datagrok

Description: Quick ways to select and filter by one column, right in the grid. Click a cell and press Shift+Enter
to select every row with the same value; Esc clears the selection. In the table's Search box (Ctrl+F), type a
condition such as AGE > 70: Shift+Enter selects the matching rows, Enter would filter them. Point at a column
header and click the small icon that appears to filter that column from its histogram; Add filter moves the
filter to the filter panel. Uses the demog demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/04-select-and-filter/13-select-and-filter-from-a-column.feature

