# Datagrok basics

Playlist: **Datagrok basics**. Features: `public/packages/Tutorials/bdd/features/howto/01-getting-started/`
(feature name `How-to · Basics`).

## Plan

| NN | Clip | One idea | Data |
|---|---|---|---|
| 01 | Find your way around Datagrok | The long tour: Home and its search, the sidebar, Browse's sections opened to show what is in them, a table opened from Files > Demo, its menus, icons, toolbox, context panel, status bar and the toggles, and Chem / Bio menus that come with molecules and sequences (every move filmed) | demog, smiles.csv, peptides.csv (all from Browse), the Biologics connection |
| 02 | Open your data | Four ways in: Browse (Files > Demo: click previews, double-click opens), the open-file icon, a dropped file, pasted text (Open text) | cars, two small local CSVs |
| 03 | Work with the grid | Move around, sort, resize, reorder, hide and show columns, row height; in a wide table, pin a column and find another by name in the Columns pane | demog, spgi |
| 04 | Format and color-code columns | Number and date formats; linear and categorical color coding | demog |
| 05 | The context panel | What it shows for the table, a column and a molecule cell; its panes | demog, smiles.csv (from Browse) |
| 06 | Arrange your workspace | Tabs switch tables, a tab dragged to a side docks two views, the splitter shares the width; the status bar toggles: toolbox, context panel, console, tabs, presentation mode | cars, iris |
| 07 | Get help and learn | Context Help on what you click, the Help menu of the sidebar, and Help > Tutorials: a tutorial whose steps tick off as you go | demog (and the tutorial's) |
| 08 | Search everything from Home | The Home search: functions, help, dashboards, files and queries; a sentence a query answers; IDs looked up as cards and opened as a table | Biologics dashboards, queries and IDs |

After the lead's 23:00 round ("17 seconds is too short"), the short clips were merged into the ones above: 09 "Work
with several tables side by side" into 06, 10 "Find a column in a wide table" into 03, 11 "Learn with interactive
tutorials" into 07, and 12 "Search everything from one box", 13 "Ask a question in the search box" and 14 "Look up
an ID" into 08.

## 01 - Find your way around Datagrok (3:34)

Title: Find your way around Datagrok

Description: A tour of the Datagrok window. It starts on Home with its search box and the sidebar, then opens each
section of Browse to show what is in it: My stuff, Spaces, the installed apps, the Dashboards gallery, a database
connection with its queries and tables, Platform with its plugins, and the demo files. Then demog.csv opens from
Files > Demo, and the tour goes round what surrounds a table: the menus and icons above the grid, the toolbox on the
left, the context panel on the right, and the status bar with the toggles that show and hide panels. At the end, a
table of molecules (smiles.csv) adds a Chem menu and a table of peptide sequences (peptides.csv) adds a Bio menu.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/01-find-your-way-around.feature`

## 02 - Open your data (0:40)

Title: Open your data in Datagrok: Browse, your computer, drop, paste

Description: Four ways to get a table into Datagrok. In Browse > Files > Demo, one click on cars.csv previews it and
a double-click opens it. The open-file icon at the top of Browse opens a file from your computer, and you can also
drag a file onto the window. Copied cells or CSV text go into the Open text view, which shows a preview and opens it
as a table when you click DONE.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/02-open-your-data.feature`

## 03 - Work with the grid (1:28)

Title: Work with the Datagrok grid: sort, resize, move, hide, find columns

Description: The everyday moves in Datagrok's grid, shown on the demog demo table. The arrow keys move the current
cell, Ctrl+End and Ctrl+Home jump to the last and first rows. A double-click on a column header sorts
by it, dragging a header's edge resizes the column, and dragging the header moves it. Right-click a header to hide
a column, and use Order or Hide Columns from the same menu to bring it back; dragging between row numbers changes
the row height. For a table wider than the screen (spgi, 100 molecules and 88 columns), Pin > Pin Column keeps a
column in sight, and "Columns" in the status bar lists them all: type to find one, click it, and the grid takes you
there.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/03-work-with-the-grid.feature`

## 04 - Format and color-code columns (0:36)

Title: Format and color-code columns in Datagrok

Description: Right-click a column header in the grid to change how its values look. Under Format, every number
format is listed with a value written in it, so you can pick by eye; a date column lists date formats instead.
Color Coding > Linear shades a numeric column from low to high, and Color Coding > Categorical gives each category
its own color. Shown on the demog demo table.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/04-format-and-color-code-columns.feature`

## 05 - The context panel (0:53)

Title: The context panel in Datagrok: details of what you click

Description: The panel on the right side of Datagrok shows whatever you clicked last, split into panes that open
with a click. Clicking the table's name in the status bar shows the table (rows, columns and more), clicking a column
header shows the column with its stats and a plot of its values, and clicking a molecule shows properties computed
for that molecule. Shown on the demog demo table and on smiles.csv, opened from Browse > Files > Demo > chem.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/05-the-context-panel.feature`

## 06 - Arrange your workspace (0:53)

Title: Arrange your Datagrok workspace: side-by-side tables and panels

Description: Every table you open gets its own tab at the top, and a click on a tab switches to it. Drag a tab to
the side of another view to see both tables at once, then drag the line between them to share the width. The small
toggles at the right end of the status bar show and hide the panels around your data: the toolbox, the context
panel (F4), the console and the view tabs. Presentation mode (F7) shows just the current view, and the "back to
design mode" link at the top right brings the panels back. Shown with the cars and iris demo tables.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/06-arrange-your-workspace.feature`

## 07 - Get help and learn (1:16)

Title: Get help in Datagrok and learn it with interactive tutorials

Description: Turn on Context Help with its toggle at the bottom right (or F1), and it shows the help page of
whatever you click: a column brings the page on tables and columns, a scatter plot the scatter plot's page; the
same toggle closes it. The Help icon in the sidebar opens a menu with the documentation, the community forum, the
YouTube channel and a feedback form. To learn by doing, pick Help > Tutorials in the menu behind the logo: a click
on a card starts a tutorial, whose bar on top tells you what to do next and ticks each step off as you do it.
Shown on the demog demo table.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/07-get-help-and-learn.feature`

## 08 - Search everything from Home (1:09)

Title: Search everything from the Datagrok Home page

Description: The search box on the Home page looks through the whole server. "scatter plot" lists functions, demos,
help pages and code samples, each kind under its own heading; a word from a name finds dashboards, files and
queries, and a double-click opens a dashboard (the house tab takes you back to Home). Some queries answer a
sentence: start typing and the box suggests them, finish one, such as "biologics assays for E. coli", and the query
runs and shows its table right under the box. An ID such as GROKADC-000001 shows its record as a card; several IDs
separated by commas show a card each, and "Open table" puts them in a table whose IDs stay links. The dashboards,
queries and IDs here come with the Biologics plugin.

Feature: `public/packages/Tutorials/bdd/features/howto/01-getting-started/08-search-everything.feature`

