# Summary columns in the grid

These scenarios check the summary columns that PowerGrid adds to the grid. A summary column is an
extra grid column that draws several values of the same row as one small picture: a line, bars, a
pie, a radar, a compact form, tags or an interval. It belongs to the grid only, so the table keeps
its own columns. PowerGrid offers eight such items in the grid's **Add > Summary Columns** menu:
Sparklines, Bar Chart, Pie Chart, VlaaiVis, Radar, Smart Form, Tags and Confidence Interval. The
three form items above them in the same menu belong to the grid itself and are not covered here.

## Setup

Run these steps before each scenario:

1. Close all views.
2. Open the demog-1000 dataset (`System:DemoFiles/demog-1000.csv`) and wait for the grid to load.
3. The grid shows the columns USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL,
   STARTED and SEVERITY, and the status bar says "Columns: 11" and "Rows: 1,000".

By default, every summary item except Tags and Confidence Interval draws the AGE, HEIGHT and WEIGHT
values of its row. Tags shows the CONTROL column: a CONTROL tag in every row whose CONTROL box is
checked. Confidence Interval draws AGE as a dot, with an interval that spans from HEIGHT to WEIGHT.

The steps below say what the Tags column should draw. On the current build it paints the whole grid
over instead, so every Tags check fails — that is GROK-20888, and these steps pass again once it is
fixed.

## Scenario 1: Each Summary Columns item adds its own summary column

1. Right-click the USUBJID cell of row 2, point at **Add**, then at **Summary Columns**.
2. The submenu lists `<item>` among the eight items that follow **Custom HTML Form...**.
3. Click `<item>`.
4. A new column named `<item>` appears at the right end of the grid, after SEVERITY.
5. The status bar still says "Columns: 11".
6. The cells of the `<item>` column show `<drawing>`.
7. Right-click the `<item>` column header and pick **Remove**.
8. The `<item>` column disappears, and the grid shows only the eleven table columns again.
9. No errors appear in the console.

| item | drawing |
|---|---|
| Sparklines | a small line through three dots, one dot each for AGE, HEIGHT and WEIGHT |
| Bar Chart | three small bars side by side, for AGE, HEIGHT and WEIGHT |
| Pie Chart | a pie of three colored sectors for AGE, HEIGHT and WEIGHT, with sector sizes that differ from row to row |
| VlaaiVis | a round pie of three colored sectors inside a faint outer ring |
| Radar | a small three-axis radar: a filled shape with a dot on each of the three axes |
| Smart Form | the labels AGE, HEIGHT and WEIGHT with the row's values; in row 2 these are AGE 30, HEIGHT 150.29 and WEIGHT 64.00 |
| Tags | a CONTROL tag in each row whose CONTROL box is checked (rows 3, 8, 11 and 19 near the top) and nothing in the other rows, while the other columns keep showing their values |
| Confidence Interval | a dot and a dashed horizontal line with a short vertical mark at each end |

## Scenario 2: The Remove selected columns icon removes a selected summary column (GROK-18256)

In GROK-18256, this icon did nothing for a summary column, while the column's own **Remove** menu
item worked.

1. Right-click the USUBJID cell of row 2 and pick **Add > Summary Columns > Sparklines**.
2. A Sparklines column appears after SEVERITY.
3. On the toolbar above the grid, the icon with the tooltip "Remove selected columns" is greyed out.
4. Hold Ctrl and click the Sparklines column header.
5. The Sparklines column is highlighted as selected, and the Remove selected columns icon becomes
   active.
6. Click the Remove selected columns icon.
7. The Sparklines column disappears, the grid shows only the eleven table columns, and the status bar
   still says "Columns: 11".
8. No errors appear in the console.

## Scenario 3: A Tags column keeps working after its source column is removed (GROK-19942)

In GROK-19942, the whole grid stopped drawing after the column that a Tags cell refers to was
removed.

1. Right-click the USUBJID cell of row 2 and pick **Add > Summary Columns > Tags**.
2. A Tags column appears after SEVERITY. It shows a CONTROL tag in each row whose CONTROL box is
   checked (rows 3, 8, 11 and 19 near the top).
3. Right-click the CONTROL column header and pick **Remove**.
4. CONTROL disappears from the grid, and the status bar says "Columns: 10". The Tags column stays in
   place, and its cells are empty.
5. Point at the grid and turn the mouse wheel down for several screens.
6. The grid keeps drawing as it scrolls: the row numbers move on to rows past 100, every column shows
   the values of those rows, and no error balloon appears.
7. Click any AGE cell and press Ctrl+Z.
8. CONTROL comes back between DEMOG and STARTED, and the status bar says "Columns: 11". The Tags
   column again shows a CONTROL tag in each row whose CONTROL box is checked.
9. Right-click the Tags column header and pick **Remove**. The Tags column disappears.
10. No errors appear in the console.

## Scenario 4: Summary columns come back from a saved layout and from a saved project (GROK-19769)

In GROK-19769, a layout that held summary columns failed to apply.

1. Add the eight summary columns one after another, in this order: Sparklines, Bar Chart, Pie Chart,
   VlaaiVis, Radar, Smart Form, Tags, Confidence Interval. For each one, right-click a cell of a
   table column (for example, the AGE cell of row 2) and pick **Add > Summary Columns** and the item.
2. Scroll the grid to the right. After SEVERITY, the grid shows the eight summary columns in the
   order they were added, and each one is drawn as described in Scenario 1.
3. In the **Toolbox**, expand **Layouts** and click **SAVE**.
4. Type demog-1000 in the **Filter by name or #tag** box of the Layouts pane. A layout card named
   demog-1000 is listed.
5. Right-click the dark sidebar on the left edge of the window and pick **Close All**. All views
   close.
6. Open demog-1000 again. The grid shows only the eleven table columns.
7. Add a Scatter plot from **Toolbox > Viewers**. The view shows the grid and a scatter plot.
8. In **Toolbox > Layouts**, type demog-1000 in the filter box and click the demog-1000 card.
9. The scatter plot is gone. After SEVERITY, the grid again shows the eight summary columns in the
   same order, each drawn as described in Scenario 1, and no error balloon appears.
10. Click **SAVE** on the toolbar above the grid. In the **Save project** dialog, replace the name
    with PgSummaryColumns and click **OK**.
11. A balloon says "Project "PgSummaryColumns" uploaded.", and a **Share PgSummaryColumns** dialog
    opens. Click **CANCEL**.
12. Right-click the dark sidebar and pick **Close All**.
13. Open **Browse > Dashboards**, type PgSummaryColumns in the search box, and double-click the
    PgSummaryColumns card.
14. The project opens the demog-1000 table. After SEVERITY, the grid shows the eight summary columns
    in the same order, each drawn as described in Scenario 1, and no error balloon appears.
15. In **Toolbox > Layouts**, type demog-1000 in the filter box, right-click the demog-1000 card,
    pick **Delete**, and confirm with **DELETE**. The card disappears.
16. Right-click the dark sidebar and pick **Close All**.
17. Open **Browse > Dashboards**, type PgSummaryColumns in the search box, right-click the
    PgSummaryColumns card, pick **Delete Project**, and confirm with **DELETE**. The card disappears.
18. No errors appear in the console.
