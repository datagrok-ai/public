# Correlation plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Correlation plot

## Double-click and cell interaction

1. Double-click an non-diagonal cell (e.g., AGE vs HEIGHT) -- a new scatter plot viewer opens
2. Close the scatter plot
3. On the Context Panel enable **Ignore Double Click**
4. Double-click the same cell -- no scatter plot should open
5. On the Context Panel disable **Ignore Double Click**
6. Click an non-diagonal cell -- context panel updates to show correlation details for the selected pair

## Column reordering

1. Click the AGE column header and drag it to the right -- column order updates
2. Right-click > **Columns** > **X Columns** -- reorder columns via the picker
3. Right-click > **Columns** > **X Columns** -- uncheck a column to hide it
4. Verify the viewer updates to reflect the new order and visibility

## Column selection

1. Right-click > **Columns** > **X Columns** -- deselect WEIGHT
2. Verify WEIGHT column disappears from X axis
3. Right-click > **Columns** > **X Columns** -- re-select WEIGHT
4. Right-click > **Columns** > **Y Columns** -- deselect HEIGHT
5. Verify HEIGHT row disappears from Y axis
6. Right-click > **Columns** > **Y Columns** -- re-select HEIGHT

## Correlation type and display options

1. On the Context Panel verify **Correlation Type** default is Pearson
2. On the Context Panel set **Correlation Type** to Spearman -- cell values update
3. On the Context Panel set **Correlation Type** back to Pearson
4. On the Context Panel verify **Show Pearson R** is enabled by default
5. On the Context Panel disable **Show Pearson R** -- values disappear from cells
6. On the Context Panel enable **Show Pearson R**
7. On the Context Panel verify **Show Tooltip** is enabled by default
8. On the Context Panel disable **Show Tooltip**
9. On the Context Panel enable **Show Tooltip**
10. On the Context Panel verify **Ignore Double Click** is disabled by default
11. On the Context Panel enable **Ignore Double Click**
12. On the Context Panel disable **Ignore Double Click**

## Row source

1. Select several rows in the grid
2. On the Context Panel verify **Row Source** default is Filtered
3. On the Context Panel set **Row Source** to Selected -- plot updates to show only selected rows
4. On the Context Panel set **Row Source** to All
5. On the Context Panel set **Row Source** back to Filtered

## Style customization

1. On the Context Panel > **Style** set **Default Cell Font** to a larger size
2. Verify cell height increases to match the new font size
3. On the Context Panel > **Style** set **Col Header Font** to bold larger size
4. On the Context Panel set **Back Color** to a light gray

## Title and description

1. On the Context Panel > **Description** enable **Show Title**
2. Set **Title** to "Correlation Analysis"
3. Set **Description** to "Shows pairwise correlations"
4. Set **Description Visibility Mode** to Always
5. Change **Description Position** to Bottom
6. Set **Description Visibility Mode** to Never

## Context menu

1. Right-click an off-diagonal cell -- menu should open
2. Verify **Show Pearson R** toggle item is present
3. Verify **Tooltip** submenu contains **Visible** and **Properties...**
4. Verify **Columns** submenu contains **X Columns** and **Y Columns** pickers
5. Close the menu
6. Right-click a diagonal cell (histogram cell) -- verify menu still opens

## Open as Table

1. Right-click an non-diagonal cell and select **Open as Table**
2. Verify a new table view opens with row/column pairs and correlation values

## Viewer filter formula

1. On the Context Panel > **Data** open **Filter**
2. Set filter formula to `${AGE} > 40` -- correlation plot recalculates on filtered rows
3. Clear filter formula

## Layout persistence

1. On the Context Panel set **Correlation Type** to Spearman
2. On the Context Panel disable **Show Pearson R**
3. Save layout via JS API
4. Close the Correlation plot viewer
5. Load saved layout via JS API -- verify Correlation Type is Spearman, Show Pearson R is disabled
6. Delete saved layout


---
{
  "order": 10,
  "datasets": ["System:DemoFiles/demog.csv"]
}
